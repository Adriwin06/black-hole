// Role: Core ray-marching function. Integrates a photon path through the curved
//       spacetime, compositing accretion disk, jet, planet, and background
//       contributions along the geodesic with Beer-Lambert transmittance.
//       Takes an initial ray direction and returns the accumulated HDR colour.

vec4 trace_ray(vec3 ray) {
    vec3 pos = cam_pos;

    {{#aberration}}
    vec3 aberration_vel = cam_vel * max(look_aberration_strength, 0.0);
    float aberration_speed = length(aberration_vel);
    if (aberration_speed > 0.999) {
        aberration_vel *= 0.999 / aberration_speed;
    }
    ray = lorentz_velocity_transformation(ray, aberration_vel);
    {{/aberration}}

    float ray_intensity = 1.0;
    float ray_doppler_factor = 1.0;
    float doppler_boost = clamp(look_doppler_boost, 0.0, 2.5);

    float gamma = 1.0/sqrt(max(1.0-dot(cam_vel,cam_vel), 0.0001));
    ray_doppler_factor = gamma*(1.0 + dot(ray,-cam_vel));
    {{#beaming}}
    {{#physical_beaming}}
    // Observer beaming is already accounted for via transfer_factor
    // in each emission source. The temperature shift of blackbody sources
    // captures the full Doppler effect. No additional ray_intensity scaling.
    {{/physical_beaming}}
    {{^physical_beaming}}
    float beaming_factor = clamp(ray_doppler_factor, 0.84, 1.16);
    ray_intensity /= pow(beaming_factor, 0.65 + 0.75*doppler_boost);
    {{/physical_beaming}}
    {{/beaming}}
    {{^doppler_shift}}
    ray_doppler_factor = 1.0;
    {{/doppler_shift}}

    float step = 0.01;
    vec4 color = vec4(0.0,0.0,0.0,1.0);
    float vol_transmittance = 1.0; // volumetric ray transmittance (Beer-Lambert)

    // initial conditions
    float u = 1.0 / length(pos), old_u;
    float u0 = u;
    float observer_speed = length(cam_vel);
    float observer_static_lapse = 1.0;
    if (interior_mode < 0.5 && observer_speed < 1e-4) {
        // Static camera rays live in the local orthonormal frame. In
        // Schwarzschild coordinates the radial basis is compressed by
        // sqrt(1 - r_s/r), so the radial component must carry that factor when
        // converted to the Binet initial condition. Without it the near-horizon
        // hover escape cone stays far too wide.
        observer_static_lapse = sqrt(max(1.0 - u, 1e-4));
    }

    vec3 normal_vec = normalize(pos);
    // Tangential component of ray (perpendicular to radial direction).
    // For near-radial rays the cross product is tiny → degenerate tangent_vec.
    // Guard with a fallback perpendicular direction.
    vec3 ray_perp = cross(cross(normal_vec, ray), normal_vec);
    float ray_perp_len = length(ray_perp);
    vec3 tangent_vec;
    if (ray_perp_len > 1e-6) {
        tangent_vec = ray_perp / ray_perp_len;
    } else {
        // Purely radial ray — pick an arbitrary perpendicular direction
        tangent_vec = abs(normal_vec.y) < 0.9
            ? normalize(cross(normal_vec, vec3(0.0, 1.0, 0.0)))
            : normalize(cross(normal_vec, vec3(1.0, 0.0, 0.0)));
        ray_perp_len = 1e-6;
    }
    vec3 spin_axis = vec3(0.0, 0.0, 1.0);
    float spin_alignment = clamp(dot(safe_normalize(cross(pos, ray)), spin_axis), -1.0, 1.0);
    float frame_drag_phase = 0.0;

    // du/dφ: rate of inverse-radius change per orbit angle.
    // For near-radial rays (small tangential component), du is large but finite.
    float radial_component = dot(ray, normal_vec) * observer_static_lapse;
    float dot_tang = dot(ray, tangent_vec);
    float du = (abs(dot_tang) > 1e-6)
        ? -radial_component / dot_tang * u
        : -sign(dot(ray, normal_vec)) * u * 200.0; // nearly radial — large du
    // Clamp: ±200 keeps near-radial rays integrable while covering all
    // physically relevant escape trajectories from inside the horizon.
    // (Critical du for escape ≈ sqrt(2|u²-u³|) which is ~124 at u=20.)
    du = clamp(du, -200.0, 200.0);
    float du0 = du;

    float phi = 0.0;

    float t = time;
    float turb_t = loopable_turbulence_time(time + turbulence_time_offset);
    float dt = 1.0;
    bool shadow_capture = false;

    {{^light_travel_time}}
    float planet_ang0 = t * PLANET_ORBITAL_ANG_VEL;
    vec3 planet_pos0 = vec3(cos(planet_ang0), sin(planet_ang0), 0)*PLANET_DISTANCE;
    {{/light_travel_time}}

    vec3 old_pos = pos;

    // ── Kerr state variables ──
    // Experimental true Kerr geodesics exist in the codebase but are not
    // exposed publicly. User-facing spin modes still trace photons with the
    // Schwarzschild Binet solver plus perturbative frame dragging.
    float kerr_a_phys = kerr_spin_a();
    {{#kerr_full_geodesic}}
    bool use_kerr = (abs(kerr_a_phys) > 0.001 && interior_mode < 0.5);
    {{/kerr_full_geodesic}}
    {{^kerr_full_geodesic}}
    bool use_kerr = false;
    {{/kerr_full_geodesic}}
    float kerr_r = 0.0, kerr_cth = 0.0, kerr_phi = 0.0;
    float kerr_pr = 0.0, kerr_pcth = 0.0;
    float kerr_xi = 0.0, kerr_eta = 0.0;
    float kerr_r_horizon = 1.0;

    if (use_kerr) {
        kerr_init(pos, ray, kerr_a_phys,
            kerr_xi, kerr_eta,
            kerr_r, kerr_cth, kerr_phi,
            kerr_pr, kerr_pcth);
        kerr_r_horizon = kerr_horizon_radius(kerr_a_phys);
    }

    // ── Interior mode: analytical escape classification ─────────────────
    // The Binet conserved energy E = (du/dφ)² + u² − u³ determines escape.
    // Potential barrier maximum at u = 2/3 (photon sphere): E_crit = 4/27.
    //   • du ≥ 0  (inward ray)           → always captured (singularity)
    //   • du < 0  and  E ≤ 4/27          → reflected by barrier → captured
    //   • du < 0  and  E > 4/27          → escapes through horizon
    // Pre-classifying captured rays avoids wasting integration steps on
    // rays that can never escape, completely eliminating the numerical
    // artefacts that produced the blue-blob rendering.
    if (interior_mode > 0.5) {
        if (du >= 0.0) {
            shadow_capture = true;  // inward → singularity
        } else {
            float E_binet = du*du + u*u*(1.0 - u);
            if (E_binet <= 4.0/27.0) {
                shadow_capture = true;  // not enough energy to escape
            }
        }
    } else if (observer_speed < 1e-4 && u > 2.0/3.0) {
        // Exact Schwarzschild escape test for a static observer inside the
        // photon sphere. This removes the numerical leakage that made the
        // visible sky stop shrinking during hover.
        if (du >= 0.0) {
            shadow_capture = true;
        } else {
            float E_binet = du*du + u*u*(1.0 - u);
            if (E_binet <= 4.0/27.0) {
                shadow_capture = true;
            }
        }
    }

    float E_binet0 = du0*du0 + u0*u0*(1.0 - u0);
    bool allow_extended_trace = !use_kerr &&
        interior_mode < 0.5 &&
        observer_speed < 1e-4 &&
        u0 > 2.0/3.0 &&
        du0 < 0.0 &&
        E_binet0 > 4.0/27.0;
    bool trace_finished = false;

    for (int pass = 0; pass < 3; pass++) {
        if (trace_finished) break;
        if (pass > 0 && !allow_extended_trace) break;

        for (int j=0; j < NSTEPS; j++) {
        if (shadow_capture) break;  // pre-classified capture → skip loop

        if (use_kerr) {
        // ── True Kerr geodesic integration ─────────────────────────────
        step = MAX_REVOLUTIONS * 2.0*M_PI / float(NSTEPS);

        // Convert angular budget to Mino-time step: dσ ≈ dφ / |dφ/dσ|.
        // At large r the azimuthal rate dφ/dσ scales with r, so without
        // this each Mino step would cover huge angles.  Near the horizon
        // dφ/dσ → ∞ (BL coordinate singularity) but the actual radial
        // traversal in Mino time is finite, so cap the scaling to prevent
        // the step from collapsing to zero before the photon can capture.
        float dphi_dsigma = abs(kerr_phi_dot(kerr_r, kerr_cth, kerr_a_phys, kerr_xi));
        float mino_scale = clamp(dphi_dsigma, 1.0, 40.0);
        step /= mino_scale;

        // Adaptive refinement near photon sphere for shadow accuracy
        float ps_dist = kerr_r - 1.5;
        step *= 1.0 - 0.5 * exp(-12.0*ps_dist*ps_dist);

        // Limit radial change per step to 30% of current r
        float dr_per_step = abs(kerr_pr) * step;
        if (dr_per_step > kerr_r * 0.3 && abs(kerr_pr) > 0.001) {
            step = kerr_r * 0.3 / abs(kerr_pr);
        }
        step = max(step, 1e-6);

        old_u = u;
        old_pos = pos;

        integrate_kerr_step(kerr_r, kerr_cth, kerr_phi,
            kerr_pr, kerr_pcth,
            step, kerr_a_phys, kerr_xi, kerr_eta);

        // Capture: photon reached the event horizon
        if (kerr_r <= kerr_r_horizon + 0.02 || kerr_r < 0.01) {
            shadow_capture = true;
            break;
        }

        // Reconstruct 3D Cartesian position from Boyer-Lindquist (r, θ, φ)
        // MUST happen before escape check so `pos` is valid for background ray
        float sth_k = sqrt(max(1.0 - kerr_cth*kerr_cth, 0.0));
        float rho_perp = sqrt(kerr_r*kerr_r + kerr_a_phys*kerr_a_phys) * sth_k;
        pos = vec3(rho_perp * cos(kerr_phi),
                   rho_perp * sin(kerr_phi),
                   kerr_r * kerr_cth);
        u = 1.0 / max(kerr_r, 0.001);

        // Escape: photon past observer and heading outward
        if (kerr_r > length(cam_pos) * 1.1 && kerr_pr > 0.0) {
            break;
        }

        } else {
        // ── Public Schwarzschild Binet integration (all exposed modes;
        //    also reused inside the horizon) ────────────────────────────
        step = MAX_REVOLUTIONS * 2.0*M_PI / float(NSTEPS);

        // adaptive step size based on rate of change
        float max_rel_u_change = max(1.0-log(max(u, 0.0001)), 0.05)*10.0 / float(NSTEPS);
        if (interior_mode > 0.5) {
            max_rel_u_change *= 3.0 + 2.0 * max(u - 1.0, 0.0);
            if (abs(du) > abs(max_rel_u_change*u) / step) {
                step = max_rel_u_change*u/abs(du);
            }
        } else {
            if ((du > 0.0 || (du0 < 0.0 && u0/u < 5.0)) && abs(du) > abs(max_rel_u_change*u) / step) {
                step = max_rel_u_change*u/abs(du);
            }
        }

        float u_photon_sphere = 0.667;
        float photon_sphere_proximity = exp(-12.0 * (u - u_photon_sphere) * (u - u_photon_sphere));
        step *= 1.0 - 0.7 * photon_sphere_proximity;

        if (interior_mode > 0.5) {
            float horizon_proximity = exp(-8.0 * (u - 1.0) * (u - 1.0));
            step *= 1.0 - 0.6 * horizon_proximity;
            step = max(step, 0.5 * M_PI / float(NSTEPS));
        }

        old_u = u;

        {{#light_travel_time}}
        {{#gravitational_time_dilation}}
        dt = sqrt(max(du*du + u*u*(1.0-u), 0.0001))/max(abs(u*u*(1.0-u)), 0.0001)*step;
        {{/gravitational_time_dilation}}
        {{/light_travel_time}}

        integrate_geodesic_step(u, du, step, spin_alignment);

        if (u < 0.0) {
            if (interior_mode > 0.5) shadow_capture = true;
            break;
        }
        if (interior_mode < 0.5) {
            if (u >= 1.0) {
                shadow_capture = true;
                break;
            }
        } else {
            if (u >= 20.0) {
                shadow_capture = true;
                break;
            }
        }

        phi += step;

        old_pos = pos;
        vec3 planar_pos = (cos(phi)*normal_vec + sin(phi)*tangent_vec)/u;

        float drag_u = clamp(u, 0.0, 1.15);
        float frame_drag_step = photon_spin_lensing_scale *
            bh_rotation_enabled * bh_spin * bh_spin_strength *
            spin_alignment * step * 0.85 * drag_u*drag_u*drag_u;
        frame_drag_phase += frame_drag_step;

        pos = rotate_about_z(planar_pos, frame_drag_phase);
        } // end if (use_kerr) / else

        ray = pos-old_pos;
        float solid_isec_t = 2.0;
        float ray_l = length(ray);

        {{#light_travel_time}}
        {{#gravitational_time_dilation}}
        float mix = smooth_step(1.0/u, 8.0);
        dt = mix*ray_l + (1.0-mix)*dt;
        {{/gravitational_time_dilation}}
        {{^gravitational_time_dilation}}
        dt = ray_l;
        {{/gravitational_time_dilation}}
        {{/light_travel_time}}

        {{#planetEnabled}}
        if (
            (
                old_pos.z * pos.z < 0.0 ||
                min(abs(old_pos.z), abs(pos.z)) < PLANET_RADIUS
            ) &&
            max(u, old_u) > 1.0/(PLANET_DISTANCE+PLANET_RADIUS) &&
            min(u, old_u) < 1.0/(PLANET_DISTANCE-PLANET_RADIUS)
        ) {

            {{#light_travel_time}}
            float planet_ang0 = t * PLANET_ORBITAL_ANG_VEL;
            vec3 planet_pos0 = vec3(cos(planet_ang0), sin(planet_ang0), 0)*PLANET_DISTANCE;
            {{/light_travel_time}}

            vec4 planet_isec = planet_intersection(old_pos, ray, t, dt,
                    planet_pos0, ray_doppler_factor);
            if (planet_isec.w > 0.0) {
                solid_isec_t = planet_isec.w;
                planet_isec.w = 1.0;
                color += planet_isec;
            }
        }
        {{/planetEnabled}}
