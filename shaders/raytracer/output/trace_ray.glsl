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
    float dot_tang = dot(ray, tangent_vec);
    float du = (abs(dot_tang) > 1e-6)
        ? -dot(ray, normal_vec) / dot_tang * u
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
    // When kerr_full_geodesic is active and spin > 0. true Carter (1968)
    // Mino-time integration is used, producing the correct D-shaped
    // shadow for spinning black holes.  Otherwise the Binet approximation
    // with frame-drag term is used — exact for Schwarzschild.
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
    }

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
        // ── Schwarzschild Binet integration (a ≈ 0 or interior) ────────
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
        float frame_drag_step = bh_rotation_enabled * bh_spin * bh_spin_strength *
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

        {{#accretion_thin_disk}}
        // Disk crossing detection with sub-step refinement.
        // Near the photon sphere (r ≈ 1.5 r_s) photon paths wind tightly
        // and a single step can skip over the z=0 plane.  When the step
        // subtends a large angle and we are close to the photon sphere,
        // subdivide the segment and test each sub-segment for a crossing
        // to reliably capture secondary and tertiary Einstein rings.
        int disk_sub_steps = 1;
        if (abs(u - 0.667) < 0.15 && step > 0.12) {
            disk_sub_steps = 4;  // 4 sub-checks near photon sphere
        }
        {
            vec3 sub_old = old_pos;
            vec3 sub_step_vec = (pos - old_pos) / float(disk_sub_steps);
            for (int ds = 0; ds < 4; ds++) {
                if (ds >= disk_sub_steps) break;
                vec3 sub_new = old_pos + sub_step_vec * float(ds + 1);
                if (sub_old.z * sub_new.z < 0.0) {
            vec3 sub_ray = sub_new - sub_old;
            float sub_ray_z = sub_ray.z;
            if (abs(sub_ray_z) < 1e-6) {
                sub_old = sub_new;
                continue;
            }
            float acc_isec_t = -sub_old.z / sub_ray_z;
            if (acc_isec_t < solid_isec_t) {
                vec3 isec = sub_old + sub_ray*acc_isec_t;

                float r = length(isec);
                {{#grmhd_enabled}}
                // GRMHD: magnetic stress allows emission inside ISCO
                // (Noble+ 2010, Penna+ 2010: plunging region contributes ~10-30%)
                float r_inner_grmhd = max(ACCRETION_MIN_R * 0.7, 1.05);
                if (r > r_inner_grmhd && r < ACCRETION_MIN_R + ACCRETION_WIDTH) {
                {{/grmhd_enabled}}
                {{^grmhd_enabled}}
                if (r > ACCRETION_MIN_R && r < ACCRETION_MIN_R + ACCRETION_WIDTH) {
                {{/grmhd_enabled}}
                    float angle = atan(isec.x, isec.y);

                    {{#grmhd_enabled}}
                    // GRMHD thin disk: Shakura-Sunyaev base + subtle GRMHD corrections.
                    // Thin disks are radiatively efficient → very mild R_high temp correction
                    // (3%) to preserve the well-established Shakura-Sunyaev color profile.
                    // Synchrotron/non-thermal add subtle modulation, not brightness boosts.
                    float gas_temp_g = accretion_temperature(r);
                    float disk_h_g = r * 0.05;
                    float beta_g = grmhd_plasma_beta(r, 0.0, disk_h_g);
                    // Add t and angle parameters needed by the updated function signature
                    float Te_corr = 0.97 + grmhd_electron_temp_ratio(beta_g) * 0.03;  // 3% R_high correction (thin disk is efficient)
                    float temperature = grmhd_electron_temperature(gas_temp_g, r, 0.0, disk_h_g, turb_t, angle) * Te_corr * gravitational_shift(r);
                    
                    // GRMHD turbulence: MRI density fluctuations with log-normal
                    // PDF + spiral arms (Sorathia+ 2012, Hawley+ 2013)
                    float r_norm_g = (r - ACCRETION_MIN_R) / ACCRETION_WIDTH;
                    float edge_fade_g = smoothstep(0.02, 0.18, r_norm_g) *
                        (1.0 - smoothstep(0.78, 1.0, r_norm_g));
                    float turbulence = grmhd_mri_turbulence(r, angle, turb_t) * edge_fade_g;

                    // Synchrotron and non-thermal corrections
                    float rho_g = grmhd_density(r, 1.0);
                    float B_g = grmhd_B_field(rho_g, gas_temp_g, beta_g);
                    float nonthermal_g = grmhd_nonthermal_boost(grmhd_electron_kappa);
                    float sync_corr_g = grmhd_synchrotron_correction(B_g, rho_g);
                    // ISCO magnetic stress: emission extends inside ISCO
                    float isco_stress_g = grmhd_isco_stress_factor(r);

                    // Attenuated GRMHD corrections: parameters are responsive but
                    // preserve the base Shakura-Sunyaev brightness profile.
                    // At defaults (κ=5, B=1): ~8% deviation from non-GRMHD.
                    float nt_mod_g = 1.0 + 0.1 * (nonthermal_g - 1.0);
                    float sync_mod_g = 1.0 + 0.1 * (sync_corr_g - 1.0);

                    float inner_glow = exp(-8.0 * (r - ACCRETION_MIN_R));
                    float accretion_intensity = ACCRETION_BRIGHTNESS * accretion_flux_profile(r) *
                        turbulence * nt_mod_g * sync_mod_g * isco_stress_g
                        * (1.0 + 0.7*inner_glow);
                    {{/grmhd_enabled}}
                    {{^grmhd_enabled}}
                    float temperature = accretion_temperature(r) * gravitational_shift(r);
                    float turbulence = accretion_emissivity(r, angle, turb_t);
                    float inner_glow = exp(-8.0 * (r - ACCRETION_MIN_R));
                    float accretion_intensity = ACCRETION_BRIGHTNESS * accretion_flux_profile(r) *
                        turbulence * (1.0 + 0.7*inner_glow);
                    {{/grmhd_enabled}}

                    // Limb darkening: Eddington approximation for an optically
                    // thick disk.  I(μ) = I(1)·(2/5 + 3/5·μ) where μ = cos(θ)
                    // is the angle between the photon and the disk normal (z-axis).
                    float cos_emission_angle = abs(ray.z) / max(ray_l, 1e-6);
                    accretion_intensity *= 0.4 + 0.6 * cos_emission_angle;

                    vec3 accretion_v;
                    {{#kerr_fast_mode}}
                    accretion_v = vec3(-isec.y, isec.x, 0.0) / (r * sqrt(2.0*(r-1.0)));
                    {{/kerr_fast_mode}}
                    {{#kerr_full_velocity}}
                    float rg_r = max(2.0*r, 1.0002); // convert r_s units to M units
                    float a_M = bh_rotation_enabled * bh_spin;
                    float omega_M = 1.0 / (pow(rg_r, 1.5) + a_M);
                    float v_phi = clamp(rg_r * omega_M, -0.995, 0.995);
                    vec2 xy = vec2(isec.x, isec.y);
                    float cyl_r = max(length(xy), 1e-4);
                    vec3 e_phi = vec3(-xy.y/cyl_r, xy.x/cyl_r, 0.0);
                    accretion_v = e_phi * v_phi;
                    {{/kerr_full_velocity}}

                    gamma = 1.0/sqrt(max(1.0-dot(accretion_v,accretion_v), 0.0001));
                    float doppler_factor = gamma*(1.0+dot(ray/ray_l,accretion_v));
                    float transfer_factor = max(ray_doppler_factor*doppler_factor, 0.05);
                    {{#beaming}}
                    {{#physical_beaming}}
                    // Liouville invariant: I_obs = D³ · I_emit.
                    // The spectrum texture stores chromaticity (normalized color),
                    // so the temperature shift only handles hue — D³ intensity
                    // boost must always be applied explicitly.
                    accretion_intensity /= pow(clamp(transfer_factor, 0.05, 20.0), 3.0);
                    {{/physical_beaming}}
                    {{^physical_beaming}}
                    // Cinematic beaming: softened for artistic rendering
                    float clamped_doppler = clamp(transfer_factor, 0.62, 1.48);
                    accretion_intensity /= pow(clamped_doppler, 1.05 + 1.10*doppler_boost);
                    {{/physical_beaming}}
                    {{/beaming}}
                    {{#doppler_shift}}
                    temperature /= transfer_factor;
                    {{/doppler_shift}}

                    vec4 thermal_color = BLACK_BODY_COLOR(temperature);
                    accretion_intensity *= look_disk_gain;
                    accretion_intensity *= 1.0 + look_glow * (0.22 + 0.78*inner_glow);
                    color += vec4(thermal_color.rgb * accretion_intensity, 1.0);
                }
            }
                }
                sub_old = sub_new;
            }
        }
        {{/accretion_thin_disk}}

        {{#accretion_thick_torus}}
        {
            // ADAF/RIAF volumetric emission: optically thin radiative transfer
            // dI/ds = j - alpha*I  (emission minus absorption)
            {{#grmhd_enabled}}
            // GRMHD: evaluate torus with stronger height modulation (2×) for
            // visible azimuthal warps and buoyant features from MRI + Parker
            // instability.  The torus is geometrically thick, so larger warps
            // are physically expected (Liska+ 2022: H/R variation ~15-30%).
            float _cr_t = max(length(pos.xy), 1e-4);
            float _a_t = atan(pos.x, pos.y);
            float _hmod_raw_t = grmhd_height_modulation(_cr_t, _a_t, turb_t);
            // Amplify: base function gives ±15%, double it for ±30% warps
            float _hmod_t = 1.0 + 2.0 * (_hmod_raw_t - 1.0);
            vec3 _warp_t = vec3(pos.xy, pos.z / max(_hmod_t, 0.2));
            float torus_j = torus_local_emissivity(_warp_t);
            {{/grmhd_enabled}}
            {{^grmhd_enabled}}
            float torus_j = torus_local_emissivity(pos);
            {{/grmhd_enabled}}
            if (torus_j > 0.001 && vol_transmittance > 0.005) {
                float path_len = length(pos - old_pos);
                float r3d = max(length(pos), 1.001);
                float cyl_r_t = max(length(pos.xy), 1e-4);
                float angle_t = atan(pos.x, pos.y);

                {{#grmhd_enabled}}
                // GRMHD two-temperature torus model (EHT-calibrated)
                // Temperature: β-dependent variation for color shifts.
                float h_torus_g = max(cyl_r_t * torus_h_ratio, 0.01);
                float gas_temp_torus = torus_temperature(r3d);
                float temperature_t = grmhd_electron_temperature(gas_temp_torus, cyl_r_t, pos.z, h_torus_g, turb_t, angle_t)
                                    * gravitational_shift_static(r3d);

                // 3D volumetric turbulence with sqrt() to moderate extremes.
                // The raw 3D FBM has high dynamic range (~100×) which produces
                // dramatic filaments in the slim disk (where Beer-Lambert
                // absorption renders surface features).  For the optically thin
                // ADAF torus, sqrt compresses this to ~10×, keeping visible 3D
                // structure without the tornado/nebula artefact.
                float raw_turb_t = grmhd_3d_density_turbulence(pos, turb_t);
                float turbulence_t = sqrt(max(raw_turb_t, 0.01));

                // Magnetic field, synchrotron, and non-thermal corrections
                float beta_torus = grmhd_plasma_beta(cyl_r_t, pos.z, h_torus_g);
                float rho_torus = grmhd_density(cyl_r_t, 1.0);
                float B_torus = grmhd_B_field(rho_torus, gas_temp_torus, beta_torus);
                float nonthermal_t = grmhd_nonthermal_boost(grmhd_electron_kappa);
                float sync_corr_t = grmhd_synchrotron_correction(B_torus, rho_torus);

                // Attenuated GRMHD corrections: parameters responsive but
                // torus brightness stays close to the non-GRMHD appearance.
                float nt_mod_t = 1.0 + 0.15 * (nonthermal_t - 1.0);
                float sync_mod_t = 1.0 + 0.15 * (sync_corr_t - 1.0);

                // Emissivity: base torus × 3D turbulence × subtle corrections
                float j_eff = ACCRETION_BRIGHTNESS * 0.05 * torus_j * turbulence_t
                            * nt_mod_t * sync_mod_t;

                // Absorption: turbulence-modulated (denser clumps are more opaque)
                float alpha_abs = torus_opacity * torus_j * pow(turbulence_t, 0.5)
                                * (0.5 + 0.5 * grmhd_density_scale);
                {{/grmhd_enabled}}
                {{^grmhd_enabled}}
                float temperature_t = torus_temperature(r3d) * gravitational_shift_static(r3d);
                float turbulence_t = accretion_turbulence(cyl_r_t, angle_t, turb_t);

                // Emission coefficient: optically thin ADAF has low emissivity
                float j_eff = ACCRETION_BRIGHTNESS * 0.05 * torus_j * turbulence_t;
                // Absorption: user-configurable opacity (default low for optically thin)
                float alpha_abs = torus_opacity * torus_j;
                {{/grmhd_enabled}}
                float tau_step = alpha_abs * path_len;
                float step_T = exp(-tau_step);

                // Sub-Keplerian gas velocity (ADAF: ~50% of Keplerian)
                vec3 accretion_v_t;
                {{#kerr_fast_mode}}
                float v_kep_t = 1.0 / sqrt(2.0 * max(cyl_r_t - 1.0, 0.01));
                float v_sub_t = clamp(0.5 * v_kep_t, 0.0, 0.95);
                accretion_v_t = vec3(-pos.y, pos.x, 0.0) / cyl_r_t * v_sub_t;
                {{/kerr_fast_mode}}
                {{#kerr_full_velocity}}
                float rg_cyl_t = max(2.0 * cyl_r_t, 1.0002);
                float a_M_t = bh_rotation_enabled * bh_spin;
                float omega_sub_t = 0.5 / (pow(rg_cyl_t, 1.5) + a_M_t);
                float v_phi_t = clamp(rg_cyl_t * omega_sub_t, -0.95, 0.95);
                vec3 e_phi_t = vec3(-pos.y, pos.x, 0.0) / cyl_r_t;
                accretion_v_t = e_phi_t * v_phi_t;
                {{/kerr_full_velocity}}

                gamma = 1.0/sqrt(max(1.0-dot(accretion_v_t,accretion_v_t), 0.0001));
                float doppler_factor_t = gamma*(1.0+dot(ray/ray_l,accretion_v_t));
                float transfer_factor_t = max(ray_doppler_factor*doppler_factor_t, 0.05);

                float torus_intensity = j_eff * path_len;
                {{#beaming}}
                {{#physical_beaming}}
                // Liouville D³: chromaticity texture doesn't encode brightness,
                // so intensity boost must always be applied explicitly.
                torus_intensity /= pow(clamp(transfer_factor_t, 0.05, 20.0), 3.0);
                {{/physical_beaming}}
                {{^physical_beaming}}
                float cd_t = clamp(transfer_factor_t, 0.62, 1.48);
                torus_intensity /= pow(cd_t, 1.05 + 1.10*doppler_boost);
                {{/physical_beaming}}
                {{/beaming}}
                {{#doppler_shift}}
                temperature_t /= transfer_factor_t;
                {{/doppler_shift}}

                vec4 thermal_color_t = BLACK_BODY_COLOR(temperature_t);
                torus_intensity *= look_disk_gain;

                // Front-to-back compositing with Beer-Lambert absorption
                color.rgb += vol_transmittance * thermal_color_t.rgb * torus_intensity;
                vol_transmittance *= step_T;
            }
        }
        {{/accretion_thick_torus}}

        {{#accretion_slim_disk}}
        {
            // Slim disk volumetric emission: optically thick radiative transfer
            // Super-Eddington flow: high absorption → surface-like rendering
            {{#grmhd_enabled}}
            // GRMHD: height modulation from MRI + radiation-driven warps
            float _cr_s = max(length(pos.xy), 1e-4);
            float _a_s = atan(pos.x, pos.y);
            float _hmod_s = grmhd_height_modulation(_cr_s, _a_s, turb_t);
            vec3 _warp_s = vec3(pos.xy, pos.z / max(_hmod_s, 0.2));
            float slim_j = slim_disk_local_emissivity(_warp_s);
            {{/grmhd_enabled}}
            {{^grmhd_enabled}}
            float slim_j = slim_disk_local_emissivity(pos);
            {{/grmhd_enabled}}
            if (slim_j > 0.001 && vol_transmittance > 0.005) {
                float path_len_s = length(pos - old_pos);
                float cyl_r_s = max(length(pos.xy), 1e-4);
                float r3d_s = max(length(pos), 1.001);
                float angle_s = atan(pos.x, pos.y);

                {{#grmhd_enabled}}
                // GRMHD slim disk: super-Eddington base model remains physically accurate.
                // GRMHD adds: mild two-temperature correction (30% R_high) + ISCO magnetic
                // stress + synchrotron/non-thermal corrections responsive to all parameters.
                float gas_temp_slim = slim_disk_temperature(cyl_r_s);
                float h_slim_g = slim_disk_height(cyl_r_s);
                float beta_slim = grmhd_plasma_beta(cyl_r_s, pos.z, h_slim_g);
                float Te_ratio_slim = grmhd_electron_temp_ratio(beta_slim);
                // Manual lerp: avoid 'mix' variable shadow from line 264
                float Te_corr_slim = 1.0 * 0.7 + Te_ratio_slim * 0.3;  // 30% correction (partially thermalized)
                float temperature_s = grmhd_electron_temperature(gas_temp_slim, cyl_r_s, pos.z, h_slim_g, turb_t, angle_s) * Te_corr_slim * gravitational_shift_static(r3d_s);

                // 3D volumetric FBM turbulence: filamentary density structure
                float turbulence_s = grmhd_3d_density_turbulence(pos, turb_t);
                
                // ISCO stress: adds ~10-30% extra luminosity from plunging region
                float isco_stress_s = grmhd_isco_stress_factor(cyl_r_s);

                // Synchrotron and non-thermal corrections (radiation-dominated
                // → gentler response than ADAF, normalized so defaults ≈ 1.0)
                float rho_slim_g = grmhd_density(cyl_r_s, 1.0);
                float B_slim_g = grmhd_B_field(rho_slim_g, gas_temp_slim, beta_slim);
                float nonthermal_slm = grmhd_nonthermal_boost(grmhd_electron_kappa);
                float sync_slm = grmhd_synchrotron_correction(B_slim_g, rho_slim_g);
                // sqrt + /2.5 normalization keeps defaults ~1.0, avoids breaking
                // the existing slim disk appearance while making all params responsive
                float grmhd_param_slim = sqrt(nonthermal_slm * sync_slm / 2.5);

                float j_eff_s = ACCRETION_BRIGHTNESS * 0.9 * slim_j * turbulence_s * isco_stress_s * grmhd_param_slim;
                // Absorption: density-modulated (clumps are more opaque)
                float alpha_abs_s = slim_opacity * slim_j * pow(turbulence_s, 0.5);
                {{/grmhd_enabled}}
                {{^grmhd_enabled}}
                float temperature_s = slim_disk_temperature(cyl_r_s) * gravitational_shift_static(r3d_s);
                float turbulence_s = accretion_turbulence(cyl_r_s, angle_s, turb_t);

                // Emission coefficient: slim disk is bright (super-Eddington luminosity)
                float j_eff_s = ACCRETION_BRIGHTNESS * 0.9 * slim_j * turbulence_s;
                // User-configurable absorption: higher = more opaque surface-like
                float alpha_abs_s = slim_opacity * slim_j;
                {{/grmhd_enabled}}
                float tau_step_s = alpha_abs_s * path_len_s;
                float step_T_s = exp(-tau_step_s);
                // Source function: emission that survives self-absorption in this step
                // S = j/alpha; contribution = (1 - exp(-tau)) * S = (1-step_T) * j/alpha
                float source_contrib = (tau_step_s > 0.001)
                    ? (1.0 - step_T_s) * (j_eff_s / alpha_abs_s)
                    : j_eff_s * path_len_s;

                // Keplerian velocity outside ISCO, capped inside
                vec3 accretion_v_s;
                {{#kerr_fast_mode}}
                float v_kep_s = 1.0 / sqrt(2.0 * max(cyl_r_s - 1.0, 0.01));
                accretion_v_s = vec3(-pos.y, pos.x, 0.0) / cyl_r_s * clamp(v_kep_s, 0.0, 0.95);
                {{/kerr_fast_mode}}
                {{#kerr_full_velocity}}
                float rg_cyl_s = max(2.0 * cyl_r_s, 1.0002);
                float a_M_s = bh_rotation_enabled * bh_spin;
                float omega_s = 1.0 / (pow(rg_cyl_s, 1.5) + a_M_s);
                float v_phi_s = clamp(rg_cyl_s * omega_s, -0.95, 0.95);
                vec3 e_phi_s = vec3(-pos.y, pos.x, 0.0) / cyl_r_s;
                accretion_v_s = e_phi_s * v_phi_s;
                {{/kerr_full_velocity}}

                gamma = 1.0/sqrt(max(1.0-dot(accretion_v_s,accretion_v_s), 0.0001));
                float doppler_factor_s = gamma*(1.0+dot(ray/ray_l,accretion_v_s));
                float transfer_factor_s = max(ray_doppler_factor*doppler_factor_s, 0.05);

                float slim_intensity = source_contrib;
                {{#beaming}}
                {{#physical_beaming}}
                // Liouville D³: chromaticity texture doesn't encode brightness,
                // so intensity boost must always be applied explicitly.
                slim_intensity /= pow(clamp(transfer_factor_s, 0.05, 20.0), 3.0);
                {{/physical_beaming}}
                {{^physical_beaming}}
                float cd_s = clamp(transfer_factor_s, 0.62, 1.48);
                slim_intensity /= pow(cd_s, 1.05 + 1.10*doppler_boost);
                {{/physical_beaming}}
                {{/beaming}}
                {{#doppler_shift}}
                temperature_s /= transfer_factor_s;
                {{/doppler_shift}}

                vec4 thermal_color_s = BLACK_BODY_COLOR(temperature_s);
                slim_intensity *= look_disk_gain;

                // Front-to-back compositing with Beer-Lambert absorption
                color.rgb += vol_transmittance * thermal_color_s.rgb * slim_intensity;
                vol_transmittance *= step_T_s;
            }
        }
        {{/accretion_slim_disk}}

        {{#jet_enabled}}
        {
            float path_len_j = length(pos - old_pos);
            for (int jet_side = 0; jet_side < 2; jet_side++) {
                float sign_z = (jet_side == 0) ? 1.0 : -1.0;

                {{#jet_physical}}
                // --- Counter-jet disk occultation ---
                // The accretion disk (optically thick near equator) blocks
                // the receding counter-jet when viewed near edge-on.
                // Model as extra absorption for the counter-jet passing
                // through the disk plane.
                float occultation = 1.0;
                if (sign_z < 0.0) {
                    // Counter-jet: check if line of sight passes through
                    // a dense region near the equatorial plane
                    float abs_z_pos = abs(pos.z);
                    float cyl_r_occ = length(pos.xy);
                    // Disk is optically thick for r > ISCO, |z| < H(r)
                    // H/R ~ 0.1-0.3 for thin/slim disks
                    float disk_H = 0.2 * cyl_r_occ;
                    if (abs_z_pos < disk_H && cyl_r_occ > 1.5 && cyl_r_occ < 15.0) {
                        // Strong absorption in disk midplane
                        float tau_disk = 8.0 * (1.0 - abs_z_pos / disk_H);
                        occultation = exp(-tau_disk);
                    }
                }
                {{/jet_physical}}

                {{#grmhd_enabled}}
                {{#jet_physical}}
                float j_jet = grmhd_jet_emissivity(pos, sign_z);
                {{/jet_physical}}
                {{#jet_simple}}
                float j_jet = jet_emissivity(pos, sign_z);
                {{/jet_simple}}
                {{/grmhd_enabled}}
                {{^grmhd_enabled}}
                float j_jet = jet_emissivity(pos, sign_z);
                {{/grmhd_enabled}}
                if (j_jet > 0.001 && vol_transmittance > 0.005) {
                    // Jet bulk velocity
                    vec3 v_jet = jet_velocity(pos, sign_z);
                    float v2_jet = dot(v_jet, v_jet);
                    float gamma_jet = 1.0 / sqrt(max(1.0 - v2_jet, 0.0001));

                    // Doppler factor: δ = 1 / (Γ(1 - β·n̂))
                    vec3 ray_dir = ray / max(ray_l, 1e-6);
                    float doppler_jet = gamma_jet * (1.0 + dot(ray_dir, v_jet));

                    // Synchrotron beaming: I_obs = δ^(2+α) * I_em
                    // α ≈ 0.7 for typical synchrotron spectral index
                    float beam_factor = 1.0 / pow(max(doppler_jet, 0.02), 2.7);

                    {{#jet_simple}}
                    // Simple mode: single-temperature blackbody approximation
                    float r_jet = max(length(pos), 1.001);
                    float g_shift_j = gravitational_shift_static(r_jet);
                    float jet_T = 18000.0 * g_shift_j;
                    {{#doppler_shift}}
                    jet_T /= max(doppler_jet, 0.05);
                    {{/doppler_shift}}
                    vec4 jet_color = BLACK_BODY_COLOR(jet_T);

                    float alpha_jet = 0.005 * j_jet;
                    float tau_jet = alpha_jet * path_len_j;
                    float step_T_jet = exp(-tau_jet);

                    float jet_emit = jet_brightness * 0.35 * j_jet * beam_factor * path_len_j;
                    jet_emit *= look_disk_gain;

                    color.rgb += vol_transmittance * jet_color.rgb * jet_emit;
                    vol_transmittance *= step_T_jet;
                    {{/jet_simple}}

                    {{#jet_physical}}
                    // Physical mode: GRMHD-calibrated emission
                    float r_jet_p = max(length(pos), 1.001);
                    float z_abs = abs(pos.z * sign_z);
                    float g_shift_j = gravitational_shift_static(r_jet_p);
                    float sigma = magnetization(z_abs);

                    {{#grmhd_enabled}}
                    // Full GRMHD jet: BZ power + kink instability + cooling break
                    float jet_T = grmhd_jet_temperature(z_abs, r_jet_p, sigma) * g_shift_j;
                    {{/grmhd_enabled}}
                    {{^grmhd_enabled}}
                    // Effective temperature with synchrotron aging
                    float jet_T = jet_effective_temperature(z_abs, r_jet_p, sigma) * g_shift_j;
                    {{/grmhd_enabled}}

                    // Doppler temperature shift for approaching/receding jet
                    jet_T /= max(doppler_jet, 0.05);

                    vec4 jet_color = BLACK_BODY_COLOR(jet_T);

                    // Absorption: higher in sheath (denser), lower in spine
                    float sigma_abs_factor = 1.0 / (1.0 + 0.5 * sigma);
                    float alpha_jet = 0.012 * j_jet * sigma_abs_factor;
                    float tau_jet = alpha_jet * path_len_j;
                    float step_T_jet = exp(-tau_jet);

                    // Emissivity: GRMHD calibration with σ-dependent efficiency
                    float sigma_emit = 1.0 / (1.0 + 0.3 * sigma);
                    float jet_emit = jet_brightness * 0.5 * j_jet * beam_factor
                                   * path_len_j * sigma_emit;
                    jet_emit *= look_disk_gain;

                    // Apply counter-jet disk occultation
                    jet_emit *= occultation;

                    color.rgb += vol_transmittance * jet_color.rgb * jet_emit;
                    vol_transmittance *= step_T_jet;
                    {{/jet_physical}}
                }
            }
        }
        {{/jet_enabled}}

        {{#light_travel_time}}
        t -= dt;
        {{/light_travel_time}}

        if (solid_isec_t <= 1.0) u = 100.0; // break (exceeds both exterior and interior thresholds)

        if (!use_kerr) {
        {{#kerr_fast_mode}}
        if (interior_mode < 0.5) {
            float capture_u = 1.0 + bh_rotation_enabled * bh_spin * bh_spin_strength *
                spin_alignment * 0.12;
            capture_u = clamp(capture_u, 0.82, 1.18);
            if (u > capture_u) {
                shadow_capture = true;
                break;
            }
        } else if (u >= 20.0) {
            shadow_capture = true;
            break;
        }
        {{/kerr_fast_mode}}
        {{#kerr_full_core}}
        // Realtime Kerr: capture at horizon or singularity
        if (interior_mode > 0.5 && u >= 20.0) {
            shadow_capture = true;
            break;
        }
        {{/kerr_full_core}}
        {{#kerr_full_geodesic}}
        // Full Kerr geodesic fallback (Binet): interior capture
        if (interior_mode > 0.5 && u >= 20.0) {
            shadow_capture = true;
            break;
        }
        {{/kerr_full_geodesic}}
        } // end if (!use_kerr)
    }

    // Background sky: show for rays that escaped to far field.
    // Interior escaped rays use analytical exit direction (Binet tangent)
    // which avoids the numerical drift of pos − old_pos over many steps.
    if (!shadow_capture && u < 1.0) {

        if (interior_mode > 0.5) {
            // Analytical exit direction from the Binet parametrisation:
            //   d(pos)/dφ  =  φ̂/u  −  r̂·(du/u²)
            // where r̂ and φ̂ are the radial and tangential unit vectors
            // rotated to the current orbit angle φ.
            vec3 r_hat  = cos(phi)*normal_vec + sin(phi)*tangent_vec;
            vec3 phi_hat = -sin(phi)*normal_vec + cos(phi)*tangent_vec;
            vec3 exit_dir = phi_hat / max(u, 0.001)
                          - r_hat  * du / max(u*u, 0.001);
            exit_dir = rotate_about_z(exit_dir, frame_drag_phase);
            ray = normalize(exit_dir);
        } else {
            ray = normalize(pos - old_pos);
        }

        vec2 tex_coord = sphere_map(ray * BG_COORDS);
        float t_coord;

        // Interior blueshift: escaped photons gain energy climbing out
        // of the gravitational potential well.  Boost ~ 1/√|1 − 1/r|.
        float interior_boost = 1.0;
        if (interior_mode > 0.5) {
            float r_obs = max(length(cam_pos), 0.08);
            interior_boost = 1.0 + 0.5 * sqrt(abs(1.0/r_obs - 1.0));
            interior_boost = min(interior_boost, 4.0);
        }

        // Gravitational blueshift for background sky: light from infinity
        // gains energy falling into the gravitational potential well.
        // For a static (hovering) observer at r:
        //   f_obs/f_emit = 1/sqrt(1 - r_s/r) = 1/grav_blueshift_factor
        // Combined with kinematic Doppler from observer orbital/dive motion.
        float bg_doppler = ray_doppler_factor * grav_blueshift_factor;

        // Liouville invariant (I_nu / nu^3 = const along a ray):
        // D = f_obs/f_emit = 1/bg_doppler is the TOTAL frequency ratio
        // from infinity to the observer, combining gravitational blueshift
        // with kinematic Doppler.  Intensity scales as D^3.
        // The spectrum texture stores chromaticity (normalized colour), not
        // absolute Planck radiance, so the temperature shift alone only
        // changes the hue.  The D^3 brightness factor must be applied
        // explicitly.  Clamping prevents blow-ups for extreme D values
        // (near-horizon freefall where T shifts far into UV anyway).
        float safe_bg_doppler = max(abs(bg_doppler), 0.01);
        float bg_D = 1.0 / safe_bg_doppler;
        float bg_boost = min(bg_D * bg_D * bg_D, 10000.0);

        vec4 star_color = texture2D(star_texture, tex_coord);
        if (star_color.r > 0.0) {
            t_coord = (STAR_MIN_TEMPERATURE +
                (STAR_MAX_TEMPERATURE-STAR_MIN_TEMPERATURE) * star_color.g)
                 / safe_bg_doppler;

            color += BLACK_BODY_COLOR(t_coord) * star_color.r * STAR_BRIGHTNESS * look_star_gain * vol_transmittance * interior_boost * bg_boost;
        }

        color += galaxy_color(tex_coord, bg_doppler) * GALAXY_BRIGHTNESS * look_galaxy_gain * vol_transmittance * interior_boost * bg_boost;
    }

    return color*ray_intensity;
}
