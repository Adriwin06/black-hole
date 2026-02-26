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
    float dt = 1.0;
    bool shadow_capture = false;

    {{^light_travel_time}}
    float planet_ang0 = t * PLANET_ORBITAL_ANG_VEL;
    vec3 planet_pos0 = vec3(cos(planet_ang0), sin(planet_ang0), 0)*PLANET_DISTANCE;
    {{/light_travel_time}}

    vec3 old_pos = pos;

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
        {{#kerr_fast_mode}}
        step = MAX_REVOLUTIONS * 2.0*M_PI / float(NSTEPS);

        // adaptive step size based on rate of change
        float max_rel_u_change = max(1.0-log(max(u, 0.0001)), 0.05)*10.0 / float(NSTEPS);
        // Interior mode: always apply adaptive stepping (the exterior
        // condition u0/u<5 breaks down when u0>>1).  Also scale the
        // allowed change with depth so rays can traverse the large
        // u-range without exhausting the step budget.
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

        // Additional step refinement near photon sphere (u ≈ 0.667 for r = 1.5)
        float u_photon_sphere = 0.667;
        float photon_sphere_proximity = exp(-12.0 * (u - u_photon_sphere) * (u - u_photon_sphere));
        step *= 1.0 - 0.7 * photon_sphere_proximity;

        // Interior mode: refine steps crossing the horizon (u ≈ 1) for stability
        if (interior_mode > 0.5) {
            float horizon_proximity = exp(-8.0 * (u - 1.0) * (u - 1.0));
            step *= 1.0 - 0.6 * horizon_proximity;
            // Floor: prevent step from becoming too tiny
            step = max(step, 0.5 * M_PI / float(NSTEPS));
        }

        old_u = u;

        {{#light_travel_time}}
        {{#gravitational_time_dilation}}
        dt = sqrt(max(du*du + u*u*(1.0-u), 0.0001))/max(abs(u*u*(1.0-u)), 0.0001)*step;
        {{/gravitational_time_dilation}}
        {{/light_travel_time}}

        integrate_geodesic_step(u, du, step, spin_alignment);

        // u < 0 is non-physical.  In exterior it means the ray escaped;
        // in interior it is a numerical artefact → treat as black.
        if (u < 0.0) {
            if (interior_mode > 0.5) shadow_capture = true;
            break;
        }
        // Interior mode: trace past horizon, stop only at singularity
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
        {{/kerr_fast_mode}}

        {{#kerr_full_core}}
        // Realtime Kerr: Binet equation with frame-drag approximation
        // Same geodesic integration as fast mode, with full Kerr disk kinematics
        step = MAX_REVOLUTIONS * 2.0*M_PI / float(NSTEPS);

        // adaptive step size based on rate of change
        float max_rel_u_change_fc = max(1.0-log(max(u, 0.0001)), 0.05)*10.0 / float(NSTEPS);
        if (interior_mode > 0.5) {
            max_rel_u_change_fc *= 3.0 + 2.0 * max(u - 1.0, 0.0);
            if (abs(du) > abs(max_rel_u_change_fc*u) / step) {
                step = max_rel_u_change_fc*u/abs(du);
            }
        } else {
            if ((du > 0.0 || (du0 < 0.0 && u0/u < 5.0)) && abs(du) > abs(max_rel_u_change_fc*u) / step) {
                step = max_rel_u_change_fc*u/abs(du);
            }
        }

        // Additional step refinement near photon sphere (u ≈ 0.667 for r = 1.5)
        float u_photon_sphere_fc = 0.667;
        float photon_sphere_proximity_fc = exp(-12.0 * (u - u_photon_sphere_fc) * (u - u_photon_sphere_fc));
        step *= 1.0 - 0.7 * photon_sphere_proximity_fc;

        // Interior mode: refine steps crossing the horizon (u ≈ 1)
        if (interior_mode > 0.5) {
            float horizon_proximity_fc = exp(-8.0 * (u - 1.0) * (u - 1.0));
            step *= 1.0 - 0.6 * horizon_proximity_fc;
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
        // Interior mode: trace past horizon, stop only at singularity
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
        vec3 planar_pos_fc = (cos(phi)*normal_vec + sin(phi)*tangent_vec)/u;

        float drag_u_fc = clamp(u, 0.0, 1.15);
        float frame_drag_step_fc = bh_rotation_enabled * bh_spin * bh_spin_strength *
            spin_alignment * step * 0.85 * drag_u_fc*drag_u_fc*drag_u_fc;
        frame_drag_phase += frame_drag_step_fc;

        pos = rotate_about_z(planar_pos_fc, frame_drag_phase);

        {{#light_travel_time}}
        {{^gravitational_time_dilation}}
        dt = length(pos - old_pos);
        {{/gravitational_time_dilation}}
        {{/light_travel_time}}
        {{/kerr_full_core}}

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
        if (old_pos.z * pos.z < 0.0) {
            // crossed plane z=0

            float acc_isec_t = -old_pos.z / ray.z;
            if (acc_isec_t < solid_isec_t) {
                vec3 isec = old_pos + ray*acc_isec_t;

                float r = length(isec);
                if (r > ACCRETION_MIN_R && r < ACCRETION_MIN_R + ACCRETION_WIDTH) {
                    float angle = atan(isec.x, isec.y);

                    float temperature = accretion_temperature(r) * gravitational_shift(r);
                    float turbulence = accretion_emissivity(r, angle, t);
                    float inner_glow = exp(-8.0 * (r - ACCRETION_MIN_R));
                    float accretion_intensity = ACCRETION_BRIGHTNESS * accretion_flux_profile(r) *
                        turbulence * (1.0 + 0.7*inner_glow);

                    vec3 accretion_v;
                    {{#kerr_fast_mode}}
                    accretion_v = vec3(-isec.y, isec.x, 0.0) / (r * sqrt(2.0*(r-1.0)));
                    {{/kerr_fast_mode}}
                    {{#kerr_full_core}}
                    float rg_r = max(2.0*r, 1.0002); // convert r_s units to M units
                    float a_M = bh_rotation_enabled * bh_spin;
                    float omega_M = 1.0 / (pow(rg_r, 1.5) + a_M);
                    float v_phi = clamp(rg_r * omega_M, -0.995, 0.995);
                    vec2 xy = vec2(isec.x, isec.y);
                    float cyl_r = max(length(xy), 1e-4);
                    vec3 e_phi = vec3(-xy.y/cyl_r, xy.x/cyl_r, 0.0);
                    accretion_v = e_phi * v_phi;
                    {{/kerr_full_core}}

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
        {{/accretion_thin_disk}}

        {{#accretion_thick_torus}}
        {
            // ADAF/RIAF volumetric emission: optically thin radiative transfer
            // dI/ds = j - alpha*I  (emission minus absorption)
            float torus_j = torus_local_emissivity(pos);
            if (torus_j > 0.001 && vol_transmittance > 0.005) {
                float path_len = length(pos - old_pos);
                float r3d = max(length(pos), 1.001);
                float cyl_r_t = max(length(pos.xy), 1e-4);
                float angle_t = atan(pos.x, pos.y);

                float temperature_t = torus_temperature(r3d) * gravitational_shift(r3d);
                float turbulence_t = accretion_turbulence(cyl_r_t, angle_t, t);

                // Emission coefficient: optically thin ADAF has low emissivity
                float j_eff = ACCRETION_BRIGHTNESS * 0.05 * torus_j * turbulence_t;
                // Absorption: user-configurable opacity (default low for optically thin)
                float alpha_abs = torus_opacity * torus_j;
                float tau_step = alpha_abs * path_len;
                float step_T = exp(-tau_step);

                // Sub-Keplerian gas velocity (ADAF: ~50% of Keplerian)
                vec3 accretion_v_t;
                {{#kerr_fast_mode}}
                float v_kep_t = 1.0 / sqrt(2.0 * max(cyl_r_t - 1.0, 0.01));
                float v_sub_t = clamp(0.5 * v_kep_t, 0.0, 0.95);
                accretion_v_t = vec3(-pos.y, pos.x, 0.0) / cyl_r_t * v_sub_t;
                {{/kerr_fast_mode}}
                {{#kerr_full_core}}
                float rg_cyl_t = max(2.0 * cyl_r_t, 1.0002);
                float a_M_t = bh_rotation_enabled * bh_spin;
                float omega_sub_t = 0.5 / (pow(rg_cyl_t, 1.5) + a_M_t);
                float v_phi_t = clamp(rg_cyl_t * omega_sub_t, -0.95, 0.95);
                vec3 e_phi_t = vec3(-pos.y, pos.x, 0.0) / cyl_r_t;
                accretion_v_t = e_phi_t * v_phi_t;
                {{/kerr_full_core}}

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
            float slim_j = slim_disk_local_emissivity(pos);
            if (slim_j > 0.001 && vol_transmittance > 0.005) {
                float path_len_s = length(pos - old_pos);
                float cyl_r_s = max(length(pos.xy), 1e-4);
                float r3d_s = max(length(pos), 1.001);
                float angle_s = atan(pos.x, pos.y);

                float temperature_s = slim_disk_temperature(cyl_r_s) * gravitational_shift(r3d_s);
                float turbulence_s = accretion_turbulence(cyl_r_s, angle_s, t);

                // Emission coefficient: slim disk is bright (super-Eddington luminosity)
                float j_eff_s = ACCRETION_BRIGHTNESS * 0.9 * slim_j * turbulence_s;
                // User-configurable absorption: higher = more opaque surface-like
                float alpha_abs_s = slim_opacity * slim_j;
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
                {{#kerr_full_core}}
                float rg_cyl_s = max(2.0 * cyl_r_s, 1.0002);
                float a_M_s = bh_rotation_enabled * bh_spin;
                float omega_s = 1.0 / (pow(rg_cyl_s, 1.5) + a_M_s);
                float v_phi_s = clamp(rg_cyl_s * omega_s, -0.95, 0.95);
                vec3 e_phi_s = vec3(-pos.y, pos.x, 0.0) / cyl_r_s;
                accretion_v_s = e_phi_s * v_phi_s;
                {{/kerr_full_core}}

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

                float j_jet = jet_emissivity(pos, sign_z);
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
                    float g_shift_j = gravitational_shift(r_jet);
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
                    float g_shift_j = gravitational_shift(r_jet_p);
                    float sigma = magnetization(z_abs);

                    // Effective temperature with synchrotron aging
                    float jet_T = jet_effective_temperature(z_abs, r_jet_p, sigma) * g_shift_j;

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

        vec4 star_color = texture2D(star_texture, tex_coord);
        if (star_color.r > 0.0) {
            t_coord = (STAR_MIN_TEMPERATURE +
                (STAR_MAX_TEMPERATURE-STAR_MIN_TEMPERATURE) * star_color.g)
                 / ray_doppler_factor;

            color += BLACK_BODY_COLOR(t_coord) * star_color.r * STAR_BRIGHTNESS * look_star_gain * vol_transmittance * interior_boost;
        }

        color += galaxy_color(tex_coord, ray_doppler_factor) * GALAXY_BRIGHTNESS * look_galaxy_gain * vol_transmittance * interior_boost;
    }

    return color*ray_intensity;
}
