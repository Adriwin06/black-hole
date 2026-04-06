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
                    float angle = equatorial_azimuth(isec.xy);

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
                    // Apply only the static Schwarzschild factor here; the
                    // separate transfer_factor below already carries the
                    // emitter gamma and line-of-sight Doppler contribution.
                    float temperature = grmhd_electron_temperature(gas_temp_g, r, 0.0, disk_h_g, turb_t, angle) * Te_corr * gravitational_shift_static(r);
                    
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
                    float temperature = accretion_temperature(r) * gravitational_shift_static(r);
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
                    accretion_v = disk_rotation_sign() *
                        vec3(-isec.y, isec.x, 0.0) / (r * sqrt(2.0*(r-1.0)));
                    {{/kerr_fast_mode}}
                    {{#kerr_inspired_velocity}}
                    float rg_r = max(2.0*r, 1.0002); // convert r_s units to M units
                    float a_M = abs(bh_rotation_enabled * bh_spin);
                    float omega_M = disk_rotation_sign() / (pow(rg_r, 1.5) + a_M);
                    float v_phi = clamp(rg_r * omega_M, -0.995, 0.995);
                    vec2 xy = vec2(isec.x, isec.y);
                    float cyl_r = max(length(xy), 1e-4);
                    vec3 e_phi = vec3(-xy.y/cyl_r, xy.x/cyl_r, 0.0);
                    accretion_v = e_phi * v_phi;
                    {{/kerr_inspired_velocity}}

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
