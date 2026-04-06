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
                    // Detailed mode: more elaborate GRMHD-inspired emissivity,
                    // still coloured via an effective-temperature proxy.
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
