        {{#accretion_slim_disk}}
        {
            // Slim disk volumetric emission: optically thick radiative transfer
            // Super-Eddington flow: high absorption → surface-like rendering
            {{#grmhd_enabled}}
            // GRMHD: height modulation from MRI + radiation-driven warps
            float _cr_s = max(length(pos.xy), 1e-4);
            float _a_s = equatorial_azimuth(pos.xy);
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
                float angle_s = equatorial_azimuth(pos.xy);

                {{#grmhd_enabled}}
                    // GRMHD-inspired slim disk: preserve the base super-Eddington
                    // morphology while adding responsive semi-analytic corrections.
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
                accretion_v_s = disk_rotation_sign() *
                    vec3(-pos.y, pos.x, 0.0) / cyl_r_s * clamp(v_kep_s, 0.0, 0.95);
                {{/kerr_fast_mode}}
                {{#kerr_inspired_velocity}}
                float rg_cyl_s = max(2.0 * cyl_r_s, 1.0002);
                float a_M_s = abs(bh_rotation_enabled * bh_spin);
                float omega_s = disk_rotation_sign() / (pow(rg_cyl_s, 1.5) + a_M_s);
                float v_phi_s = clamp(rg_cyl_s * omega_s, -0.95, 0.95);
                vec3 e_phi_s = vec3(-pos.y, pos.x, 0.0) / cyl_r_s;
                accretion_v_s = e_phi_s * v_phi_s;
                {{/kerr_inspired_velocity}}

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
