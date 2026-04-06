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
            float _a_t = equatorial_azimuth(pos.xy);
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
                float angle_t = equatorial_azimuth(pos.xy);

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
                accretion_v_t = disk_rotation_sign() *
                    vec3(-pos.y, pos.x, 0.0) / cyl_r_t * v_sub_t;
                {{/kerr_fast_mode}}
                {{#kerr_inspired_velocity}}
                float rg_cyl_t = max(2.0 * cyl_r_t, 1.0002);
                float a_M_t = abs(bh_rotation_enabled * bh_spin);
                float omega_sub_t = 0.5 * disk_rotation_sign() / (pow(rg_cyl_t, 1.5) + a_M_t);
                float v_phi_t = clamp(rg_cyl_t * omega_sub_t, -0.95, 0.95);
                vec3 e_phi_t = vec3(-pos.y, pos.x, 0.0) / cyl_r_t;
                accretion_v_t = e_phi_t * v_phi_t;
                {{/kerr_inspired_velocity}}

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
