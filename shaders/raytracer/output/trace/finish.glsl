        {{#light_travel_time}}
        t -= dt;
        {{/light_travel_time}}

        if (solid_isec_t <= 1.0) u = 100.0; // break (exceeds both exterior and interior thresholds)

        if (!use_kerr) {
        {{#kerr_fast_mode}}
        if (interior_mode < 0.5) {
            float capture_u = 1.0 + photon_spin_lensing_scale *
                bh_rotation_enabled * bh_spin * bh_spin_strength *
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
        {{#kerr_inspired_mode}}
        // Kerr-inspired mode: same Binet photon solver, same interior capture
        if (interior_mode > 0.5 && u >= 20.0) {
            shadow_capture = true;
            break;
        }
        {{/kerr_inspired_mode}}
        {{#kerr_full_geodesic}}
        // Full Kerr geodesic fallback (Binet): interior capture
        if (interior_mode > 0.5 && u >= 20.0) {
            shadow_capture = true;
            break;
        }
        {{/kerr_full_geodesic}}
        } // end if (!use_kerr)
    }

        trace_finished = shadow_capture || u < 1.0 || use_kerr || !allow_extended_trace || u > 20.0;
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
