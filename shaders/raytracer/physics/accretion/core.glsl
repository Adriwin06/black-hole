// Role: Accretion disk emissivity and temperature models — thin Shakura-Sunyaev disk,
//       ADAF/RIAF thick torus, and super-Eddington slim disk. Also computes the
//       gravitational redshift factor and the planet irradiation temperature.

float loopable_turbulence_time(float t) {
    // Always wrap to prevent float32 precision loss in fbm/value_noise.
    // Without wrapping, large t values cause fract()/floor() to lose
    // sub-grid-cell resolution in the highest FBM octaves, flattening
    // fine turbulence detail over time (most visible after ~1 hour).
    // When the explicit loop is off, use a long prime period (10007 s
    // ≈ 2.78 h) so the wrap is imperceptible; when on, use the user
    // period for seamless video loops.
    float period = (turbulence_loop_enabled < 0.5)
        ? 10007.0
        : max(turbulence_loop_seconds, 1e-4);
    float wrapped = mod(t, period);
    if (wrapped < 0.0) wrapped += period;
    return wrapped;
}

float disk_rotation_sign() {
    // Keep the existing prograde default when spin is zero/off, but flip the
    // co-rotating flow for explicitly negative a/M.
    return (bh_rotation_enabled > 0.5 && bh_spin < 0.0) ? -1.0 : 1.0;
}

float equatorial_azimuth(vec2 xy) {
    // GLSL's two-argument atan is atan(y, x); swapping the arguments reverses
    // the phase advection direction.
    return atan(xy.y, xy.x);
}

float accretion_turbulence(float radius, float angle, float t) {
    float orbit_phase = angle - disk_rotation_sign() * 0.45*t / pow(max(radius, 1.001), 1.5);
    vec2 orbit_unit = vec2(cos(orbit_phase), sin(orbit_phase));

    // Use periodic angular coordinates to avoid seam artifacts at angle wrap.
    float large_scale = fbm(vec2(
        radius*2.5 + orbit_unit.x*4.5,
        orbit_unit.y*4.5 + t*0.08
    ));
    float small_scale = fbm(vec2(
        radius*12.0 + orbit_unit.x*15.0 + 2.5*large_scale,
        orbit_unit.y*15.0 - t*0.18
    ));
    // High-frequency chaos layer, moving faster and concentrated near the middle (ISCO)
    float micro_scale = fbm(vec2(
        radius*25.0 + orbit_unit.x*35.0 - small_scale*4.0,
        orbit_unit.y*35.0 + t*0.55
    ));
    
    // Distort the regular spiral using small_scale noise so it looks chaotic instead of uniform stripes
    float swirl_phase = 12.0*orbit_phase + 8.0*log(max(radius, 1.001)) + small_scale * 4.0;
    float swirl = sin(swirl_phase);
    
    float isco_factor = exp(-0.8 * max(radius - ACCRETION_MIN_R, 0.0));
    float filaments = 0.6 + 0.4*swirl;
    
    // Mix scales, giving micro_scale more weight near the center
    float base_plasma = mix(large_scale, small_scale, 0.5);
    float plasma = mix(base_plasma, micro_scale, 0.3 + 0.5 * isco_factor);

    return max((0.4 + 1.25*plasma) * (0.8 + 0.2*filaments), 0.02);
}

float accretion_emissivity(float radius, float angle, float t) {
    float r_norm = (radius - ACCRETION_MIN_R) / ACCRETION_WIDTH;
    float edge_fade = smoothstep(0.02, 0.18, r_norm) *
        (1.0 - smoothstep(0.78, 1.0, r_norm));
    return accretion_turbulence(radius, angle, t) * edge_fade;
}

// --- Disk Self-Irradiation (Returning Radiation) ---
// Heuristic Cunningham-inspired inner-disk brightening. This is NOT a
// ray-traced returning-radiation calculation or a tabulated Cunningham transfer
// function; it is a spin-scaled enhancement localized near the ISCO so the
// renderer can mimic the qualitative extra heating of strongly lensed inner
// disk emission.
float accretion_returning_radiation_enhancement(float radius) {
{{#disk_self_irradiation_enabled}}
    float r_norm = max(radius / ACCRETION_MIN_R, 1.0001);
    // Higher spin moves the ISCO inward, so we allow a stronger heuristic boost.
    float spin_a = bh_rotation_enabled * bh_spin;
    float peak_enhancement = 0.2 + 1.2 * abs(spin_a);
    // Rapid outward decay keeps the effect concentrated near the inner disk.
    float enhancement = peak_enhancement * pow(1.0 / r_norm, 3.5);
    return 1.0 + enhancement;
{{/disk_self_irradiation_enabled}}
{{^disk_self_irradiation_enabled}}
    return 1.0;
{{/disk_self_irradiation_enabled}}
}

float accretion_flux_profile(float radius) {
    float x = max(radius / ACCRETION_MIN_R, 1.0001);
    float inner_edge = max(1.0 - sqrt(1.0 / x), 0.0);
    float flux = inner_edge / (x*x*x);
    // Multiply by returning radiation enhancement to capture self-irradiation heating
    return flux * 18.0 * accretion_returning_radiation_enhancement(radius);
}

float accretion_temperature(float radius) {
    // Normalize Shakura-Sunyaev profile so disk_temperature corresponds
    // to the peak effective temperature (at x = 49/36).
    const float SS_PEAK_NORMALIZATION = 2.04910267;
    float x = max(radius / ACCRETION_MIN_R, 1.0001);
    float inner_edge = max(1.0 - sqrt(1.0 / x), 0.02);
    
    float t_base = disk_temperature * SS_PEAK_NORMALIZATION *
        pow(1.0 / x, 0.75) * pow(inner_edge, 0.25);
        
    // Stefan-Boltzmann law: Flux prop T^4. Thus, T = T_base * (F / F_base)^0.25
    return t_base * pow(accretion_returning_radiation_enhancement(radius), 0.25);
}

// Observer metric factor shared by the Schwarzschild redshift helpers.
// Uses the static Schwarzschild factor because the observer's own
// kinematic Doppler is applied separately via cam_vel aberration.
float _observer_metric_factor() {
    float r_obs = max(length(cam_pos), 0.05);
    return max(abs(1.0 - 1.0 / r_obs), 0.001);
}

float gravitational_shift_static(float emission_radius) {
    // Static Schwarzschild redshift only. Moving-source kinematic Doppler
    // (including transverse time dilation via gamma) is applied separately in
    // trace_ray(), which avoids double counting the emitter gamma factor.
    float r_emit = max(emission_radius, 1.0001);
    float emission_term = max(abs(1.0 - 1.0 / r_emit), 0.0001);
    return sqrt(emission_term / _observer_metric_factor());
}
