// Role: Accretion disk emissivity and temperature models — thin Shakura-Sunyaev disk,
//       ADAF/RIAF thick torus, and super-Eddington slim disk. Also computes the
//       gravitational redshift factor and the planet irradiation temperature.

float accretion_turbulence(float radius, float angle, float t) {
    float orbit_phase = angle - 0.45*t / pow(max(radius, 1.001), 1.5);
    vec2 orbit_unit = vec2(cos(orbit_phase), sin(orbit_phase));
    float swirl = sin(18.0*orbit_phase + 10.0*log(max(radius, 1.001)));

    // Use periodic angular coordinates to avoid seam artifacts at angle wrap.
    float large_scale = fbm(vec2(
        radius*1.5 + orbit_unit.x*2.7,
        orbit_unit.y*2.7 + t*0.05
    ));
    float small_scale = fbm(vec2(
        radius*7.2 + orbit_unit.x*9.0 + 1.5*large_scale,
        orbit_unit.y*9.0 - t*0.11
    ));
    float filaments = 0.6 + 0.4*swirl;
    float plasma = mix(large_scale, small_scale, 0.6);

    return max((0.45 + 1.1*plasma) * (0.8 + 0.2*filaments), 0.02);
}

float accretion_emissivity(float radius, float angle, float t) {
    float r_norm = (radius - ACCRETION_MIN_R) / ACCRETION_WIDTH;
    float edge_fade = smoothstep(0.02, 0.18, r_norm) *
        (1.0 - smoothstep(0.78, 1.0, r_norm));
    return accretion_turbulence(radius, angle, t) * edge_fade;
}

float accretion_flux_profile(float radius) {
    float x = max(radius / ACCRETION_MIN_R, 1.0001);
    float inner_edge = max(1.0 - sqrt(1.0 / x), 0.0);
    float flux = inner_edge / (x*x*x);
    return flux * 18.0;
}

float accretion_temperature(float radius) {
    // Normalize Shakura-Sunyaev profile so disk_temperature corresponds
    // to the peak effective temperature (at x = 49/36).
    const float SS_PEAK_NORMALIZATION = 2.04910267;
    float x = max(radius / ACCRETION_MIN_R, 1.0001);
    float inner_edge = max(1.0 - sqrt(1.0 / x), 0.02);
    return disk_temperature * SS_PEAK_NORMALIZATION *
        pow(1.0 / x, 0.75) * pow(inner_edge, 0.25);
}

float gravitational_shift(float emission_radius) {
    float observer_term = max(1.0 - 1.0 / max(length(cam_pos), 1.0001), 0.0001);
    float emission_term = max(1.0 - 1.0 / max(emission_radius, 1.0001), 0.0001);
    return sqrt(emission_term / observer_term);
}

{{#accretion_thick_torus}}
// ADAF/RIAF thick torus: geometrically thick, optically thin accretion flow
// Models low-luminosity AGN like M87* and Sgr A* (EHT targets)
float torus_local_emissivity(vec3 p) {
    float cyl_r = length(p.xy);
    float r0 = max(torus_r0, 1.5);
    float h = max(cyl_r * torus_h_ratio, 0.01);
    float z_norm = abs(p.z) / h;

    if (z_norm > 3.0 || cyl_r < 0.9) return 0.0;

    // Gaussian vertical profile (hydrostatic equilibrium)
    float vert = exp(-0.5 * z_norm * z_norm);

    // Radial emissivity: bremsstrahlung j ~ n^2 T^(1/2)
    // Self-similar ADAF (Narayan & Yi 1994): n ~ r^(-3/2), T ~ r^(-1)
    // j_ff ~ r^(-3) * r^(-1/2) = r^(-3.5)
    float radial = pow(r0 / max(cyl_r, 1.0), 3.5);

    // Smooth cutoff at event horizon
    float inner = smoothstep(0.9, 1.5, cyl_r);

    // Outer edge falloff
    float outer_r = max(r0 * 3.5, ACCRETION_MIN_R + ACCRETION_WIDTH);
    float outer = 1.0 - smoothstep(outer_r * 0.75, outer_r, cyl_r);

    return vert * radial * inner * outer;
}

float torus_temperature(float r) {
    float r0 = max(torus_r0, 1.5);
    // ADAF electron temperature: T_e ~ r^(-0.5) for two-temperature ADAF
    // (Narayan & Yi 1995, Esin et al. 1997)
    return disk_temperature * pow(r0 / max(r, 1.0), 0.5);
}
{{/accretion_thick_torus}}

{{#accretion_slim_disk}}
// Slim disk: super-Eddington accretion, extends inside ISCO
// Radiation-pressure supported, geometrically thicker than thin disk
float slim_disk_height(float cyl_r) {
    float base_h_ratio = 0.12;
    // Puffs up near and inside ISCO due to radiation pressure
    float isco_proximity = exp(-1.5 * max(cyl_r - ACCRETION_MIN_R, 0.0));
    return max(cyl_r * base_h_ratio * (1.0 + 2.5 * isco_proximity), 0.01);
}

float slim_disk_local_emissivity(vec3 p) {
    float cyl_r = length(p.xy);
    float h = slim_disk_height(cyl_r);
    float z_norm = abs(p.z) / h;

    if (z_norm > 3.0 || cyl_r < 0.9) return 0.0;

    float vert = exp(-0.5 * z_norm * z_norm);

    // Extends inside ISCO with plunging-region emission
    float radial;
    if (cyl_r >= ACCRETION_MIN_R) {
        radial = accretion_flux_profile(cyl_r);
    } else {
        // Plunging region: conserved specific energy, decreasing efficiency
        float f = cyl_r / max(ACCRETION_MIN_R, 1.0);
        radial = accretion_flux_profile(ACCRETION_MIN_R) * pow(f, 2.5);
    }

    float inner = smoothstep(0.9, 1.3, cyl_r);
    float outer_r = ACCRETION_MIN_R + ACCRETION_WIDTH;
    float outer = 1.0 - smoothstep(outer_r * 0.8, outer_r, cyl_r);

    return vert * radial * inner * outer;
}

float slim_disk_temperature(float cyl_r) {
    if (cyl_r >= ACCRETION_MIN_R) {
        return accretion_temperature(cyl_r);
    } else {
        // Inside ISCO: advection-dominated, T rises as r^(-0.5)
        float t_isco = accretion_temperature(ACCRETION_MIN_R);
        return t_isco * pow(ACCRETION_MIN_R / max(cyl_r, 1.0), 0.5);
    }
}
{{/accretion_slim_disk}}

float planet_irradiation_temperature() {
    {{#accretion_thin_disk}}
    float r1 = ACCRETION_MIN_R * 1.8;
    float r2 = ACCRETION_MIN_R * 2.8;
    float r3 = ACCRETION_MIN_R * 4.2;

    float w1 = accretion_flux_profile(r1);
    float w2 = accretion_flux_profile(r2);
    float w3 = accretion_flux_profile(r3);
    float wsum = max(w1 + w2 + w3, 1e-4);

    float t1 = accretion_temperature(r1);
    float t2 = accretion_temperature(r2);
    float t3 = accretion_temperature(r3);
    t1 *= gravitational_shift(r1);
    t2 *= gravitational_shift(r2);
    t3 *= gravitational_shift(r3);
    return (w1*t1 + w2*t2 + w3*t3) / wsum;
    {{/accretion_thin_disk}}

    {{#accretion_thick_torus}}
    float r_t = max(torus_r0 * 1.25, ACCRETION_MIN_R + 0.5);
    return torus_temperature(r_t) * gravitational_shift(r_t);
    {{/accretion_thick_torus}}

    {{#accretion_slim_disk}}
    float r_s = ACCRETION_MIN_R * 1.6;
    return slim_disk_temperature(r_s) * gravitational_shift(r_s);
    {{/accretion_slim_disk}}

    return disk_temperature;
}
