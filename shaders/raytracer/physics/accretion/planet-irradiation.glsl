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
    // Use the static Schwarzschild factor here. The detailed, angle-dependent
    // disk orbital Doppler seen by the planet is not resolved in this helper.
    t1 *= gravitational_shift_static(r1);
    t2 *= gravitational_shift_static(r2);
    t3 *= gravitational_shift_static(r3);
    return (w1*t1 + w2*t2 + w3*t3) / wsum;
    {{/accretion_thin_disk}}

    {{#accretion_thick_torus}}
    float r_t = max(torus_r0 * 1.25, ACCRETION_MIN_R + 0.5);
    return torus_temperature(r_t) * gravitational_shift_static(r_t);
    {{/accretion_thick_torus}}

    {{#accretion_slim_disk}}
    float r_s = ACCRETION_MIN_R * 1.6;
    return slim_disk_temperature(r_s) * gravitational_shift_static(r_s);
    {{/accretion_slim_disk}}

    return disk_temperature;
}
