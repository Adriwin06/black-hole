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
    // User-configurable power law index (default ~2.5, physically 2-4)
    // Using r0 as the reference radius for the peak emissivity region
    float falloff = max(torus_radial_falloff, 0.5);
    float radial;
    if (cyl_r < r0) {
        // Inside torus center: emissivity rises but not as steeply
        radial = pow(r0 / max(cyl_r, 1.0), falloff * 0.6);
    } else {
        // Outside torus center: standard power-law decay
        radial = pow(r0 / cyl_r, falloff);
    }
    // Normalize so peak is at r0 with a smooth profile
    float torus_profile = exp(-1.2 * pow((cyl_r - r0) / max(r0 * 0.6, 1.0), 2.0));
    radial = mix(radial, radial * torus_profile, 0.5);

    // Smooth cutoff at event horizon
    float inner = smoothstep(0.9, 1.5, cyl_r);

    // Outer edge falloff using configurable outer radius
    float outer_r = max(r0 * torus_outer_radius, ACCRETION_MIN_R + ACCRETION_WIDTH);
    float outer = 1.0 - smoothstep(outer_r * 0.7, outer_r, cyl_r);

    return vert * radial * inner * outer;
}

float torus_temperature(float r) {
    float r0 = max(torus_r0, 1.5);
    // ADAF electron temperature: T_e ~ r^(-0.5) for two-temperature ADAF
    // (Narayan & Yi 1995, Esin et al. 1997)
    return disk_temperature * pow(r0 / max(r, 1.0), 0.5);
}
{{/accretion_thick_torus}}
