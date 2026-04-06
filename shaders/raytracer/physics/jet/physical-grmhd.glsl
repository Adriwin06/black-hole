{{#grmhd_enabled}}
// ═══════════════════════════════════════════════════════════════════
// GRMHD JET ENHANCEMENTS
// ═══════════════════════════════════════════════════════════════════
// Additional physics for the GRMHD-inspired jet mode:
//  - Blandford-Znajek power scaling: P_BZ ∝ Φ² Ω_H² f(Ω_H)
//  - MAD magnetic flux saturation: Φ_BH ~ 50 √(Ṁ r_g² c)
//  - Current-driven kink instability (CDI) modulation
//  - Synchrotron cooling spectral break
//  - Enhanced spine/sheath contrast from force-free σ(r,θ) profile

// --- Blandford-Znajek jet power factor ---
// P_BZ = (κ/4πc) Φ_BH² Ω_H² f(Ω_H) where Ω_H = a/(2r_+)
// In MAD state, Φ_BH saturates → P_BZ ∝ a² Ṁ (Tchekhovskoy+2011)
float grmhd_bz_power_factor() {
    float a = clamp(bh_spin * bh_rotation_enabled, -0.99, 0.99);
    float r_plus = 0.5 + 0.5 * sqrt(max(1.0 - a * a, 0.0));
    float omega_H = a / (2.0 * r_plus);
    // P_BZ ∝ Φ² Ω_H² ≈ a² for MAD; for SANE, Φ scales weaker
    float phi_factor = mix(0.3, 1.0, grmhd_mad_flux); // MAD has higher Φ
    return phi_factor * omega_H * omega_H * 40.0 + 0.1;
}

// --- Current-driven kink instability (CDI) ---
// The helical magnetic field in the jet is subject to CDI for
// wavelengths > jet radius. Creates m = 1 helical distortions.
// (Bromberg & Tchekhovskoy 2016)
float grmhd_kink_instability(float z, float cyl_r, float jet_r) {
    float growth_length = max(jet_r * 4.0, 5.0);
    float kink_amp = smoothstep(growth_length, growth_length * 2.0, z);
    // Helical m = 1 distortion
    float phase = z * 0.8 / max(jet_r, 0.5);
    float kink = 1.0 + 0.25 * kink_amp * sin(phase) * grmhd_mad_flux;
    return kink;
}

// --- Synchrotron cooling break ---
// Electrons lose energy as t_cool ∝ 1/(B² γ_e) → spectral steepening
// at high frequencies. Creates a characteristic "cooling break" in the
// spectrum. Visual effect: jet becomes progressively redder with distance.
float grmhd_cooling_break_temp_modifier(float z, float sigma) {
    // Cooling is faster in high-B (low-σ) regions near the base
    float cool_rate = 1.0 + 0.5 / max(sigma, 0.1);
    // Temperature decreases with distance due to radiative losses
    float cooling = exp(-0.05 * cool_rate * max(z - 3.0, 0.0));
    return max(cooling, 0.3);
}

// --- Enhanced GRMHD jet emissivity ---
// Replaces the standard jet_emissivity when GRMHD mode is active.
// Adds BZ power scaling, kink instability, and density fluctuations.
float grmhd_jet_emissivity(vec3 p, float sign_z) {
    float z = p.z * sign_z;
    if (z < 0.8) return 0.0;

    float cyl_r = length(p.xy);
    float r3d = length(p);
    float jet_r = jet_boundary_radius(z);

    // BZ power factor modulates overall jet luminosity
    float bz = grmhd_bz_power_factor();

    if (cyl_r > jet_r * 1.15) {
        if (z < 5.0 && r3d < 5.0) {
            return corona_base_emission(p, z, cyl_r) * bz;
        }
        return 0.0;
    }

    float sigma = magnetization(z);
    float profile = jet_transverse_profile(cyl_r, jet_r, sigma);
    float z_onset = smoothstep(0.8, 2.5, z);
    float z_decay = pow(max(z, 1.0), -2.2);
    float z_cutoff = 1.0 - smoothstep(jet_length * 0.7, jet_length, z);
    float knots = reconfinement_knots(z, jet_r);

    // GRMHD enhancements
    float kink = grmhd_kink_instability(z, cyl_r, jet_r);
    float nonthermal = grmhd_nonthermal_boost(grmhd_electron_kappa);

    // GRMHD turbulent density structure in the jet
    float angle_j = equatorial_azimuth(p.xy);
    float jet_turb = grmhd_mri_turbulence(r3d, angle_j, loopable_turbulence_time(time));

    float j = profile * z_onset * z_decay * z_cutoff * knots
            * kink * nonthermal * jet_turb * bz;

    j += corona_base_emission(p, z, cyl_r) * bz;
    return j;
}

// --- GRMHD jet effective temperature ---
// Incorporates synchrotron cooling break and non-thermal broadening.
float grmhd_jet_temperature(float z, float r3d, float sigma) {
    float T_base = jet_effective_temperature(z, r3d, sigma);
    float cooling = grmhd_cooling_break_temp_modifier(z, sigma);
    return T_base * cooling;
}
{{/grmhd_enabled}}
{{/jet_physical}}
{{/jet_enabled}}
