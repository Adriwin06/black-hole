{{#grmhd_enabled}}
// ═══════════════════════════════════════════════════════════════════
// GRMHD-INSPIRED PHYSICS MODULE
// ═══════════════════════════════════════════════════════════════════
// Semi-analytic models inspired by GRMHD simulation results from:
//  - Event Horizon Telescope Collaboration et al. (2019, Paper V)
//  - Mościbrodzka, Falcke & Shiokawa (2016): two-temperature GRMHD
//  - Narayan, Sądowski et al. (2012): GRMHD accretion flows
//  - Fishbone & Moncrief (1976): equilibrium torus initial data
//  - Mahadevan, Narayan & Yi (1996): synchrotron fitting functions
//  - Howes (2010): turbulent electron heating prescription
//  - Pandya et al. (2016): synchrotron emissivities & absorptivities
//  - Tchekhovskoy, Narayan & McKinney (2011): MAD simulations

// --- Plasma beta profile ---
// β = P_gas / P_mag varies with position in GRMHD simulations:
//  Disk midplane: β >> 1 (gas-pressure dominated, MRI-turbulent)
//  Disk corona:   β ~ 1  (magnetic/gas equipartition)
//  Jet funnel:    β << 1 (magnetically dominated Poynting flux)
//  MAD state:     β lower everywhere (arrested magnetic flux)
float grmhd_plasma_beta(float cyl_r, float z, float disk_h) {
    float z_norm = abs(z) / max(disk_h, 0.01);
    // Midplane β from the user parameter
    float beta_mid = grmhd_magnetic_beta;
    // β decreases away from midplane (magnetic pressure rises)
    float beta = beta_mid * exp(-1.5 * z_norm * z_norm);
    // MAD state: β reduced by magnetic flux saturation
    beta *= mix(1.0, 0.08, grmhd_mad_flux);
    // β increases outward (field falls faster than gas at large r)
    beta *= pow(max(cyl_r / 4.0, 0.5), 0.5);
    return max(beta, 0.005);
}

// --- Two-temperature electron temperature (EHT R_high prescription) ---
// The physical prescription is usually written for R = T_i / T_e:
//   R(β) = R_high × β² / (1 + β²) + R_low / (1 + β²)
// where R_low = 1 (electrons ≈ ions at equipartition when β → 0).
// This naturally creates:
//   • Hot electrons in the jet funnel (low β → R ≈ 1 → T_e ≈ T_i)
//   • Cool electrons in the disk body (high β → R ≈ R_high → T_e ≈ T_i/R_high)
// When R_high is high, emission becomes concentrated near the jet wall
// and polar regions, producing the crescent morphology seen in EHT images.
// The helper below returns a renderer-side weight 1/(1+R), not the literal T_e/T_i ratio.
float grmhd_electron_temp_ratio(float beta) {
    float R = grmhd_r_high * beta * beta / (1.0 + beta * beta)
            + 1.0 / (1.0 + beta * beta);
    return 1.0 / (1.0 + R);
}

// --- GRMHD-inspired electron temperature ---
// For visual rendering the disk model temperature IS the nominal electron
// temperature (not the GRMHD ion/virial temperature which would be ~10^12 K).
// Applying the full R_high suppression would tank T_e to infrared (~100 K),
// making the disk invisible.  Instead, β modulates temperature mildly (±25%):
// hotter at low-β (magnetically dominated corona/funnel wall), cooler at
// high-β (gas-pressure dominated disk body).  The dominant R_high emission
// morphology is applied separately via grmhd_r_high_emissivity().
float grmhd_electron_temperature(float gas_temp, float cyl_r, float z, float disk_h, float t, float angle) {
    float beta = grmhd_plasma_beta(cyl_r, z, disk_h);
    // Low β → temp_mod ≈ 1.25 (hotter), high β → temp_mod ≈ 0.78 (cooler)
    float temp_mod = 0.75 + 0.5 / (1.0 + 0.15 * beta);
    
    // Smooth time-dependent variation to avoid sharp stripes
    float orbit_phase = angle - disk_rotation_sign() * 0.3 * t / pow(max(cyl_r, 1.0), 1.5);
    vec2 orbit_unit = vec2(cos(orbit_phase), sin(orbit_phase));
    float t_noise = fbm(vec2(
        cyl_r * 2.0 + orbit_unit.x * 2.5,
        orbit_unit.y * 2.5 + t * 0.15
    ));
    // Reduced amplitude to 10% because T^4 causes massive brightness swings
    float t_var = 1.0 + 0.1 * (t_noise - 0.5) * exp(-0.4 * max(cyl_r - ACCRETION_MIN_R, 0.0));
    
    return gas_temp * temp_mod * t_var;
}

// --- R_high emissivity morphology (EHT crescent) ---
// Electron heating efficiency depends on local plasma β via R_high.
// Low-β regions (funnel wall, corona) emit efficiently (factor → 1);
// high-β regions (disk body) are suppressed proportionally to R_high.
// This creates the crescent/ring morphology seen in EHT images.
// Uses sqrt to maintain visual sensitivity across the full R_high range.
float grmhd_r_high_emissivity(float cyl_r, float z, float disk_h) {
    float beta = grmhd_plasma_beta(cyl_r, z, disk_h);
    float Te_ratio = grmhd_electron_temp_ratio(beta);
    // Renderer weight ranges from ~0.5 (β→0) to ~1/(1+R_high) (β→∞)
    // Normalize so low-β regions get factor ≈ 1.0
    return sqrt(Te_ratio * 2.0);
}

// --- Magnetic field strength (dimensionless code units) ---
// B ∝ √(2 ρ T_norm / β) where T_norm = T/T_disk is the dimensionless
// temperature ratio.  This keeps B in the range O(0.1–3) so that
// downstream factors like (1 + c·B) stay well-behaved.
float grmhd_B_field(float density, float temperature, float beta) {
    float T_norm = max(temperature, 1.0) / max(disk_temperature, 1.0);
    return grmhd_magnetic_field_str * sqrt(2.0 * density * T_norm / max(beta, 0.01));
}

// --- Synchrotron spectral correction factor (Mahadevan et al. 1996 inspired) ---
// Synchrotron emissivity scales with B² and electron density.  Uses magnetic
// field strength and density (both O(1) in code units) rather than the
// suppressed T_e ratio, ensuring B-field and density parameters remain
// responsive across their full slider range.
// Returns a modulation factor in the range ~[0.3, 3.0].
float grmhd_synchrotron_correction(float B_norm, float density_norm) {
    float efficiency = B_norm * (0.2 + 0.8 * density_norm);
    return 0.3 + 2.7 * efficiency / (1.0 + efficiency);
}

// --- Synchrotron self-absorption (dimensionless scale factor) ---
// α_ν ∝ n_e / (T_e^(5/2) B) (Kirchhoff's law).  Normalized so that
// typical disk/torus conditions produce O(0.01–1) opacity contributions.
float grmhd_synchrotron_absorption(float density, float T_e_norm, float B_norm) {
    return grmhd_density_scale * density /
        (pow(max(T_e_norm, 0.1), 1.5) * max(B_norm, 0.1) + 0.1);
}

// --- GRMHD density profile ---
// Power-law decay with radius, modified by model-specific normalization.
// For ADAF: ρ ∝ r^(-3/2+s), s ≈ 0.3 (self-similar wind parameter).
// For thin disk: Σ ∝ r^(-3/5) (Shakura-Sunyaev).
float grmhd_density(float r, float base_density) {
    return base_density * pow(4.0 / max(r, 1.0), 1.2) * grmhd_density_scale;
}

// --- MRI turbulence structure (GRMHD-inspired) ---
// MRI generates non-axisymmetric fluctuations with predominantly
// m = 1, 2 azimuthal modes (spiral arms). In MAD state, a dominant
// m = 1 mode arises from the arrested flux tube.
// Density PDF is log-normal with σ_ln(ρ) ≈ 0.5–1.0.
float grmhd_mri_turbulence(float r, float angle, float t) {
    float base = accretion_turbulence(r, angle, t);
    float orbital_freq = 1.0 / pow(max(r, 1.5), 1.5);
    float spiral_phase = angle - disk_rotation_sign() * orbital_freq * t;

    // MAD spiral arm: dominant m = 1 mode from flux tube interaction
    float mad_spiral = 0.0;
    if (grmhd_mad_flux > 0.1) {
        float sp1 = sin(spiral_phase + 2.5 * log(max(r, 1.0)));
        float sp2 = sin(2.0 * spiral_phase + 4.0 * log(max(r, 1.0)) + 1.3);
        mad_spiral = grmhd_mad_flux * (0.3 * max(sp1, 0.0) + 0.15 * max(sp2, 0.0));
    }

    // Log-normal fluctuations calibrated to GRMHD σ_ln(ρ) ≈ 0.5–1.0
    // accretion_turbulence averages ~0.85, center log-normal there
    // so the mean output ≈ 1.0 (no systematic brightness bias)
    float sigma_ln = 0.4 + 0.6 * grmhd_turbulence_amp;
    float log_fluct = exp(sigma_ln * (base - 0.85) * 2.0);
    log_fluct = clamp(log_fluct, 0.2, 5.0);

    return log_fluct * (1.0 + mad_spiral);
}

// --- Non-thermal electron boost (κ-distribution) ---
// Magnetic reconnection in GRMHD produces non-thermal power-law
// tails on the electron energy distribution (Ball+ 2018, Werner+ 2018).
// κ → ∞ is thermal; κ ~ 3.5 has significant non-thermal synchrotron.
float grmhd_nonthermal_boost(float kappa) {
    float k = max(kappa, 2.5);
    return 1.0 + 2.0 / (k - 2.0);
}

// --- Fishbone-Moncrief equilibrium torus profile ---
// Standard initial condition for GRMHD simulations (FM76).
// The enthalpy distribution creates a peaked density at the
// pressure maximum radius r0.
float grmhd_fishbone_moncrief(float cyl_r, float z, float r0) {
    float r3d = sqrt(cyl_r * cyl_r + z * z);
    float sin_theta = cyl_r / max(r3d, 0.01);
    // Pseudo-Newtonian potential (Paczyński-Wiita approximation)
    float W = -1.0 / (2.0 * max(r3d - 1.0, 0.01))
            + 0.5 / (cyl_r * cyl_r + 0.01)
            * pow(cyl_r * sin_theta, 2.0);
    float W0 = -1.0 / (2.0 * max(r0 - 1.0, 0.01))
             + 0.5 / (r0 * r0 + 0.01) * r0 * r0;
    float dW = W - W0;
    // Only positive enthalpy region is inside the torus
    float rho = max(-dW * 3.0, 0.0);
    return pow(rho, 1.5) * grmhd_density_scale;
}

// --- Magnetic stress at ISCO (non-zero-torque boundary) ---
// GRMHD simulations show that magnetic stress at the ISCO is NOT zero
// (contradicting the Novikov-Thorne "no-torque" condition).
// This adds ~10-30% extra luminosity from the plunging region.
// (Noble et al. 2010, Penna et al. 2010)
float grmhd_isco_stress_factor(float r) {
    if (r >= ACCRETION_MIN_R) return 1.0;
    // Inside ISCO: decaying emission (not zero!)
    float f = r / max(ACCRETION_MIN_R, 1.0);
    // Magnetic stress persists into plunging region
    float stress = pow(f, 1.5) * (1.0 + 0.5 * grmhd_mad_flux);
    return stress;
}

// --- GRMHD 3D volumetric density turbulence ---
// Multi-scale FBM in co-rotating 3D coordinates, producing the filamentary,
// clumpy density structure characteristic of MRI-turbulent accretion flows.
// Calibrated to log-normal density PDFs from numerical GRMHD simulations:
//   σ_ln(ρ) ≈ 0.5–1.0 in the turbulent disk body (Hawley & Balbus 2002,
//   Beckwith et al. 2008, Hogg & Reynolds 2018).
// Spiral arm structure from dominant m=1,2 MRI modes; in MAD state,
// an additional m=1 spiral from the arrested magnetic flux tube
// (Tchekhovskoy, Narayan & McKinney 2011; Ripperda+ 2022).
float grmhd_3d_density_turbulence(vec3 p, float t) {
    float cyl_r = length(p.xy);
    float angle = equatorial_azimuth(p.xy);

    // Co-moving orbital advection (Keplerian shear)
    float orbit_phase = angle - disk_rotation_sign() * 0.45 * t / pow(max(cyl_r, 1.001), 1.5);
    vec2 orbit_unit = vec2(cos(orbit_phase), sin(orbit_phase));
    float z_sc = p.z * 1.5;

    // MRI spiral arms: dominant m=1,2 azimuthal modes
    // (Balbus & Hawley 1998, Beckwith et al. 2011)
    float sp1 = sin(orbit_phase + 3.0 * log(max(cyl_r, 1.0)));
    float sp2 = sin(2.0 * orbit_phase + 5.5 * log(max(cyl_r, 1.0)) + 1.7);
    float spirals = 0.45 * smoothstep(-0.3, 1.0, sp1)
                  + 0.35 * smoothstep(-0.3, 1.0, sp2);

    // MAD state: dominant m=1 spiral from arrested magnetic flux tube
    if (grmhd_mad_flux > 0.1) {
        float sp_mad = sin(orbit_phase + 2.5 * log(max(cyl_r, 1.0)));
        spirals += grmhd_mad_flux * 0.7 * smoothstep(-0.1, 1.0, sp_mad);
    }

    // Multi-scale FBM in co-rotating frame (2 octaves → 8 noise evals)
    // Large scale: ~ few r_g blobs and density waves
    float n1 = fbm(vec2(
        cyl_r * 2.0 + orbit_unit.x * 4.0,
        orbit_unit.y * 4.0 + z_sc * 1.5 + t * 0.05
    ));
    // Medium scale: filaments, streams, magnetic flux bundles
    float n2 = fbm(vec2(
        cyl_r * 8.0 + orbit_unit.x * 12.0 + n1 * 2.0,
        orbit_unit.y * 12.0 + z_sc * 5.0 - t * 0.15
    ));
    // Small scale (chaos): rapid turbulent shredding, especially intense near ISCO
    // and highly dynamic (fast t modulation).
    float isco_factor = exp(-0.8 * max(cyl_r - ACCRETION_MIN_R, 0.0));
    float n3 = fbm(vec2(
        cyl_r * 20.0 + orbit_unit.x * 35.0 - n2 * 4.0,
        orbit_unit.y * 35.0 + z_sc * 12.0 + t * 0.55
    ));

    // Log-normal density fluctuations: high dynamic range for the optically
    // thick slim disk where Beer-Lambert absorption renders surface features.
    // (Only the slim disk uses this function; torus uses grmhd_mri_turbulence.)
    float sigma_ln = 0.8 + 1.5 * grmhd_turbulence_amp;
    float combined = 0.35 * n1 + 0.35 * n2 + (0.15 + 0.4 * isco_factor) * n3;
    float log_density = exp(sigma_ln * (combined - 0.5) * 4.5);

    return log_density * (0.65 + 1.2 * spirals);
}

// --- GRMHD disk height modulation ---
// MRI turbulence and Parker instability cause the scale height to vary
// with azimuthal angle, creating warps, buoyant loops, and vertically
// extended features (Zhu & Stone 2018, Liska et al. 2022).
// Returns a multiplicative factor for the local disk scale height.
// Variation is ±15% at default turbulence (calibrated to Liska+ 2022
// Fig. 5: azimuthal H/R variation ~10-20% around the mean).
float grmhd_height_modulation(float cyl_r, float angle, float t) {
    float orbit_phase = angle - disk_rotation_sign() * 0.45 * t / pow(max(cyl_r, 1.001), 1.5);
    vec2 orbit_unit = vec2(cos(orbit_phase), sin(orbit_phase));

    float warp = fbm(vec2(
        cyl_r * 0.9 + orbit_unit.x * 2.0,
        orbit_unit.y * 2.0 + t * 0.05
    ));

    // Centered at 1.0: height varies by ±(15% * turbulence_amp)
    // warp ∈ [0,1], so (warp - 0.5) ∈ [-0.5, 0.5]
    float variation = 0.5 * (warp - 0.5) * grmhd_turbulence_amp;
    return 1.0 + variation;
}
{{/grmhd_enabled}}
