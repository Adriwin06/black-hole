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

float accretion_turbulence(float radius, float angle, float t) {
    float orbit_phase = angle - 0.45*t / pow(max(radius, 1.001), 1.5);
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
// Radiation emitted from the inner disk is strongly lensed by the black hole's
// gravity, causing a significant fraction of it to strike the disk again.
// Based on Cunningham (1976), this returning radiation enhances the
// local flux and temperature, peaking near the ISCO and scaling with spin.
float accretion_returning_radiation_enhancement(float radius) {
{{#disk_self_irradiation_enabled}}
    float r_norm = max(radius / ACCRETION_MIN_R, 1.0001);
    // Closer ISCO (higher spin) means the potential well is deeper, bending more light back.
    // For a=0 (Schwarzschild), ~20% enhancement. For a=0.99 (Extreme Kerr), > 100% enhancement.
    float spin_a = bh_rotation_enabled * bh_spin;
    float peak_enhancement = 0.2 + 1.2 * abs(spin_a);
    // The returning radiation flux decays rapidly outward as ~ r^-3.5
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

// Observer metric factor shared by both shift functions.
// Uses the static Schwarzschild factor because the observer's own
// kinematic Doppler is applied separately via cam_vel aberration.
float _observer_metric_factor() {
    float r_obs = max(length(cam_pos), 0.05);
    return max(abs(1.0 - 1.0 / r_obs), 0.001);
}

float gravitational_shift(float emission_radius) {
    // Combined gravitational + orbital redshift for matter on a circular
    // Keplerian orbit in Schwarzschild geometry.  The emitter redshift is:
    //   g_emit = sqrt(1 - 3M/r) = sqrt(1 - 3/(2r))   (r_s = 1, M = 0.5)
    // This accounts for BOTH the gravitational time dilation AND the
    // transverse Doppler effect of the orbiting disk material.
    // Below the photon sphere (r < 1.5 r_s) the quantity goes negative
    // (no stable circular orbits) — clamp with a small positive floor.
    float r_emit = max(emission_radius, 1.0001);
    float emission_term = max(1.0 - 1.5 / r_emit, 0.0001);
    return sqrt(emission_term / _observer_metric_factor());
}

float gravitational_shift_static(float emission_radius) {
    // Gravitational redshift only — for emitters that are NOT on circular
    // Keplerian orbits (ADAF torus, slim-disk plunging region, jets).
    // Uses the static Schwarzschild metric component sqrt(|1 - r_s/r|).
    // The emitter's kinematic Doppler is computed separately in trace_ray.
    float r_emit = max(emission_radius, 1.0001);
    float emission_term = max(abs(1.0 - 1.0 / r_emit), 0.0001);
    return sqrt(emission_term / _observer_metric_factor());
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

{{#accretion_slim_disk}}
// Slim disk: super-Eddington accretion, extends inside ISCO
// Radiation-pressure supported, geometrically thicker than thin disk
float slim_disk_height(float cyl_r) {
    // User-configurable base H/R ratio (default ~0.15, range 0.05-0.5)
    float base_h = max(slim_h_ratio, 0.05);
    // Puffs up near and inside ISCO due to radiation pressure
    float puff = max(slim_puff_factor, 0.0);
    float isco_proximity = exp(-1.5 * max(cyl_r - ACCRETION_MIN_R, 0.0));
    return max(cyl_r * base_h * (1.0 + puff * isco_proximity), 0.01);
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
// The electron-to-ion temperature ratio depends on the local plasma β:
//   R(β) = R_high × β² / (1 + β²) + R_low / (1 + β²)
// where R_low = 1 (electrons ≈ ions at equipartition when β → 0).
// This naturally creates:
//   • Hot electrons in the jet funnel (low β → R ≈ 1 → T_e ≈ T_i/2)
//   • Cool electrons in the disk body (high β → R ≈ R_high → T_e << T_i)
// When R_high is high, emission becomes concentrated near the jet wall
// and polar regions, producing the crescent morphology seen in EHT images.
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
    float orbit_phase = angle - 0.3 * t / pow(max(cyl_r, 1.0), 1.5);
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
    // Te_ratio ranges from ~0.5 (β→0) to ~1/(1+R_high) (β→∞)
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
    float spiral_phase = angle - orbital_freq * t;

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
    float angle = atan(p.x, p.y);

    // Co-moving orbital advection (Keplerian shear)
    float orbit_phase = angle - 0.45 * t / pow(max(cyl_r, 1.001), 1.5);
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
    float orbit_phase = angle - 0.45 * t / pow(max(cyl_r, 1.001), 1.5);
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
    return torus_temperature(r_t) * gravitational_shift_static(r_t);
    {{/accretion_thick_torus}}

    {{#accretion_slim_disk}}
    float r_s = ACCRETION_MIN_R * 1.6;
    return slim_disk_temperature(r_s) * gravitational_shift_static(r_s);
    {{/accretion_slim_disk}}

    return disk_temperature;
}
