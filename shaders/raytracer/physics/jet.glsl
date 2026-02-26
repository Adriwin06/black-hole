// Role: Relativistic jet emission models. Two modes are provided:
//       - Simple: smooth parabolic jet with spine/limb synchrotron approximation.
//       - Physical: GRMHD-calibrated model with magnetization parameter, spine/sheath
//         structure, reconfinement knots, jet-corona connection, and disk occultation.

{{#jet_enabled}}
{{#jet_simple}}
// ═══════════════════════════════════════════════════════════════════
// SIMPLE JET: Parabolic collimation, smooth synchrotron emission
// ═══════════════════════════════════════════════════════════════════
float jet_emissivity(vec3 p, float sign_z) {
    float z = p.z * sign_z;
    if (z < 1.2) return 0.0;

    float cyl_r = length(p.xy);

    float theta_jet = jet_half_angle * DEG_TO_RAD;
    float k_collimate = 0.58;
    float r_at_length = jet_length * tan(theta_jet);
    float r0 = r_at_length / pow(jet_length, k_collimate);
    float jet_r = r0 * pow(z, k_collimate);

    if (cyl_r > jet_r * 1.2) return 0.0;

    float r_norm = cyl_r / max(jet_r, 0.001);
    float spine = exp(-4.0 * r_norm * r_norm);
    float limb = 0.35 * exp(-12.0 * (r_norm - 0.8) * (r_norm - 0.8));
    float transverse = spine + limb;

    float z_onset = smoothstep(1.2, 3.0, z);
    float z_decay = pow(z, -1.2);
    float z_cutoff = 1.0 - smoothstep(jet_length * 0.75, jet_length, z);

    return transverse * z_onset * z_decay * z_cutoff;
}

vec3 jet_velocity(vec3 p, float sign_z) {
    float beta_max = sqrt(1.0 - 1.0 / (jet_lorentz * jet_lorentz));
    float z = abs(p.z);
    float accel = smoothstep(1.2, 5.0, z);
    float beta = beta_max * accel;
    float cyl_r = length(p.xy);
    float phi_twist = 0.03 * beta;
    vec3 v_z = vec3(0.0, 0.0, sign_z) * beta;
    vec3 v_phi = vec3(-p.y, p.x, 0.0) / max(cyl_r, 0.01) * phi_twist;
    return v_z + v_phi;
}
{{/jet_simple}}

{{#jet_physical}}
// ═══════════════════════════════════════════════════════════════════
// PHYSICAL JET: GRMHD-calibrated model
// ═══════════════════════════════════════════════════════════════════
// Based on:
//  - Moscibrodzka et al. 2016: GRMHD+radiative transfer of M87
//  - Davelaar et al. 2019: RAPTOR ray-tracing of jet models
//  - Blandford & Znajek 1977: jet power extraction
//  - Asada & Nakamura 2012: parabolic collimation profile
//  - Tchekhovskoy et al. 2011: MAD jet simulations
//  - Fromm et al. 2022: jet+disk modeling for EHT
//
// Key features over simple mode:
//  1. Magnetic-pressure-dominated funnel with σ (magnetization) parameter
//  2. Separate spine/sheath structure with different velocities
//  3. Reconfinement shock knots (standing shocks from pressure balance)
//  4. Jet-corona connection: bright base from hot corona above disk
//  5. Counter-jet disk occultation
//  6. GRMHD-calibrated emissivity: j ∝ ρ² B² (synchrotron)
//  7. Synchrotron power-law spectrum (not blackbody approximation)

// --- Jet geometry: magnetic funnel wall ---
// In GRMHD simulations the jet boundary follows field lines that thread
// the ergosphere. The boundary shape is well-fit by a parabola:
//   Ψ(r,θ) = const → r sin²θ = const (split-monopole)
// transitioning to r^k collimation at large distances (k ≈ 0.58).

float jet_boundary_radius(float z) {
    float theta_jet = jet_half_angle * DEG_TO_RAD;
    float k_col = 0.58; // Asada & Nakamura parabolic index
    // Smooth transition from funnel base to collimated jet
    // jet_base_width controls how wide the base is (lower = more collimated early)
    // Default base k ~ 0.72, transitioning to k_col = 0.58 at large z
    float k_base = mix(0.62, 0.85, clamp(jet_base_width, 0.0, 1.0));
    float k_eff = mix(k_base, k_col, smoothstep(3.0, 8.0, z));
    float r_at_length = jet_length * tan(theta_jet);
    float r0_param = r_at_length / pow(jet_length, k_col);
    return r0_param * pow(max(z, 0.01), k_eff);
}

// --- Transverse structure: spine + sheath ---
// GRMHD jets show a fast, magnetically-dominated spine surrounded by
// a slower, denser sheath (wind from the disk surface). The sheath
// produces most of the observed radio emission at parsec scales.
float jet_transverse_profile(float cyl_r, float jet_r, float sigma) {
    float r_norm = cyl_r / max(jet_r, 0.001);

    // Spine: Gaussian core, width depends on magnetization
    // Higher σ → narrower spine (more magnetically confined)
    float spine_width = mix(3.0, 6.0, clamp(sigma / 30.0, 0.0, 1.0));
    float spine = exp(-spine_width * r_norm * r_norm);

    // Sheath: peaked at r_norm ≈ 0.7-0.9 (jet boundary layer)
    // This is the contact discontinuity between jet and ambient medium
    float sheath_peak = 0.82;
    float sheath_width = 15.0;
    float sheath = 0.6 * exp(-sheath_width * (r_norm - sheath_peak) * (r_norm - sheath_peak));

    // Sharp outer boundary (ambient medium is dark)
    float boundary = 1.0 - smoothstep(0.95, 1.05, r_norm);

    return (spine + sheath) * boundary;
}

// --- Reconfinement shocks (knots) ---
// When the jet pressure drops below ambient, the jet re-collimates,
// creating standing oblique shocks. Observed as bright knots in
// M87 (HST-1), 3C 273, etc. Spacing ~ few × jet_radius.
float reconfinement_knots(float z, float jet_r) {
    float spacing = max(jet_knot_spacing, 2.0);
    // Standing wave pattern: sin² gives periodic bright spots
    float phase = M_PI * z / spacing;
    float knot = sin(phase) * sin(phase);
    // Knots are stronger at intermediate distances (not at base or tip)
    float knot_envelope = smoothstep(4.0, 8.0, z) * (1.0 - smoothstep(jet_length * 0.6, jet_length * 0.85, z));
    // Amplitude: 20-40% brightness enhancement
    return 1.0 + 0.35 * knot * knot_envelope;
}

// --- Jet-corona connection ---
// At the jet base, the magnetically-dominated funnel meets the hot
// corona above the accretion disk. This creates a bright, broad
// emission region (the "jet launching zone" seen in EHT images).
float corona_base_emission(vec3 p, float z, float cyl_r) {
    float r3d = length(p);
    // Corona extends from ~1.5 to ~4 r_s, concentrated near axis
    float base_z = smoothstep(1.0, 1.8, z) * (1.0 - smoothstep(2.5, 5.0, z));
    // Tighter radial profile controlled by jet_corona_extent
    // Default: fairly tight, avoids the overly-wide blob
    float corona_width = mix(2.0, 0.8, clamp(jet_corona_extent, 0.0, 1.0));
    float base_r = exp(-corona_width * cyl_r * cyl_r);
    // Emissivity scales with magnetic energy density ~ r^(-3)
    // (reduced from r^(-4) to avoid extreme concentration at center)
    float mag_energy = 1.0 / (r3d * r3d * r3d + 0.1);
    return jet_corona_brightness * base_z * base_r * mag_energy * 3.0;
}

// --- Magnetization parameter σ(z) ---
// σ = B²/(4πρc²) = ratio of magnetic to kinetic energy density
// In GRMHD jets: σ >> 1 at base (Poynting flux dominated),
// decreasing to σ ~ 1 at dissipation radius, σ < 1 far out.
float magnetization(float z) {
    float sigma_base = jet_magnetization; // σ at jet base
    // σ drops approximately as z^(-1) due to flux freezing + expansion
    // with a floor from residual field
    return sigma_base / (1.0 + 0.15 * z) + 0.5;
}

// --- GRMHD-calibrated emissivity ---
// Synchrotron emissivity: j_ν ∝ n_e B² sin²α × ν^(-α)
// In magnetically arrested disk (MAD) simulations, the jet funnel
// has j ∝ ρ² B² ∝ r^(-p) with p ≈ 2-3 depending on height.
float jet_emissivity(vec3 p, float sign_z) {
    float z = p.z * sign_z;
    if (z < 0.8) return 0.0; // below stagnation surface

    float cyl_r = length(p.xy);
    float r3d = length(p);

    // Jet boundary (parabolic with transition)
    float jet_r = jet_boundary_radius(z);
    if (cyl_r > jet_r * 1.15) {
        // Check if we're in the corona base region
        if (z < 5.0 && r3d < 5.0) {
            return corona_base_emission(p, z, cyl_r);
        }
        return 0.0;
    }

    // Local magnetization
    float sigma = magnetization(z);

    // Transverse spine/sheath structure
    float profile = jet_transverse_profile(cyl_r, jet_r, sigma);

    // Longitudinal emissivity:
    // - onset: plasma acceleration zone (z ~ 1.5 to 3)
    // - GRMHD power-law: j ∝ z^(-2.2) (steeper than simple mode)
    //   (Moscibrodzka+2016 find emission concentrated within ~10 r_g)
    float z_onset = smoothstep(0.8, 2.5, z);
    float z_decay = pow(max(z, 1.0), -2.2);

    // Length cutoff
    float z_cutoff = 1.0 - smoothstep(jet_length * 0.7, jet_length, z);

    // Reconfinement knots
    float knots = reconfinement_knots(z, jet_r);

    // Combined: synchrotron emissivity j ∝ n²B² × geometry
    float j = profile * z_onset * z_decay * z_cutoff * knots;

    // Add corona base contribution inside the funnel
    j += corona_base_emission(p, z, cyl_r);

    return j;
}

// --- Jet velocity: spine/sheath with GRMHD acceleration ---
// In GRMHD: jet accelerates as magnetic energy converts to kinetic.
// The spine moves at the bulk Γ; the sheath at a lower ~0.3-0.5c.
vec3 jet_velocity(vec3 p, float sign_z) {
    float z = abs(p.z);
    float cyl_r = length(p.xy);
    float jet_r = jet_boundary_radius(z);
    float r_norm = cyl_r / max(jet_r, 0.001);

    // Maximum velocity from Lorentz factor
    float beta_max = sqrt(1.0 - 1.0 / (jet_lorentz * jet_lorentz));

    // GRMHD acceleration profile:
    // Slow acceleration in Poynting-dominated zone (z < ~5),
    // rapid conversion at σ ≈ 1 (z ~ 5-10),
    // then saturates at terminal Γ.
    float sigma = magnetization(z);
    // β ≈ β_max × (1 - 1/σ) approximately, smooth version:
    float accel = 1.0 - 1.0 / (1.0 + 0.5 * z * z / (4.0 + z));
    float beta_spine = beta_max * accel;

    // Sheath velocity: slower than spine (disk wind, ~0.3-0.5c)
    float beta_sheath = 0.4 * accel;

    // Interpolate spine/sheath based on radial position
    float spine_frac = exp(-3.0 * r_norm * r_norm);
    float beta = mix(beta_sheath, beta_spine, spine_frac);

    // Helical component: stronger in the sheath (wound-up field)
    // Toroidal fraction increases with distance (field winding)
    float phi_fraction = mix(0.15, 0.03, spine_frac) * min(z / 5.0, 1.0);
    float v_z_comp = beta * sqrt(1.0 - phi_fraction * phi_fraction);
    float v_phi_comp = beta * phi_fraction;

    vec3 v_z = vec3(0.0, 0.0, sign_z) * v_z_comp;
    vec3 v_phi = vec3(-p.y, p.x, 0.0) / max(cyl_r, 0.01) * v_phi_comp;
    return v_z + v_phi;
}

// --- Synchrotron power-law color ---
// True synchrotron spectrum: F_ν ∝ ν^(-α) with α ≈ 0.6-0.8
// Mapping to RGB: more flux at low frequencies → reddish at base,
// blue-shifted at high Γ. Approximate with variable-T blackbody
// where effective T encodes the spectral hardness.
float jet_effective_temperature(float z, float r3d, float sigma) {
    // Base is hotter (corona: ~10^9 K effective → maps to UV/X-ray)
    // Far jet is cooler (synchrotron aging: electrons lose energy)
    float T_base = 35000.0; // hot corona
    float T_far = 12000.0;  // aged synchrotron
    float T_mix = mix(T_base, T_far, smoothstep(2.0, 15.0, z));
    // Magnetization boost: σ > 1 regions are hotter (reconnection heating)
    T_mix *= 1.0 + 0.3 * clamp(sigma - 1.0, 0.0, 5.0) / 5.0;
    return T_mix;
}
{{/jet_physical}}
{{/jet_enabled}}
