#define M_PI 3.141592653589793238462643383279
#define R_SQRT_2 0.7071067811865475
#define DEG_TO_RAD (M_PI/180.0)
#define SQ(x) ((x)*(x))

#define ROT_Y(a) mat3(0, cos(a), sin(a), 1, 0, 0, 0, sin(a), -cos(a))


// spectrum texture lookup helper macros
const float BLACK_BODY_TEXTURE_COORD = 1.0;
const float SINGLE_WAVELENGTH_TEXTURE_COORD = 0.5;
const float TEMPERATURE_LOOKUP_RATIO_TEXTURE_COORD = 0.0;

// black-body texture metadata
const float SPECTRUM_TEX_TEMPERATURE_RANGE = 65504.0;
const float SPECTRUM_TEX_WAVELENGTH_RANGE = 2048.0;
const float SPECTRUM_TEX_RATIO_RANGE = 6.48053329012;

// multi-line macros don't seem to work in WebGL :(
#define BLACK_BODY_COLOR(t) texture2D(spectrum_texture, vec2((t) / SPECTRUM_TEX_TEMPERATURE_RANGE, BLACK_BODY_TEXTURE_COORD))
#define SINGLE_WAVELENGTH_COLOR(lambda) texture2D(spectrum_texture, vec2((lambda) / SPECTRUM_TEX_WAVELENGTH_RANGE, SINGLE_WAVELENGTH_TEXTURE_COORD))
#define TEMPERATURE_LOOKUP(ratio) (texture2D(spectrum_texture, vec2((ratio) / SPECTRUM_TEX_RATIO_RANGE, TEMPERATURE_LOOKUP_RATIO_TEXTURE_COORD)).r * SPECTRUM_TEX_TEMPERATURE_RANGE)

uniform vec2 resolution;
uniform float time;

uniform vec3 cam_pos;
uniform vec3 cam_x;
uniform vec3 cam_y;
uniform vec3 cam_z;
uniform vec3 cam_vel;

uniform float planet_distance, planet_radius;
uniform float disk_temperature;
uniform float bh_spin, bh_spin_strength, bh_rotation_enabled;
uniform float accretion_inner_r;
uniform float look_exposure, look_disk_gain, look_glow, look_doppler_boost;
uniform float look_aberration_strength;
uniform float look_star_gain, look_galaxy_gain;
uniform float look_tonemap_mode;
uniform float torus_r0, torus_h_ratio;
uniform float jet_half_angle, jet_lorentz, jet_brightness, jet_length;
{{#jet_physical}}
uniform float jet_magnetization, jet_knot_spacing, jet_corona_brightness;
{{/jet_physical}}

uniform sampler2D galaxy_texture, star_texture,
    planet_texture, spectrum_texture;

// stepping and anti-aliasing parameters
const int NSTEPS = {{n_steps}};
const int SAMPLE_COUNT = {{sample_count}};
const float MAX_REVOLUTIONS = float({{max_revolutions}});

// ACCRETION_MIN_R is now a uniform (accretion_inner_r) that varies with black hole spin
// Using ISCO from Bardeen-Press-Teukolsky formula: r_ISCO = 3 r_s for Schwarzschild
const float ACCRETION_WIDTH = 12.0;
#define ACCRETION_MIN_R accretion_inner_r
const float ACCRETION_BRIGHTNESS = 0.95;

const float STAR_MIN_TEMPERATURE = 4000.0;
const float STAR_MAX_TEMPERATURE = 15000.0;

const float STAR_BRIGHTNESS = 0.52;
const float GALAXY_BRIGHTNESS = 0.14;
const float GLOBAL_EXPOSURE = 0.60;
const float GALAXY_DOPPLER_STRENGTH = 1.0;
const float GALAXY_MAX_BOOST = 10.0;

const float PLANET_AMBIENT = 0.0;
const float PLANET_LIGHTNESS = 1.5;
const float PLANET_LIGHT_REF_DISTANCE = 14.0;

// background texture coordinate system
mat3 BG_COORDS = ROT_Y(45.0 * DEG_TO_RAD);

// planet texture coordinate system
const float PLANET_AXIAL_TILT = 30.0 * DEG_TO_RAD;
mat3 PLANET_COORDS = ROT_Y(PLANET_AXIAL_TILT);

const float FOV_ANGLE_DEG = 90.0;
float FOV_MULT = 1.0 / tan(DEG_TO_RAD * FOV_ANGLE_DEG*0.5);

// derived "constants" (from uniforms)
float PLANET_RADIUS,
    PLANET_DISTANCE,
    PLANET_ORBITAL_ANG_VEL,
    PLANET_ROTATION_ANG_VEL,
    PLANET_GAMMA;

vec2 sphere_map(vec3 p) {
    return vec2(atan(p.x,p.y)/M_PI*0.5+0.5, asin(p.z)/M_PI+0.5);
}

float smooth_step(float x, float threshold) {
    const float STEEPNESS = 1.0;
    return 1.0 / (1.0 + exp(-(x-threshold)*STEEPNESS));
}

vec3 lorentz_velocity_transformation(vec3 moving_v, vec3 frame_v) {
    float v = length(frame_v);
    if (v > 0.0) {
        vec3 v_axis = -frame_v / v;
        float gamma = 1.0/sqrt(1.0 - v*v);

        float moving_par = dot(moving_v, v_axis);
        vec3 moving_perp = moving_v - v_axis*moving_par;

        float denom = 1.0 + v*moving_par;
        return (v_axis*(moving_par+v)+moving_perp/gamma)/denom;
    }
    return moving_v;
}

vec3 contract(vec3 x, vec3 d, float mult) {
    float par = dot(x,d);
    return (x-par*d) + d*par*mult;
}

vec3 safe_normalize(vec3 v) {
    float l = length(v);
    if (l > 1e-6) return v/l;
    return vec3(0.0, 0.0, 0.0);
}

vec3 rotate_about_z(vec3 p, float angle) {
    float c = cos(angle);
    float s = sin(angle);
    return vec3(
        c*p.x - s*p.y,
        s*p.x + c*p.y,
        p.z
    );
}

const float KERR_M = 0.5; // r_s = 1 => M = r_s/2

float kerr_spin_a() {
    return bh_rotation_enabled * bh_spin * KERR_M;
}

float kerr_delta(float r, float a) {
    return r*r - r + a*a; // r^2 - 2Mr + a^2, with 2M = 1
}

float kerr_sigma(float r, float theta, float a) {
    float c = cos(theta);
    return r*r + a*a*c*c;
}

float kerr_horizon_radius(float a) {
    return KERR_M + sqrt(max(KERR_M*KERR_M - a*a, 0.0));
}

float geodesic_accel(float u, float spin_alignment) {
    // Schwarzschild photon geodesic (Binet equation): d²u/dφ² = -u + (3/2)u² in units where r_s = 1
    // Photon sphere at u = 2/3 (r = 1.5 r_s), see e.g. MTW "Gravitation" or
    // https://en.wikipedia.org/wiki/Schwarzschild_geodesics#Bending_of_light_by_gravity
    float schwarzschild_accel = -u + 1.5*u*u;
    
    // Improved Kerr frame-dragging approximation
    // In true Kerr, frame dragging affects photon trajectories via the metric term g_tφ
    // The effect scales as a/r³ where a is the spin parameter
    // This approximation captures the qualitative behavior: 
    // - Prograde light bends less (aligned with rotation)
    // - Retrograde light bends more (against rotation)
    // Using u³ term (1/r³ dependence) which is more physical than u⁴
    float frame_drag_term = bh_rotation_enabled * bh_spin * bh_spin_strength *
        spin_alignment * 0.8 * u*u*u;
    
    return schwarzschild_accel + frame_drag_term;
}

void integrate_geodesic_step(inout float u, inout float du, float step,
        float spin_alignment) {
    {{#rk4_integration}}
    float k1_u = du;
    float k1_du = geodesic_accel(u, spin_alignment);

    float u2 = u + 0.5*step*k1_u;
    float du2 = du + 0.5*step*k1_du;
    float k2_u = du2;
    float k2_du = geodesic_accel(u2, spin_alignment);

    float u3 = u + 0.5*step*k2_u;
    float du3 = du + 0.5*step*k2_du;
    float k3_u = du3;
    float k3_du = geodesic_accel(u3, spin_alignment);

    float u4 = u + step*k3_u;
    float du4 = du + step*k3_du;
    float k4_u = du4;
    float k4_du = geodesic_accel(u4, spin_alignment);

    u += (step/6.0) * (k1_u + 2.0*k2_u + 2.0*k3_u + k4_u);
    du += (step/6.0) * (k1_du + 2.0*k2_du + 2.0*k3_du + k4_du);
    {{/rk4_integration}}
    {{^rk4_integration}}
    u += du*step;
    du += geodesic_accel(u, spin_alignment)*step;
    {{/rk4_integration}}
}

float hash12(vec2 p) {
    vec3 p3 = fract(vec3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

float value_noise(vec2 p) {
    vec2 i = floor(p);
    vec2 f = fract(p);

    float a = hash12(i);
    float b = hash12(i + vec2(1.0, 0.0));
    float c = hash12(i + vec2(0.0, 1.0));
    float d = hash12(i + vec2(1.0, 1.0));

    vec2 u = f*f*(3.0 - 2.0*f);
    return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

float fbm(vec2 p) {
    float sum = 0.0;
    float amp = 0.5;
    for (int i = 0; i < 4; i++) {
        sum += amp * value_noise(p);
        p = p*2.03 + vec2(17.13, -11.70);
        amp *= 0.5;
    }
    return sum;
}

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
    // Smooth transition from wide base to collimated jet
    // At z < 3: wider funnel (split-monopole-like, k~1)
    // At z > 5: parabolic collimation (k=0.58)
    float k_eff = mix(0.9, k_col, smoothstep(2.0, 6.0, z));
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
    // Corona extends from ~1.5 to ~4 r_s, broad and bright
    float base_z = smoothstep(1.0, 1.8, z) * (1.0 - smoothstep(2.5, 5.0, z));
    // Wider than the jet itself at the base
    float base_r = exp(-0.6 * cyl_r * cyl_r);
    // Emissivity scales with magnetic energy density ~ r^(-4)
    float mag_energy = 1.0 / (r3d * r3d * r3d * r3d + 0.01);
    return jet_corona_brightness * base_z * base_r * mag_energy * 8.0;
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

vec2 sample_offset(int i, vec2 pixel) {
    if (SAMPLE_COUNT <= 1) return vec2(0.0, 0.0);
    float fi = float(i);
    float sample_count_f = float(SAMPLE_COUNT);
    float radius = 0.5 * sqrt((fi + 0.5) / sample_count_f);
    float base_angle = 2.0*M_PI*fract(fi*0.61803398875 + hash12(pixel*0.5));
    return radius * vec2(cos(base_angle), sin(base_angle));
}

// --- ACES Filmic tonemapper ---
vec3 aces_filmic(vec3 x) {
    return clamp((x*(2.51*x + 0.03)) / (x*(2.43*x + 0.59) + 0.14), 0.0, 1.0);
}

// --- AGX tonemapper (Troy Sobotka) ---
// Better handling of high-saturation, high-luminance colors;
// avoids the hue shifts and excessive desaturation of ACES in
// extreme HDR scenes like accretion disks.
vec3 agx_default_contrast_approx(vec3 x) {
    vec3 x2 = x * x;
    vec3 x4 = x2 * x2;
    return + 15.5     * x4 * x2
           - 40.14    * x4 * x
           + 31.96    * x4
           - 6.868    * x2 * x
           + 0.4298   * x2
           + 0.1191   * x
           - 0.00232;
}

vec3 agx_tonemap(vec3 val) {
    // AGX input transform (column-major)
    const mat3 agx_mat = mat3(
        0.842479062253094,  0.0784335999999992, 0.0792237451477643,
        0.0423282422610123, 0.878468636469772,  0.0791661274605434,
        0.0423756549057051, 0.0784336,          0.879142973793104
    );

    // Standard AGX EV range (16.5 stops).
    const float min_ev = -12.47393;
    const float max_ev = 4.026069;

    val = agx_mat * max(val, vec3(1e-10));
    val = clamp(log2(val), min_ev, max_ev);
    val = (val - min_ev) / (max_ev - min_ev);
    val = agx_default_contrast_approx(val);
    return val;
}

vec3 agx_eotf(vec3 val) {
    // AGX inverse output transform (column-major)
    const mat3 agx_mat_inv = mat3(
        1.19687900512017,   -0.0980208811401368, -0.0990297440797205,
        -0.0528968517574562, 1.15190312990417,   -0.0989611768448433,
        -0.0529716355144438,-0.0980434501171241,  1.15107367264116
    );
    val = agx_mat_inv * val;
    return pow(max(val, vec3(0.0)), vec3(2.2));
}

vec3 agx_look_punchy(vec3 val) {
    // Photographic contrast around middle-gray pivot (0.18).
    // Raises highlights and deepens shadows while keeping mid-gray
    // anchored — this amplifies the bright/dim beaming ratio without
    // blowing out the highlights or inverting anything.
    const float contrast = 1.25;
    const float pivot    = 0.18;
    val = pow(max(val / pivot, 1e-6), vec3(contrast)) * pivot;
    val = clamp(val, 0.0, 1.0);

    // Saturation boost to make Doppler hue shifts more visible
    float luma = dot(val, vec3(0.2126, 0.7152, 0.0722));
    val = clamp(mix(vec3(luma), val, 1.35), 0.0, 1.0);
    return val;
}

// --- Inferno colormap (matplotlib) polynomial approximation ---
// Used for scientific false-color rendering
vec3 inferno(float t) {
    t = clamp(t, 0.0, 1.0);
    const vec3 c0 = vec3(0.0002189403691192265, 0.001651004631001012, -0.01948089843709184);
    const vec3 c1 = vec3(0.1065134194856116, 0.5639564367884091, 3.932712388889277);
    const vec3 c2 = vec3(11.60249308247187, -3.972853965665698, -15.9423941062914);
    const vec3 c3 = vec3(-41.70399613139459, 17.43639888205313, 44.35414519872813);
    const vec3 c4 = vec3(77.16275788913913, -33.40235894210092, -81.80730925738993);
    const vec3 c5 = vec3(-71.31942824499214, 32.62606426397723, 73.20951985803202);
    const vec3 c6 = vec3(25.13112622477341, -12.24266895238567, -23.07032500287172);
    return c0+t*(c1+t*(c2+t*(c3+t*(c4+t*(c5+t*c6)))));
}

// --- Logarithmic scientific tonemapper ---
// Mimics how scientific papers (EHT, Luminet, GRMHD simulations)
// display accretion disk data: log-scale intensity mapped through
// the inferno false-color palette.
//
// High sensitivity: typical post-exposure disk luminances (~0.01–0.5)
// are spread across the full purple → red → orange → yellow → white
// range.  The exposure slider shifts the saturation point.
vec3 tonemap_log_scientific(vec3 color, float exposure) {
    float lum = dot(color, vec3(0.2126, 0.7152, 0.0722));
    // Very high log_k so even faint disk emission is visible.
    // With exposure=1:  log_k=80, denom anchor = 0.15 * 80 = 12.
    //   lum ≈ 0.005 → t ≈ 0.13  (near-black)
    //   lum ≈ 0.02  → t ≈ 0.35  (dark red)
    //   lum ≈ 0.05  → t ≈ 0.55  (orange)
    //   lum ≈ 0.10  → t ≈ 0.73  (yellow)
    //   lum ≈ 0.15  → t ≈ 0.85  (bright yellow)
    //   lum ≈ 0.25+ → t → 1.0   (white)
    float log_k = max(exposure * 80.0, 0.01);
    float mapped = log2(1.0 + lum * log_k) / log2(1.0 + log_k * 0.15);
    mapped = clamp(mapped, 0.0, 1.0);
    return inferno(mapped);
}

float screen_dither() {
    return (hash12(gl_FragCoord.xy) - 0.5) / 255.0;
}

vec4 finalize_color(vec4 color) {
    {{#cinematic_tonemap}}
    vec3 exposed = max(color.rgb * (GLOBAL_EXPOSURE * look_exposure), vec3(0.0));
    vec3 mapped;

    if (look_tonemap_mode < 0.5) {
        // Mode 0 — ACES Filmic
        mapped = aces_filmic(exposed);
        float lum = dot(mapped, vec3(0.2126, 0.7152, 0.0722));
        mapped += mapped * smoothstep(0.82, 1.0, lum) * 0.06;
        mapped = pow(mapped, vec3(1.0/2.2));
    } else if (look_tonemap_mode < 1.5) {
        // Mode 1 — AGX
        // 1.15× nudge: AGX maps middle-gray slightly darker than ACES
        mapped = agx_tonemap(exposed * 1.15);
        mapped = agx_look_punchy(mapped);
        mapped = agx_eotf(mapped);
        mapped = pow(mapped, vec3(1.0/2.2));
    } else {
        // Mode 2 — Logarithmic Scientific (false-color)
        mapped = tonemap_log_scientific(exposed, look_exposure);
    }

    mapped += vec3(screen_dither());
    mapped = clamp(mapped, 0.0, 1.0);
    return vec4(mapped, 1.0);
    {{/cinematic_tonemap}}
    {{^cinematic_tonemap}}
    vec3 mapped = color.rgb + vec3(screen_dither());
    return vec4(clamp(mapped, 0.0, 1.0), color.a);
    {{/cinematic_tonemap}}
}

vec4 planet_intersection(vec3 old_pos, vec3 ray, float t, float dt,
        vec3 planet_pos0, float ray_doppler_factor) {

    vec4 ret = vec4(0,0,0,0);
    ray = ray/dt;

    vec3 planet_dir = vec3(planet_pos0.y, -planet_pos0.x, 0.0) / PLANET_DISTANCE;

    {{#light_travel_time}}
    float planet_ang1 = (t-dt) * PLANET_ORBITAL_ANG_VEL;
    vec3 planet_pos1 = vec3(cos(planet_ang1), sin(planet_ang1), 0)*PLANET_DISTANCE;
    vec3 planet_vel = (planet_pos1-planet_pos0)/dt;

    // transform to moving planet coordinate system
    ray = ray - planet_vel;
    {{/light_travel_time}}
    {{^light_travel_time}}
    vec3 planet_vel = planet_dir * PLANET_ORBITAL_ANG_VEL * PLANET_DISTANCE;
    {{/light_travel_time}}

    // ray-sphere intersection
    vec3 d = old_pos - planet_pos0;

    {{#lorentz_contraction}}
    ray = contract(ray, planet_dir, PLANET_GAMMA);
    d = contract(d, planet_dir, PLANET_GAMMA);
    {{/lorentz_contraction}}

    float dotp = dot(d,ray);
    float c_coeff = dot(d,d) - SQ(PLANET_RADIUS);
    float ray2 = dot(ray, ray);
    float discr = dotp*dotp - ray2*c_coeff;

    if (discr < 0.0) return ret;
    float isec_t = (-dotp - sqrt(discr)) / ray2;

    float MIN_ISEC_DT = 0.0;
    {{#lorentz_contraction}}
    MIN_ISEC_DT = -dt;
    {{/lorentz_contraction}}

    if (isec_t < MIN_ISEC_DT || isec_t > dt) return ret;

    vec3 surface_point = (d + isec_t*ray) / PLANET_RADIUS;
    isec_t = isec_t/dt;

    vec3 light_dir = planet_pos0;
    float rot_phase = t;

    {{#light_travel_time}}
    light_dir += planet_vel*isec_t*dt;
    rot_phase -= isec_t*dt;
    {{/light_travel_time}}

    rot_phase = rot_phase * PLANET_ROTATION_ANG_VEL*0.5/M_PI;
    light_dir = light_dir / PLANET_DISTANCE;

    {{#light_travel_time}}
    light_dir = light_dir - planet_vel;
    {{/light_travel_time}}

    vec3 surface_normal = surface_point;
    {{#lorentz_contraction}}
    light_dir = contract(light_dir, planet_dir, PLANET_GAMMA);
    {{/lorentz_contraction}}
    light_dir = normalize(light_dir);

    vec2 tex_coord = sphere_map(surface_point * PLANET_COORDS);
    tex_coord.x = mod(tex_coord.x + rot_phase, 1.0);

    float diffuse = max(0.0, dot(surface_normal, -light_dir));
    float distance_attenuation = SQ(
        PLANET_LIGHT_REF_DISTANCE / max(PLANET_DISTANCE, PLANET_RADIUS + 1.0)
    );
    distance_attenuation = clamp(distance_attenuation, 0.05, 4.0);

    float lightness = ((1.0-PLANET_AMBIENT)*diffuse + PLANET_AMBIENT) *
        PLANET_LIGHTNESS * distance_attenuation;

    float light_temperature = planet_irradiation_temperature();
    float transfer_factor = max(ray_doppler_factor, 0.05);
    {{#doppler_shift}}
    float doppler_factor = SQ(PLANET_GAMMA) *
        (1.0 + dot(planet_vel, light_dir)) *
        (1.0 - dot(planet_vel, normalize(ray)));
    transfer_factor = max(doppler_factor * ray_doppler_factor, 0.05);
    light_temperature /= transfer_factor;
    {{/doppler_shift}}
    {{#beaming}}
    {{#physical_beaming}}
    lightness /= pow(clamp(transfer_factor, 0.05, 20.0), 3.0);
    {{/physical_beaming}}
    {{^physical_beaming}}
    float clamped_planet_doppler = clamp(transfer_factor, 0.62, 1.48);
    lightness /= pow(clamped_planet_doppler, 1.05 + 1.10*doppler_boost);
    {{/physical_beaming}}
    {{/beaming}}

    vec3 blackbody_rgb = BLACK_BODY_COLOR(light_temperature).rgb;
    float bb_max = max(max(blackbody_rgb.r, blackbody_rgb.g), blackbody_rgb.b);
    vec3 physical_tint = vec3(1.0, 1.0, 1.0);
    if (bb_max > 1e-5) {
        physical_tint = blackbody_rgb / bb_max;
    }

    // Use physical blackbody tint for planet illumination
    vec3 light_tint = physical_tint;

    ret = texture2D(planet_texture, tex_coord) * lightness;
    ret.rgb *= light_tint;
    if (isec_t < 0.0) isec_t = 0.5;
    ret.w = isec_t;

    return ret;
}

vec4 galaxy_color(vec2 tex_coord, float doppler_factor) {

    vec4 base_color = texture2D(galaxy_texture, tex_coord);
    vec4 color = base_color;
    {{^observerMotion}}
    return color;
    {{/observerMotion}}

    {{#observerMotion}}
    vec4 ret = vec4(0.0,0.0,0.0,0.0);
    float red = max(0.0, color.r - color.g);

    const float H_ALPHA_RATIO = 0.1;
    const float TEMPERATURE_BIAS = 0.95;

    color.r -= red*H_ALPHA_RATIO;

    float i1 = max(color.r, max(color.g, color.b));
    float ratio = (color.g+color.b) / color.r;

    if (i1 > 0.0 && color.r > 0.0) {

        float temperature = TEMPERATURE_LOOKUP(ratio) * TEMPERATURE_BIAS;
        color = BLACK_BODY_COLOR(temperature);

        float i0 = max(color.r, max(color.g, color.b));
        if (i0 > 0.0) {
            temperature /= max(doppler_factor, 0.75);
            float remap_gain = clamp(i1 / max(i0, 0.18), 0.0, GALAXY_MAX_BOOST);
            ret = BLACK_BODY_COLOR(temperature) * remap_gain;
        }
    }

    ret += SINGLE_WAVELENGTH_COLOR(656.28 * doppler_factor) * red / 0.214 * H_ALPHA_RATIO;

    ret = mix(base_color, ret, GALAXY_DOPPLER_STRENGTH);
    ret.rgb = min(ret.rgb, base_color.rgb * GALAXY_MAX_BOOST + vec3(0.03));
    return ret;
    {{/observerMotion}}
}

vec4 trace_ray(vec3 ray) {
    vec3 pos = cam_pos;

    {{#aberration}}
    vec3 aberration_vel = cam_vel * max(look_aberration_strength, 0.0);
    float aberration_speed = length(aberration_vel);
    if (aberration_speed > 0.999) {
        aberration_vel *= 0.999 / aberration_speed;
    }
    ray = lorentz_velocity_transformation(ray, aberration_vel);
    {{/aberration}}

    float ray_intensity = 1.0;
    float ray_doppler_factor = 1.0;
    float doppler_boost = clamp(look_doppler_boost, 0.0, 2.5);

    float gamma = 1.0/sqrt(max(1.0-dot(cam_vel,cam_vel), 0.0001));
    ray_doppler_factor = gamma*(1.0 + dot(ray,-cam_vel));
    {{#beaming}}
    {{#physical_beaming}}
    // Observer beaming is already accounted for via transfer_factor
    // in each emission source. The temperature shift of blackbody sources
    // captures the full Doppler effect. No additional ray_intensity scaling.
    {{/physical_beaming}}
    {{^physical_beaming}}
    float beaming_factor = clamp(ray_doppler_factor, 0.84, 1.16);
    ray_intensity /= pow(beaming_factor, 0.65 + 0.75*doppler_boost);
    {{/physical_beaming}}
    {{/beaming}}
    {{^doppler_shift}}
    ray_doppler_factor = 1.0;
    {{/doppler_shift}}

    float step = 0.01;
    vec4 color = vec4(0.0,0.0,0.0,1.0);
    float vol_transmittance = 1.0; // volumetric ray transmittance (Beer-Lambert)

    // initial conditions
    float u = 1.0 / length(pos), old_u;
    float u0 = u;

    vec3 normal_vec = normalize(pos);
    vec3 tangent_vec = normalize(cross(cross(normal_vec, ray), normal_vec));
    vec3 spin_axis = vec3(0.0, 0.0, 1.0);
    float spin_alignment = clamp(dot(safe_normalize(cross(pos, ray)), spin_axis), -1.0, 1.0);
    float frame_drag_phase = 0.0;

    float du = -dot(ray,normal_vec) / dot(ray,tangent_vec) * u;
    float du0 = du;

    float phi = 0.0;

    float t = time;
    float dt = 1.0;
    bool shadow_capture = false;

    {{^light_travel_time}}
    float planet_ang0 = t * PLANET_ORBITAL_ANG_VEL;
    vec3 planet_pos0 = vec3(cos(planet_ang0), sin(planet_ang0), 0)*PLANET_DISTANCE;
    {{/light_travel_time}}

    vec3 old_pos = pos;

    for (int j=0; j < NSTEPS; j++) {
        {{#kerr_fast_mode}}
        step = MAX_REVOLUTIONS * 2.0*M_PI / float(NSTEPS);

        // adaptive step size based on rate of change
        float max_rel_u_change = (1.0-log(max(u, 0.0001)))*10.0 / float(NSTEPS);
        if ((du > 0.0 || (du0 < 0.0 && u0/u < 5.0)) && abs(du) > abs(max_rel_u_change*u) / step) {
            step = max_rel_u_change*u/abs(du);
        }

        // Additional step refinement near photon sphere (u ≈ 0.667 for r = 1.5)
        float u_photon_sphere = 0.667;
        float photon_sphere_proximity = exp(-12.0 * (u - u_photon_sphere) * (u - u_photon_sphere));
        step *= 1.0 - 0.7 * photon_sphere_proximity;

        old_u = u;

        {{#light_travel_time}}
        {{#gravitational_time_dilation}}
        dt = sqrt(du*du + u*u*(1.0-u))/(u*u*(1.0-u))*step;
        {{/gravitational_time_dilation}}
        {{/light_travel_time}}

        integrate_geodesic_step(u, du, step, spin_alignment);

        if (u < 0.0) break;
        if (u >= 1.0) {
            shadow_capture = true;
            break;
        }

        phi += step;

        old_pos = pos;
        vec3 planar_pos = (cos(phi)*normal_vec + sin(phi)*tangent_vec)/u;

        float drag_u = clamp(u, 0.0, 1.15);
        float frame_drag_step = bh_rotation_enabled * bh_spin * bh_spin_strength *
            spin_alignment * step * 0.85 * drag_u*drag_u*drag_u;
        frame_drag_phase += frame_drag_step;

        pos = rotate_about_z(planar_pos, frame_drag_phase);
        {{/kerr_fast_mode}}

        {{#kerr_full_core}}
        // Realtime Kerr: Binet equation with frame-drag approximation
        // Same geodesic integration as fast mode, with full Kerr disk kinematics
        step = MAX_REVOLUTIONS * 2.0*M_PI / float(NSTEPS);

        // adaptive step size based on rate of change
        float max_rel_u_change_fc = (1.0-log(max(u, 0.0001)))*10.0 / float(NSTEPS);
        if ((du > 0.0 || (du0 < 0.0 && u0/u < 5.0)) && abs(du) > abs(max_rel_u_change_fc*u) / step) {
            step = max_rel_u_change_fc*u/abs(du);
        }

        // Additional step refinement near photon sphere (u ≈ 0.667 for r = 1.5)
        float u_photon_sphere_fc = 0.667;
        float photon_sphere_proximity_fc = exp(-12.0 * (u - u_photon_sphere_fc) * (u - u_photon_sphere_fc));
        step *= 1.0 - 0.7 * photon_sphere_proximity_fc;

        old_u = u;

        {{#light_travel_time}}
        {{#gravitational_time_dilation}}
        dt = sqrt(du*du + u*u*(1.0-u))/(u*u*(1.0-u))*step;
        {{/gravitational_time_dilation}}
        {{/light_travel_time}}

        integrate_geodesic_step(u, du, step, spin_alignment);

        if (u < 0.0) break;
        if (u >= 1.0) {
            shadow_capture = true;
            break;
        }

        phi += step;

        old_pos = pos;
        vec3 planar_pos_fc = (cos(phi)*normal_vec + sin(phi)*tangent_vec)/u;

        float drag_u_fc = clamp(u, 0.0, 1.15);
        float frame_drag_step_fc = bh_rotation_enabled * bh_spin * bh_spin_strength *
            spin_alignment * step * 0.85 * drag_u_fc*drag_u_fc*drag_u_fc;
        frame_drag_phase += frame_drag_step_fc;

        pos = rotate_about_z(planar_pos_fc, frame_drag_phase);

        {{#light_travel_time}}
        {{^gravitational_time_dilation}}
        dt = length(pos - old_pos);
        {{/gravitational_time_dilation}}
        {{/light_travel_time}}
        {{/kerr_full_core}}

        ray = pos-old_pos;
        float solid_isec_t = 2.0;
        float ray_l = length(ray);

        {{#light_travel_time}}
        {{#gravitational_time_dilation}}
        float mix = smooth_step(1.0/u, 8.0);
        dt = mix*ray_l + (1.0-mix)*dt;
        {{/gravitational_time_dilation}}
        {{^gravitational_time_dilation}}
        dt = ray_l;
        {{/gravitational_time_dilation}}
        {{/light_travel_time}}

        {{#planetEnabled}}
        if (
            (
                old_pos.z * pos.z < 0.0 ||
                min(abs(old_pos.z), abs(pos.z)) < PLANET_RADIUS
            ) &&
            max(u, old_u) > 1.0/(PLANET_DISTANCE+PLANET_RADIUS) &&
            min(u, old_u) < 1.0/(PLANET_DISTANCE-PLANET_RADIUS)
        ) {

            {{#light_travel_time}}
            float planet_ang0 = t * PLANET_ORBITAL_ANG_VEL;
            vec3 planet_pos0 = vec3(cos(planet_ang0), sin(planet_ang0), 0)*PLANET_DISTANCE;
            {{/light_travel_time}}

            vec4 planet_isec = planet_intersection(old_pos, ray, t, dt,
                    planet_pos0, ray_doppler_factor);
            if (planet_isec.w > 0.0) {
                solid_isec_t = planet_isec.w;
                planet_isec.w = 1.0;
                color += planet_isec;
            }
        }
        {{/planetEnabled}}

        {{#accretion_thin_disk}}
        if (old_pos.z * pos.z < 0.0) {
            // crossed plane z=0

            float acc_isec_t = -old_pos.z / ray.z;
            if (acc_isec_t < solid_isec_t) {
                vec3 isec = old_pos + ray*acc_isec_t;

                float r = length(isec);
                if (r > ACCRETION_MIN_R && r < ACCRETION_MIN_R + ACCRETION_WIDTH) {
                    float angle = atan(isec.x, isec.y);

                    float temperature = accretion_temperature(r) * gravitational_shift(r);
                    float turbulence = accretion_emissivity(r, angle, t);
                    float inner_glow = exp(-8.0 * (r - ACCRETION_MIN_R));
                    float accretion_intensity = ACCRETION_BRIGHTNESS * accretion_flux_profile(r) *
                        turbulence * (1.0 + 0.7*inner_glow);

                    vec3 accretion_v;
                    {{#kerr_fast_mode}}
                    accretion_v = vec3(-isec.y, isec.x, 0.0) / (r * sqrt(2.0*(r-1.0)));
                    {{/kerr_fast_mode}}
                    {{#kerr_full_core}}
                    float rg_r = max(2.0*r, 1.0002); // convert r_s units to M units
                    float a_M = bh_rotation_enabled * bh_spin;
                    float omega_M = 1.0 / (pow(rg_r, 1.5) + a_M);
                    float v_phi = clamp(rg_r * omega_M, -0.995, 0.995);
                    vec2 xy = vec2(isec.x, isec.y);
                    float cyl_r = max(length(xy), 1e-4);
                    vec3 e_phi = vec3(-xy.y/cyl_r, xy.x/cyl_r, 0.0);
                    accretion_v = e_phi * v_phi;
                    {{/kerr_full_core}}

                    gamma = 1.0/sqrt(max(1.0-dot(accretion_v,accretion_v), 0.0001));
                    float doppler_factor = gamma*(1.0+dot(ray/ray_l,accretion_v));
                    float transfer_factor = max(ray_doppler_factor*doppler_factor, 0.05);
                    {{#beaming}}
                    {{#physical_beaming}}
                    // Liouville invariant: I_obs = D³ · I_emit.
                    // The spectrum texture stores chromaticity (normalized color),
                    // so the temperature shift only handles hue — D³ intensity
                    // boost must always be applied explicitly.
                    accretion_intensity /= pow(clamp(transfer_factor, 0.05, 20.0), 3.0);
                    {{/physical_beaming}}
                    {{^physical_beaming}}
                    // Cinematic beaming: softened for artistic rendering
                    float clamped_doppler = clamp(transfer_factor, 0.62, 1.48);
                    accretion_intensity /= pow(clamped_doppler, 1.05 + 1.10*doppler_boost);
                    {{/physical_beaming}}
                    {{/beaming}}
                    {{#doppler_shift}}
                    temperature /= transfer_factor;
                    {{/doppler_shift}}

                    vec4 thermal_color = BLACK_BODY_COLOR(temperature);
                    accretion_intensity *= look_disk_gain;
                    accretion_intensity *= 1.0 + look_glow * (0.22 + 0.78*inner_glow);
                    color += vec4(thermal_color.rgb * accretion_intensity, 1.0);
                }
            }
        }
        {{/accretion_thin_disk}}

        {{#accretion_thick_torus}}
        {
            // ADAF/RIAF volumetric emission: optically thin radiative transfer
            // dI/ds = j - alpha*I  (emission minus absorption)
            float torus_j = torus_local_emissivity(pos);
            if (torus_j > 0.001 && vol_transmittance > 0.005) {
                float path_len = length(pos - old_pos);
                float r3d = max(length(pos), 1.001);
                float cyl_r_t = max(length(pos.xy), 1e-4);
                float angle_t = atan(pos.x, pos.y);

                float temperature_t = torus_temperature(r3d) * gravitational_shift(r3d);
                float turbulence_t = accretion_turbulence(cyl_r_t, angle_t, t);

                // Emission coefficient: optically thin ADAF has low emissivity
                float j_eff = ACCRETION_BRIGHTNESS * 0.05 * torus_j * turbulence_t;
                // Absorption: very small for optically thin bremsstrahlung plasma
                float alpha_abs = 0.012 * torus_j;
                float tau_step = alpha_abs * path_len;
                float step_T = exp(-tau_step);

                // Sub-Keplerian gas velocity (ADAF: ~50% of Keplerian)
                vec3 accretion_v_t;
                {{#kerr_fast_mode}}
                float v_kep_t = 1.0 / sqrt(2.0 * max(cyl_r_t - 1.0, 0.01));
                float v_sub_t = clamp(0.5 * v_kep_t, 0.0, 0.95);
                accretion_v_t = vec3(-pos.y, pos.x, 0.0) / cyl_r_t * v_sub_t;
                {{/kerr_fast_mode}}
                {{#kerr_full_core}}
                float rg_cyl_t = max(2.0 * cyl_r_t, 1.0002);
                float a_M_t = bh_rotation_enabled * bh_spin;
                float omega_sub_t = 0.5 / (pow(rg_cyl_t, 1.5) + a_M_t);
                float v_phi_t = clamp(rg_cyl_t * omega_sub_t, -0.95, 0.95);
                vec3 e_phi_t = vec3(-pos.y, pos.x, 0.0) / cyl_r_t;
                accretion_v_t = e_phi_t * v_phi_t;
                {{/kerr_full_core}}

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

        {{#accretion_slim_disk}}
        {
            // Slim disk volumetric emission: optically thick radiative transfer
            // Super-Eddington flow: high absorption → surface-like rendering
            float slim_j = slim_disk_local_emissivity(pos);
            if (slim_j > 0.001 && vol_transmittance > 0.005) {
                float path_len_s = length(pos - old_pos);
                float cyl_r_s = max(length(pos.xy), 1e-4);
                float r3d_s = max(length(pos), 1.001);
                float angle_s = atan(pos.x, pos.y);

                float temperature_s = slim_disk_temperature(cyl_r_s) * gravitational_shift(r3d_s);
                float turbulence_s = accretion_turbulence(cyl_r_s, angle_s, t);

                // Emission coefficient: slim disk is bright (super-Eddington luminosity)
                float j_eff_s = ACCRETION_BRIGHTNESS * 0.9 * slim_j * turbulence_s;
                // Moderate absorption: optically thick → surface-like but not completely opaque
                float alpha_abs_s = 0.8 * slim_j;
                float tau_step_s = alpha_abs_s * path_len_s;
                float step_T_s = exp(-tau_step_s);
                // Source function: emission that survives self-absorption in this step
                // S = j/alpha; contribution = (1 - exp(-tau)) * S = (1-step_T) * j/alpha
                float source_contrib = (tau_step_s > 0.001)
                    ? (1.0 - step_T_s) * (j_eff_s / alpha_abs_s)
                    : j_eff_s * path_len_s;

                // Keplerian velocity outside ISCO, capped inside
                vec3 accretion_v_s;
                {{#kerr_fast_mode}}
                float v_kep_s = 1.0 / sqrt(2.0 * max(cyl_r_s - 1.0, 0.01));
                accretion_v_s = vec3(-pos.y, pos.x, 0.0) / cyl_r_s * clamp(v_kep_s, 0.0, 0.95);
                {{/kerr_fast_mode}}
                {{#kerr_full_core}}
                float rg_cyl_s = max(2.0 * cyl_r_s, 1.0002);
                float a_M_s = bh_rotation_enabled * bh_spin;
                float omega_s = 1.0 / (pow(rg_cyl_s, 1.5) + a_M_s);
                float v_phi_s = clamp(rg_cyl_s * omega_s, -0.95, 0.95);
                vec3 e_phi_s = vec3(-pos.y, pos.x, 0.0) / cyl_r_s;
                accretion_v_s = e_phi_s * v_phi_s;
                {{/kerr_full_core}}

                gamma = 1.0/sqrt(max(1.0-dot(accretion_v_s,accretion_v_s), 0.0001));
                float doppler_factor_s = gamma*(1.0+dot(ray/ray_l,accretion_v_s));
                float transfer_factor_s = max(ray_doppler_factor*doppler_factor_s, 0.05);

                float slim_intensity = source_contrib;
                {{#beaming}}
                {{#physical_beaming}}
                // Liouville D³: chromaticity texture doesn't encode brightness,
                // so intensity boost must always be applied explicitly.
                slim_intensity /= pow(clamp(transfer_factor_s, 0.05, 20.0), 3.0);
                {{/physical_beaming}}
                {{^physical_beaming}}
                float cd_s = clamp(transfer_factor_s, 0.62, 1.48);
                slim_intensity /= pow(cd_s, 1.05 + 1.10*doppler_boost);
                {{/physical_beaming}}
                {{/beaming}}
                {{#doppler_shift}}
                temperature_s /= transfer_factor_s;
                {{/doppler_shift}}

                vec4 thermal_color_s = BLACK_BODY_COLOR(temperature_s);
                slim_intensity *= look_disk_gain;

                // Front-to-back compositing with Beer-Lambert absorption
                color.rgb += vol_transmittance * thermal_color_s.rgb * slim_intensity;
                vol_transmittance *= step_T_s;
            }
        }
        {{/accretion_slim_disk}}

        {{#jet_enabled}}
        {
            float path_len_j = length(pos - old_pos);
            for (int jet_side = 0; jet_side < 2; jet_side++) {
                float sign_z = (jet_side == 0) ? 1.0 : -1.0;

                {{#jet_physical}}
                // --- Counter-jet disk occultation ---
                // The accretion disk (optically thick near equator) blocks
                // the receding counter-jet when viewed near edge-on.
                // Model as extra absorption for the counter-jet passing
                // through the disk plane.
                float occultation = 1.0;
                if (sign_z < 0.0) {
                    // Counter-jet: check if line of sight passes through
                    // a dense region near the equatorial plane
                    float abs_z_pos = abs(pos.z);
                    float cyl_r_occ = length(pos.xy);
                    // Disk is optically thick for r > ISCO, |z| < H(r)
                    // H/R ~ 0.1-0.3 for thin/slim disks
                    float disk_H = 0.2 * cyl_r_occ;
                    if (abs_z_pos < disk_H && cyl_r_occ > 1.5 && cyl_r_occ < 15.0) {
                        // Strong absorption in disk midplane
                        float tau_disk = 8.0 * (1.0 - abs_z_pos / disk_H);
                        occultation = exp(-tau_disk);
                    }
                }
                {{/jet_physical}}

                float j_jet = jet_emissivity(pos, sign_z);
                if (j_jet > 0.001 && vol_transmittance > 0.005) {
                    // Jet bulk velocity
                    vec3 v_jet = jet_velocity(pos, sign_z);
                    float v2_jet = dot(v_jet, v_jet);
                    float gamma_jet = 1.0 / sqrt(max(1.0 - v2_jet, 0.0001));

                    // Doppler factor: δ = 1 / (Γ(1 - β·n̂))
                    vec3 ray_dir = ray / max(ray_l, 1e-6);
                    float doppler_jet = gamma_jet * (1.0 + dot(ray_dir, v_jet));

                    // Synchrotron beaming: I_obs = δ^(2+α) * I_em
                    // α ≈ 0.7 for typical synchrotron spectral index
                    float beam_factor = 1.0 / pow(max(doppler_jet, 0.02), 2.7);

                    {{#jet_simple}}
                    // Simple mode: single-temperature blackbody approximation
                    float r_jet = max(length(pos), 1.001);
                    float g_shift_j = gravitational_shift(r_jet);
                    float jet_T = 18000.0 * g_shift_j;
                    {{#doppler_shift}}
                    jet_T /= max(doppler_jet, 0.05);
                    {{/doppler_shift}}
                    vec4 jet_color = BLACK_BODY_COLOR(jet_T);

                    float alpha_jet = 0.005 * j_jet;
                    float tau_jet = alpha_jet * path_len_j;
                    float step_T_jet = exp(-tau_jet);

                    float jet_emit = jet_brightness * 0.35 * j_jet * beam_factor * path_len_j;
                    jet_emit *= look_disk_gain;

                    color.rgb += vol_transmittance * jet_color.rgb * jet_emit;
                    vol_transmittance *= step_T_jet;
                    {{/jet_simple}}

                    {{#jet_physical}}
                    // Physical mode: GRMHD-calibrated emission
                    float r_jet_p = max(length(pos), 1.001);
                    float z_abs = abs(pos.z * sign_z);
                    float g_shift_j = gravitational_shift(r_jet_p);
                    float sigma = magnetization(z_abs);

                    // Effective temperature with synchrotron aging
                    float jet_T = jet_effective_temperature(z_abs, r_jet_p, sigma) * g_shift_j;

                    // Doppler temperature shift for approaching/receding jet
                    // (synchrotron is non-thermal, but this approximates the
                    // spectral shift in our blackbody color mapping)
                    jet_T /= max(doppler_jet, 0.05);

                    vec4 jet_color = BLACK_BODY_COLOR(jet_T);

                    // Absorption: higher in sheath (denser), lower in spine
                    // Synchrotron self-absorption: α_ν ∝ n_e B^(α+2) ν^(-(α+4)/2)
                    // For σ >> 1 (spine): low density → low absorption
                    // For σ ~ 1 (sheath): higher density → more absorption
                    float sigma_abs_factor = 1.0 / (1.0 + 0.5 * sigma);
                    float alpha_jet = 0.012 * j_jet * sigma_abs_factor;
                    float tau_jet = alpha_jet * path_len_j;
                    float step_T_jet = exp(-tau_jet);

                    // Emissivity: GRMHD calibration with σ-dependent
                    // efficiency. Low σ (sheath) emits more synchrotron;
                    // high σ (spine) emits less but is more beamed.
                    float sigma_emit = 1.0 / (1.0 + 0.3 * sigma);
                    float jet_emit = jet_brightness * 0.5 * j_jet * beam_factor
                                   * path_len_j * sigma_emit;
                    jet_emit *= look_disk_gain;

                    // Apply counter-jet disk occultation
                    jet_emit *= occultation;

                    color.rgb += vol_transmittance * jet_color.rgb * jet_emit;
                    vol_transmittance *= step_T_jet;
                    {{/jet_physical}}
                }
            }
        }
        {{/jet_enabled}}

        {{#light_travel_time}}
        t -= dt;
        {{/light_travel_time}}

        if (solid_isec_t <= 1.0) u = 2.0; // break

        {{#kerr_fast_mode}}
        float capture_u = 1.0 + bh_rotation_enabled * bh_spin * bh_spin_strength *
            spin_alignment * 0.12;
        capture_u = clamp(capture_u, 0.82, 1.18);
        if (u > capture_u) {
            shadow_capture = true;
            break;
        }
        {{/kerr_fast_mode}}
        {{#kerr_full_core}}
        // Realtime Kerr: capture at horizon (u >= 1.0 already handled above)
        {{/kerr_full_core}}
    }

    // the event horizon is at u = 1
    if (!shadow_capture && u < 1.0) {
        ray = normalize(pos - old_pos);
        vec2 tex_coord = sphere_map(ray * BG_COORDS);
        float t_coord;

        vec4 star_color = texture2D(star_texture, tex_coord);
        if (star_color.r > 0.0) {
            t_coord = (STAR_MIN_TEMPERATURE +
                (STAR_MAX_TEMPERATURE-STAR_MIN_TEMPERATURE) * star_color.g)
                 / ray_doppler_factor;

            color += BLACK_BODY_COLOR(t_coord) * star_color.r * STAR_BRIGHTNESS * look_star_gain * vol_transmittance;
        }

        color += galaxy_color(tex_coord, ray_doppler_factor) * GALAXY_BRIGHTNESS * look_galaxy_gain * vol_transmittance;
    }

    return color*ray_intensity;
}

void main() {

    {{#planetEnabled}}
    // "constants" derived from uniforms
    PLANET_RADIUS = planet_radius;
    PLANET_DISTANCE = max(planet_distance,planet_radius+1.5);
    float planet_orbital_v = 1.0 / sqrt(2.0*(PLANET_DISTANCE-1.0));
    PLANET_ORBITAL_ANG_VEL = -planet_orbital_v *
        sqrt(max(1.0 - 1.0/PLANET_DISTANCE, 0.0)) / PLANET_DISTANCE;
    float MAX_PLANET_ROT = max((1.0 + PLANET_ORBITAL_ANG_VEL*PLANET_DISTANCE) / PLANET_RADIUS,0.0);
    PLANET_ROTATION_ANG_VEL = -PLANET_ORBITAL_ANG_VEL + MAX_PLANET_ROT * 0.5;
    PLANET_GAMMA = 1.0/sqrt(1.0-SQ(PLANET_ORBITAL_ANG_VEL*PLANET_DISTANCE));
    {{/planetEnabled}}

    vec4 accumulated = vec4(0.0, 0.0, 0.0, 0.0);

    for (int sample_index = 0; sample_index < SAMPLE_COUNT; sample_index++) {
        vec2 jitter = sample_offset(sample_index, gl_FragCoord.xy);
        vec2 p = -1.0 + 2.0 * (gl_FragCoord.xy + jitter) / resolution.xy;
        p.y *= resolution.y / resolution.x;

        vec3 ray = normalize(p.x*cam_x + p.y*cam_y + FOV_MULT*cam_z);
        accumulated += trace_ray(ray);
    }

    vec4 color = accumulated / float(SAMPLE_COUNT);
    gl_FragColor = finalize_color(color);
}
