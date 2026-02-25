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
uniform float torus_r0, torus_h_ratio;
uniform float jet_half_angle, jet_lorentz, jet_brightness, jet_length;

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

const float PLANET_AMBIENT = 0.1;
const float PLANET_LIGHTNESS = 1.5;

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

void cartesian_to_spherical(vec3 p, out float r, out float theta, out float phi) {
    r = max(length(p), 1e-5);
    theta = acos(clamp(p.z / r, -1.0, 1.0));
    phi = atan(p.y, p.x);
}

vec3 spherical_to_cartesian(float r, float theta, float phi) {
    float st = sin(theta);
    return vec3(
        r * st * cos(phi),
        r * st * sin(phi),
        r * cos(theta)
    );
}

void spherical_basis(float theta, float phi, out vec3 e_r, out vec3 e_theta, out vec3 e_phi) {
    float st = sin(theta);
    float ct = cos(theta);
    float sp = sin(phi);
    float cp = cos(phi);

    e_r = vec3(st*cp, st*sp, ct);
    e_theta = vec3(ct*cp, ct*sp, -st);
    e_phi = vec3(-sp, cp, 0.0);
}

float kerr_radial_potential(float r, float a, float Lz, float Q) {
    float Delta = kerr_delta(r, a);
    float P = r*r + a*a - a*Lz;
    float B = (Lz - a)*(Lz - a) + Q;
    return P*P - Delta*B;
}

float kerr_polar_potential(float theta, float a, float Lz, float Q) {
    float st = max(sin(theta), 1e-4);
    float ct = cos(theta);
    return Q - ct*ct * (Lz*Lz/(st*st) - a*a);
}

// Simplified Kerr ray step using coordinate velocity approach
// Uses effective potential derivative for photon geodesics
void kerr_accel(vec3 pos, vec3 vel, float a, out vec3 accel) {
    float r = max(length(pos), 0.5);  // Prevent divide by zero
    float r2 = r*r;
    float r3 = r2*r;
    float r4 = r2*r2;
    float a2 = a*a;
    float M = 0.5;  // M = rs/2 where rs = 1
    float rs = 1.0; // Schwarzschild radius
    
    vec3 r_hat = pos / r;
    vec3 z_hat = vec3(0.0, 0.0, 1.0);
    
    // Angular momentum magnitude squared
    vec3 L_vec = cross(pos, vel);
    float h2 = dot(L_vec, L_vec);  // specific angular momentum squared
    
    // From geodesic equation, radial acceleration for photons:
    // a_r = -M/r² + h²/r³ - 3Mh²/r⁴
    // The last term is the key GR correction causing photon sphere at r=3M=1.5rs
    float radial_accel = -M/r2 + h2/r3 - 3.0*M*h2/r4;
    
    accel = r_hat * radial_accel;
    
    // Transverse acceleration to maintain geodesic (angular momentum evolution)
    // For non-radial motion, there's also theta/phi acceleration
    vec3 v_radial = r_hat * dot(vel, r_hat);
    vec3 v_perp = vel - v_radial;
    float v_perp_mag = length(v_perp);
    if (v_perp_mag > 0.001) {
        vec3 theta_hat = v_perp / v_perp_mag;
        // Centripetal-like term from spherical coordinates
        float v_r = dot(vel, r_hat);
        accel -= theta_hat * v_r * v_perp_mag / r;
    }
    
    // Frame dragging for Kerr
    float cos_theta = clamp(dot(r_hat, z_hat), -1.0, 1.0);
    float sin_theta = sqrt(max(1.0 - cos_theta*cos_theta, 0.0001));
    if (abs(a) > 0.001 && sin_theta > 0.001) {
        float omega_fd = 2.0*M*a*r / (r4 + a2*r2 + 2.0*M*a2*r);
        vec3 phi_hat = normalize(cross(z_hat, pos));
        accel += phi_hat * omega_fd * sin_theta;
    }
}

void integrate_kerr_simple_step(inout vec3 pos, inout vec3 vel, float h, float a) {
    vec3 k1_v, k1_a;
    vec3 k2_v, k2_a;
    vec3 k3_v, k3_a;
    vec3 k4_v, k4_a;
    
    kerr_accel(pos, vel, a, k1_a);
    k1_v = vel;
    
    kerr_accel(pos + 0.5*h*k1_v, vel + 0.5*h*k1_a, a, k2_a);
    k2_v = vel + 0.5*h*k1_a;
    
    kerr_accel(pos + 0.5*h*k2_v, vel + 0.5*h*k2_a, a, k3_a);
    k3_v = vel + 0.5*h*k2_a;
    
    kerr_accel(pos + h*k3_v, vel + h*k3_a, a, k4_a);
    k4_v = vel + h*k3_a;
    
    pos += (h/6.0) * (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v);
    vel += (h/6.0) * (k1_a + 2.0*k2_a + 2.0*k3_a + k4_a);
    
    // Re-normalize velocity to maintain null geodesic (light speed = 1)
    vel = normalize(vel);
}

void kerr_constants_from_ray(vec3 pos, vec3 ray, float a,
        out float Lz, out float Q, out float sign_r, out float sign_theta) {

    float r, theta, phi;
    cartesian_to_spherical(pos, r, theta, phi);

    vec3 e_r, e_theta, e_phi;
    spherical_basis(theta, phi, e_r, e_theta, e_phi);

    // ray points in the direction we trace (toward scene, backward along photon path)
    float p_r = dot(ray, e_r);
    float p_theta = dot(ray, e_theta) * r;
    float p_phi = dot(ray, e_phi) * r * max(sin(theta), 1e-4);

    Lz = p_phi;

    float st = max(sin(theta), 1e-4);
    float ct = cos(theta);
    Q = p_theta*p_theta + ct*ct * (Lz*Lz/(st*st) - a*a);
    Q = max(Q, 0.0);

    sign_r = (p_r >= 0.0) ? 1.0 : -1.0;
    sign_theta = (dot(ray, e_theta) >= 0.0) ? 1.0 : -1.0;
}

void kerr_derivatives(float r, float theta, float a, float Lz, float Q,
        float sign_r, float sign_theta,
        out float dr, out float dtheta, out float dphi) {

    float Sigma = max(kerr_sigma(r, theta, a), 1e-5);
    float Delta = max(kerr_delta(r, a), 1e-5);
    float st = max(sin(theta), 1e-4);

    float R = max(kerr_radial_potential(r, a, Lz, Q), 0.0);
    float T = max(kerr_polar_potential(theta, a, Lz, Q), 0.0);
    float P = r*r + a*a - a*Lz;

    dr = sign_r * sqrt(R) / Sigma;
    dtheta = sign_theta * sqrt(T) / Sigma;
    dphi = (Lz/(st*st) - a + a*P/Delta) / Sigma;
}

void integrate_kerr_bl_step(inout float r, inout float theta, inout float phi,
        float h, float a, float Lz, float Q,
        inout float sign_r, inout float sign_theta) {

    float k1_r, k1_t, k1_p;
    float k2_r, k2_t, k2_p;
    float k3_r, k3_t, k3_p;
    float k4_r, k4_t, k4_p;

    kerr_derivatives(r, theta, a, Lz, Q, sign_r, sign_theta, k1_r, k1_t, k1_p);
    kerr_derivatives(r + 0.5*h*k1_r, theta + 0.5*h*k1_t, a, Lz, Q,
        sign_r, sign_theta, k2_r, k2_t, k2_p);
    kerr_derivatives(r + 0.5*h*k2_r, theta + 0.5*h*k2_t, a, Lz, Q,
        sign_r, sign_theta, k3_r, k3_t, k3_p);
    kerr_derivatives(r + h*k3_r, theta + h*k3_t, a, Lz, Q,
        sign_r, sign_theta, k4_r, k4_t, k4_p);

    r += (h/6.0) * (k1_r + 2.0*k2_r + 2.0*k3_r + k4_r);
    theta += (h/6.0) * (k1_t + 2.0*k2_t + 2.0*k3_t + k4_t);
    phi += (h/6.0) * (k1_p + 2.0*k2_p + 2.0*k3_p + k4_p);

    theta = clamp(theta, 1e-4, M_PI - 1e-4);

    float R = kerr_radial_potential(r, a, Lz, Q);
    float T = kerr_polar_potential(theta, a, Lz, Q);

    if (R < 1e-7) sign_r *= -1.0;
    if (T < 1e-7) sign_theta *= -1.0;
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
    float x = max(radius / ACCRETION_MIN_R, 1.0001);
    float inner_edge = max(1.0 - sqrt(1.0 / x), 0.02);
    return disk_temperature * pow(1.0 / x, 0.75) * pow(inner_edge, 0.25);
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

{{#jet_enabled}}
// Relativistic jet: Blandford-Znajek powered bipolar outflow
// Conical geometry along spin axis (±z), synchrotron emission,
// relativistic beaming with bulk Lorentz factor Γ
// Returns: emissivity at point p for a given jet cone
float jet_emissivity(vec3 p, float sign_z) {
    float z = p.z * sign_z; // flip for counter-jet
    if (z < 0.5) return 0.0; // jet starts above ~0.5 r_s from BH

    float cyl_r = length(p.xy);
    float r3d = length(p);

    // Cone half-angle in radians
    float theta_jet = jet_half_angle * DEG_TO_RAD;
    float tan_half = tan(theta_jet);

    // Cone boundary: cyl_r < z * tan(theta)
    float cone_r = z * tan_half;
    if (cyl_r > cone_r * 1.5) return 0.0;

    // Transverse Gaussian profile across cone
    float r_norm = cyl_r / max(cone_r, 0.01);
    float transverse = exp(-2.0 * r_norm * r_norm);

    // Longitudinal: synchrotron emissivity ~ z^(-2) (adiabatic expansion)
    // with smooth onset above the launching region
    float z_onset = smoothstep(0.5, 1.5, z);
    float z_decay = 1.0 / (1.0 + z * z * 0.01);

    // Length cutoff
    float z_cutoff = 1.0 - smoothstep(jet_length * 0.8, jet_length, z);

    // Limb-brightening: enhanced emission at cone boundary
    // (magnetic field compression at sheath, typical in VLBI observations)
    float limb = 1.0 + 1.5 * exp(-8.0 * (r_norm - 0.85) * (r_norm - 0.85));

    return transverse * z_onset * z_decay * z_cutoff * limb;
}

// Jet bulk velocity vector (along ±z axis with slight helical twist)
vec3 jet_velocity(vec3 p, float sign_z) {
    float beta = sqrt(1.0 - 1.0 / (jet_lorentz * jet_lorentz));
    // Predominantly along spin axis with small toroidal component
    // (helical magnetic field → helical velocity)
    float cyl_r = length(p.xy);
    float phi_twist = 0.05 * beta; // small helical fraction
    vec3 v_z = vec3(0.0, 0.0, sign_z) * beta;
    vec3 v_phi = vec3(-p.y, p.x, 0.0) / max(cyl_r, 0.01) * phi_twist;
    return v_z + v_phi;
}
{{/jet_enabled}}

vec2 sample_offset(int i, vec2 pixel) {
    if (SAMPLE_COUNT <= 1) return vec2(0.0, 0.0);
    float fi = float(i);
    float sample_count_f = float(SAMPLE_COUNT);
    float radius = 0.5 * sqrt((fi + 0.5) / sample_count_f);
    float base_angle = 2.0*M_PI*fract(fi*0.61803398875 + hash12(pixel*0.5));
    return radius * vec2(cos(base_angle), sin(base_angle));
}

vec3 aces_filmic(vec3 x) {
    return clamp((x*(2.51*x + 0.03)) / (x*(2.43*x + 0.59) + 0.14), 0.0, 1.0);
}

float screen_dither() {
    return (hash12(gl_FragCoord.xy) - 0.5) / 255.0;
}

vec4 finalize_color(vec4 color) {
    {{#cinematic_tonemap}}
    vec3 mapped = aces_filmic(max(color.rgb * (GLOBAL_EXPOSURE * look_exposure), vec3(0.0)));
    float lum = dot(mapped, vec3(0.2126, 0.7152, 0.0722));
    mapped += mapped * smoothstep(0.82, 1.0, lum) * 0.06;
    mapped = pow(mapped, vec3(1.0/2.2));
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
    float lightness = ((1.0-PLANET_AMBIENT)*diffuse + PLANET_AMBIENT) *
        PLANET_LIGHTNESS;

    float light_temperature = disk_temperature;
    {{#doppler_shift}}
    float doppler_factor = SQ(PLANET_GAMMA) *
        (1.0 + dot(planet_vel, light_dir)) *
        (1.0 - dot(planet_vel, normalize(ray)));
    light_temperature /= max(doppler_factor * ray_doppler_factor, 0.05);
    {{/doppler_shift}}

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
    float kerr_r = length(pos);
    float kerr_theta = acos(clamp(pos.z/max(kerr_r, 1e-5), -1.0, 1.0));
    float kerr_phi = atan(pos.y, pos.x);
    float kerr_a = kerr_spin_a();
    float kerr_horizon = kerr_horizon_radius(kerr_a);
    float kerr_Lz = 0.0;
    float kerr_Q = 0.0;
    float kerr_sign_r = 1.0;
    float kerr_sign_theta = 1.0;
    kerr_constants_from_ray(pos, ray, kerr_a, kerr_Lz, kerr_Q, kerr_sign_r, kerr_sign_theta);

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
        {{^kerr_offline}}
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
        kerr_r = length(pos);

        {{#light_travel_time}}
        dt = length(pos - old_pos);
        {{/light_travel_time}}
        {{/kerr_offline}}

        {{#kerr_offline}}
        // Offline: Exact Boyer-Lindquist Kerr geodesic integration
        // Uses conserved quantities Lz, Q from Carter (1968)

        float bl_dr_est, bl_dtheta_est, bl_dphi_est;
        kerr_derivatives(kerr_r, kerr_theta, kerr_a, kerr_Lz, kerr_Q,
            kerr_sign_r, kerr_sign_theta, bl_dr_est, bl_dtheta_est, bl_dphi_est);

        // Use angular rate to set step for consistent coverage
        float phi_rate = max(abs(bl_dphi_est), 0.001);
        float phi_step_bl = MAX_REVOLUTIONS * 2.0*M_PI / float(NSTEPS);
        step = phi_step_bl / phi_rate;

        // Conservative limits: never exceed 5% of current radius
        step = min(step, kerr_r * 0.05);
        step = max(step, 0.0001);

        // Refine near photon sphere and horizon
        float ps_dist_bl = kerr_r - 1.5;
        float ps_prox_bl = exp(-4.0 * ps_dist_bl * ps_dist_bl);
        step *= max(1.0 - 0.9 * ps_prox_bl, 0.05);

        float horizon_margin = max((kerr_r - kerr_horizon) / kerr_r, 0.0);
        step *= clamp(horizon_margin * 5.0, 0.01, 1.0);

        old_u = u;
        old_pos = pos;

        integrate_kerr_bl_step(kerr_r, kerr_theta, kerr_phi,
            step, kerr_a, kerr_Lz, kerr_Q,
            kerr_sign_r, kerr_sign_theta);

        if (kerr_r <= kerr_horizon + 0.005) {
            shadow_capture = true;
            break;
        }

        pos = spherical_to_cartesian(kerr_r, kerr_theta, kerr_phi);
        u = 1.0 / max(kerr_r, 0.01);

        {{#light_travel_time}}
        dt = length(pos - old_pos);
        {{/light_travel_time}}
        {{/kerr_offline}}
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
                    // For blackbody: B_ν(D·T) = D³·B_{ν/D}(T), so shifting temperature
                    // already accounts for the full D³ Liouville invariant.
                    // Only apply D³ if doppler_shift is disabled (no temperature shift).
                    {{^doppler_shift}}
                    accretion_intensity /= pow(clamp(transfer_factor, 0.05, 20.0), 3.0);
                    {{/doppler_shift}}
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
                // Temperature shift handles D³ for blackbody (no double-counting)
                {{^doppler_shift}}
                torus_intensity /= pow(clamp(transfer_factor_t, 0.05, 20.0), 3.0);
                {{/doppler_shift}}
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
                // Temperature shift handles D³ for blackbody (no double-counting)
                {{^doppler_shift}}
                slim_intensity /= pow(clamp(transfer_factor_s, 0.05, 20.0), 3.0);
                {{/doppler_shift}}
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
            // Relativistic bipolar jet: volumetric synchrotron emission
            // Process both jet (+z) and counter-jet (-z)
            float path_len_j = length(pos - old_pos);
            for (int jet_side = 0; jet_side < 2; jet_side++) {
                float sign_z = (jet_side == 0) ? 1.0 : -1.0;
                float j_jet = jet_emissivity(pos, sign_z);
                if (j_jet > 0.001 && vol_transmittance > 0.005) {
                    // Jet bulk velocity
                    vec3 v_jet = jet_velocity(pos, sign_z);
                    float v2_jet = dot(v_jet, v_jet);
                    float gamma_jet = 1.0 / sqrt(max(1.0 - v2_jet, 0.0001));

                    // Doppler factor: δ = 1 / (Γ(1 - β·n̂))
                    // where n̂ is the direction FROM emitter TO observer (= -ray_dir)
                    vec3 ray_dir = ray / max(ray_l, 1e-6);
                    float doppler_jet = gamma_jet * (1.0 + dot(ray_dir, v_jet));

                    // Synchrotron beaming: I_obs = δ^(2+α) * I_em
                    // α ≈ 0.7 for synchrotron spectral index → exponent ≈ 2.7
                    // Use δ^2 (continuous jet) for balanced side/face-on visibility
                    // Synchrotron beaming: I_obs = δ^(2+α) * I_em
                    // α ≈ 0.7 for typical synchrotron spectral index
                    float beam_factor = 1.0 / pow(max(doppler_jet, 0.02), 2.7);

                    // Synchrotron emission: approximate color as high-T blackbody
                    // No Doppler temperature shift (synchrotron is non-thermal;
                    // the D^2.7 beaming already handles the full relativistic effect)
                    float r_jet = max(length(pos), 1.001);
                    float g_shift_j = gravitational_shift(r_jet);
                    float jet_T = 25000.0 * g_shift_j;
                    vec4 jet_color = BLACK_BODY_COLOR(jet_T);

                    // Small absorption (optically thin jet)
                    float alpha_jet = 0.008 * j_jet;
                    float tau_jet = alpha_jet * path_len_j;
                    float step_T_jet = exp(-tau_jet);

                    float jet_emit = jet_brightness * 0.25 * j_jet * beam_factor * path_len_j;
                    jet_emit *= look_disk_gain;

                    color.rgb += vol_transmittance * jet_color.rgb * jet_emit;
                    vol_transmittance *= step_T_jet;
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
        {{^kerr_offline}}
        // Capture when photon reaches Kerr horizon radius
        float capture_u_fc = 1.0 / max(kerr_horizon, 0.01);
        if (u > capture_u_fc) {
            shadow_capture = true;
            break;
        }
        {{/kerr_offline}}
        {{#kerr_offline}}
        // BL integrator: capture handled above during integration
        if (kerr_r <= kerr_horizon + 0.005) {
            shadow_capture = true;
            break;
        }
        {{/kerr_offline}}
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
    PLANET_ORBITAL_ANG_VEL = -1.0 / sqrt(2.0*(PLANET_DISTANCE-1.0)) / PLANET_DISTANCE;
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
