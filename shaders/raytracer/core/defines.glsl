// Role: Global macro definitions, uniform declarations, and scene constants used
//       across all shader stages. Loaded first in the shader concatenation order.

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
uniform vec2 cam_pan;
uniform vec2 taa_jitter;

uniform float interior_mode; // 0.0 = exterior (capture at horizon), 1.0 = interior (trace past horizon)

uniform float planet_distance, planet_radius;
uniform float disk_temperature;
uniform float bh_spin, bh_spin_strength, bh_rotation_enabled;
uniform float accretion_inner_r;
uniform float look_exposure, look_disk_gain, look_glow, look_doppler_boost;
uniform float look_aberration_strength;
uniform float look_star_gain, look_galaxy_gain;
uniform float look_tonemap_mode;
uniform float torus_r0, torus_h_ratio, torus_radial_falloff, torus_opacity, torus_outer_radius;
uniform float slim_h_ratio, slim_opacity, slim_puff_factor;
uniform float jet_half_angle, jet_lorentz, jet_brightness, jet_length;
{{#jet_physical}}
uniform float jet_magnetization, jet_knot_spacing, jet_corona_brightness;
uniform float jet_base_width, jet_corona_extent;
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
