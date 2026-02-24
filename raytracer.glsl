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

uniform sampler2D galaxy_texture, star_texture,
    planet_texture, spectrum_texture;

// stepping and anti-aliasing parameters
const int NSTEPS = {{n_steps}};
const int SAMPLE_COUNT = {{sample_count}};
const float MAX_REVOLUTIONS = float({{max_revolutions}});

const float ACCRETION_MIN_R = 1.5;
const float ACCRETION_WIDTH = 5.0;
const float ACCRETION_BRIGHTNESS = 1.2;

const float STAR_MIN_TEMPERATURE = 4000.0;
const float STAR_MAX_TEMPERATURE = 15000.0;

const float STAR_BRIGHTNESS = 0.7;
const float GALAXY_BRIGHTNESS = 0.22;
const float GLOBAL_EXPOSURE = 0.72;
const float GALAXY_DOPPLER_STRENGTH = 0.45;
const float GALAXY_MAX_BOOST = 1.30;

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

float geodesic_accel(float u, float spin_alignment) {
    float schwarzschild_accel = -u*(1.0 - 1.5*u*u);
    float frame_drag_term = bh_rotation_enabled * bh_spin * bh_spin_strength *
        spin_alignment * 0.55 * u*u*u*u;
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

float accretion_emissivity(float radius, float angle, float t) {
    float r_norm = (radius - ACCRETION_MIN_R) / ACCRETION_WIDTH;
    float edge_fade = smoothstep(0.02, 0.18, r_norm) *
        (1.0 - smoothstep(0.78, 1.0, r_norm));

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

    return max((0.45 + 1.1*plasma) * (0.8 + 0.2*filaments) * edge_fade, 0.02);
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
    vec3 mapped = aces_filmic(max(color.rgb * GLOBAL_EXPOSURE, vec3(0.0)));
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

    vec4 light_color = BLACK_BODY_COLOR(light_temperature);
    ret = texture2D(planet_texture, tex_coord) * lightness * light_color;
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
    ray = lorentz_velocity_transformation(ray, cam_vel);
    {{/aberration}}

    float ray_intensity = 1.0;
    float ray_doppler_factor = 1.0;

    float gamma = 1.0/sqrt(max(1.0-dot(cam_vel,cam_vel), 0.0001));
    ray_doppler_factor = gamma*(1.0 + dot(ray,-cam_vel));
    {{#beaming}}
    float beaming_factor = clamp(ray_doppler_factor, 0.78, 1.22);
    ray_intensity /= beaming_factor*beaming_factor*beaming_factor;
    {{/beaming}}
    {{^doppler_shift}}
    ray_doppler_factor = 1.0;
    {{/doppler_shift}}

    float step = 0.01;
    vec4 color = vec4(0.0,0.0,0.0,1.0);

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

        step = MAX_REVOLUTIONS * 2.0*M_PI / float(NSTEPS);

        // adaptive step size
        float max_rel_u_change = (1.0-log(max(u, 0.0001)))*10.0 / float(NSTEPS);
        if ((du > 0.0 || (du0 < 0.0 && u0/u < 5.0)) && abs(du) > abs(max_rel_u_change*u) / step) {
            step = max_rel_u_change*u/abs(du);
        }

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

        {{#accretion_disk}}
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

                    vec3 accretion_v = vec3(-isec.y, isec.x, 0.0) / sqrt(2.0*(r-1.0)) / (r*r);
                    gamma = 1.0/sqrt(max(1.0-dot(accretion_v,accretion_v), 0.0001));
                    float doppler_factor = gamma*(1.0+dot(ray/ray_l,accretion_v));
                    {{#beaming}}
                    float clamped_doppler = max(doppler_factor, 0.05);
                    accretion_intensity /= clamped_doppler*clamped_doppler*clamped_doppler;
                    {{/beaming}}
                    {{#doppler_shift}}
                    temperature /= max(ray_doppler_factor*doppler_factor, 0.05);
                    {{/doppler_shift}}

                    vec4 thermal_color = BLACK_BODY_COLOR(temperature);
                    color += vec4(thermal_color.rgb * accretion_intensity, 1.0);
                }
            }
        }
        {{/accretion_disk}}

        {{#light_travel_time}}
        t -= dt;
        {{/light_travel_time}}

        if (solid_isec_t <= 1.0) u = 2.0; // break

        float capture_u = 1.0 + bh_rotation_enabled * bh_spin * bh_spin_strength *
            spin_alignment * 0.12;
        capture_u = clamp(capture_u, 0.82, 1.18);
        if (u > capture_u) {
            shadow_capture = true;
            break;
        }
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

            color += BLACK_BODY_COLOR(t_coord) * star_color.r * STAR_BRIGHTNESS;
        }

        color += galaxy_color(tex_coord, ray_doppler_factor) * GALAXY_BRIGHTNESS;
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
