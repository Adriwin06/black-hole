// Role: GLSL entry point. Drives the multi-sample anti-aliasing loop — for
//       each sample it computes a jittered ray, calls trace_ray(), accumulates
//       the result, then passes the averaged colour to finalize_color().
//       Also initialises the mutable planet orbital constants from uniforms
//       (PLANET_RADIUS, PLANET_DISTANCE, etc.) before any ray is traced.

// Halton-like low-discrepancy jitter pattern for MSAA sub-pixel sampling.
// Returns a sub-pixel offset in [−0.5, 0.5]² for sample index i.
vec2 sample_offset(int i, vec2 pixel) {
    if (SAMPLE_COUNT <= 1) return vec2(0.0, 0.0);
    float fi = float(i);
    float sample_count_f = float(SAMPLE_COUNT);
    float radius = 0.5 * sqrt((fi + 0.5) / sample_count_f);
    float base_angle = 2.0*M_PI*fract(fi*0.61803398875 + hash12(pixel*0.5));
    return radius * vec2(cos(base_angle), sin(base_angle));
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
        vec2 p_cam = p + cam_pan;

        vec3 ray = normalize(p_cam.x*cam_x + p_cam.y*cam_y + FOV_MULT*cam_z);
        accumulated += trace_ray(ray);
    }

    vec4 color = accumulated / float(SAMPLE_COUNT);
    gl_FragColor = finalize_color(color);
}
