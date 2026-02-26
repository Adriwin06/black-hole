// Role: Background sky rendering — Milky Way galaxy texture with relativistic
//       Doppler color shift (temperature lookup + blackbody remap) and H-alpha
//       emission handling for moving observer modes.

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
