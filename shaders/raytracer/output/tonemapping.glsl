// Role: Tone-mapping operators and final color output. Provides ACES Filmic,
//       AGX (with punchy look), Logarithmic Scientific (false-color via Inferno),
//       plus screen-space dithering and the finalize_color() compositor.

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
