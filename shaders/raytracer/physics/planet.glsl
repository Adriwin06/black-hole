// Role: Ray-sphere intersection and shading for the orbiting reference planet.
//       Handles Lorentz contraction, light-travel-time retardation, Doppler
//       temperature shift, and relativistic beaming of the disk illumination.

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
