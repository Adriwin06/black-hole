// Role: Kerr metric helpers and photon geodesic integration. Implements the
//       Schwarzschild/Kerr Binet equation (u = 1/r) with optional RK4 stepping
//       and frame-dragging from spin.

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
