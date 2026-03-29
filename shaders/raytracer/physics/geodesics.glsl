// Role: Kerr metric helpers and photon geodesic integration.
//       Two integration modes:
//         1. Schwarzschild Binet equation (u = 1/r) — exact for a = 0.
//         2. Experimental Kerr geodesic helpers in Mino time — Carter (1968)
//            separated equations for the full Kerr metric. These helpers are
//            present in the codebase but not exposed in the public UI.

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

// ═══════════════════════════════════════════════════════════════════
// 1. Schwarzschild Binet equation (used when spin ≈ 0 or interior)
// ═══════════════════════════════════════════════════════════════════

float geodesic_accel(float u, float spin_alignment) {
    float schwarzschild_accel = -u + 1.5*u*u;
    float u_drag = min(u, 1.2);
    float frame_drag_term = bh_rotation_enabled * bh_spin * bh_spin_strength *
        spin_alignment * 0.8 * u_drag*u_drag*u_drag;
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
    du += 0.5 * geodesic_accel(u, spin_alignment) * step;
    u  += du * step;
    du += 0.5 * geodesic_accel(u, spin_alignment) * step;
    {{/rk4_integration}}
}

// ═══════════════════════════════════════════════════════════════════
// 2. True Kerr geodesics — Mino-time second-order equations
// ═══════════════════════════════════════════════════════════════════
//
// Carter (1968) showed the Kerr geodesic equations separate into
// decoupled ODEs in Mino time σ (where dλ = Σ dσ):
//
//   d²r/dσ²      = R'(r)/2
//   d²(cosθ)/dσ² = Θ̃'(cosθ)/2
//   dφ/dσ        = Φ_r(r) + Φ_θ(θ)
//
// R(r) = P² − Δ K,   P = r²+a²−aξ,  K = (ξ−a)²+η
// Θ̃(c) = η(1−c²) + a²c²(1−c²) − ξ²c²   [c = cosθ]
// Constants: ξ = L_z/E,  η = Q/E²  (with E normalised to 1).
//
// These ODEs are polynomials — cheap to evaluate per step and
// produce the exact Kerr shadow (D-shaped for a > 0).

// Mino-time radial acceleration:
//   d²r/dσ² = R'(r)/2 = 2rP − (r − M)K,  M = 0.5
float kerr_r_accel(float r, float a, float xi, float eta) {
    float a2 = a * a;
    float P = r*r + a2 - a*xi;
    float K = (xi - a)*(xi - a) + eta;
    return 2.0*r*P - (r - 0.5)*K;
}

// Mino-time polar acceleration:
//   d²c/dσ² = Θ̃'(c)/2 = −c(η+ξ²−a²) − 2a²c³
float kerr_cth_accel(float cth, float a, float xi, float eta) {
    float a2 = a * a;
    return -cth*(eta + xi*xi - a2) - 2.0*a2*cth*cth*cth;
}

// Mino-time azimuthal velocity:
//   dφ/dσ = aP/Δ + ξ/sin²θ − a
float kerr_phi_dot(float r, float cth, float a, float xi) {
    float a2 = a * a;
    float Delta = max(r*r - r + a2, 0.001); // clamp near horizon
    float P = r*r + a2 - a*xi;
    float sth2 = max(1.0 - cth*cth, 1e-10);
    return a*P/Delta + xi/sth2 - a;
}

// Initialise Kerr constants of motion and Mino-time state from the
// camera position and ray direction.  Uses flat-space expressions for
// ξ and η (valid when the observer is far from the black hole) and
// projects the initial momenta onto the Kerr constraint surface
// (dr/dσ)² = R(r), (dcosθ/dσ)² = Θ̃(cosθ) so that H = 0 exactly.
void kerr_init(vec3 pos, vec3 dir, float a,
        out float xi, out float eta,
        out float r, out float cth, out float phi,
        out float pr, out float pcth) {

    float a2 = a * a;

    // Solve for Boyer-Lindquist r from Cartesian position:
    // x²+y²+z² = r² + a²sin²θ  →  r⁴ − (ρ²−a²)r² − a²z² = 0
    float rho2 = dot(pos, pos);
    float z2 = pos.z * pos.z;
    float w = rho2 - a2;
    r = sqrt(max((w + sqrt(max(w*w + 4.0*a2*z2, 0.0))) * 0.5, 0.001));

    cth = pos.z / r;
    phi = atan(pos.y, pos.x);

    float c2 = cth * cth;
    float sth2 = max(1.0 - c2, 1e-10);
    float sth  = sqrt(sth2);

    // Flat-space velocity components
    float rdot = dot(pos, dir) / r;             // dr/dλ
    float cdot = (dir.z - cth * rdot) / r;      // d(cosθ)/dλ

    // L_z = (pos × dir)·ẑ  (exact in flat space with E = 1)
    xi = pos.x * dir.y - pos.y * dir.x;

    // Carter constant Q = p_θ² + cos²θ (ξ²/sin²θ − a²)
    float theta_dot = -cdot / max(sth, 1e-6);   // dθ/dλ
    float p_theta   = r * r * theta_dot;         // p_θ = r²θ̇
    eta = p_theta*p_theta + xi*xi*c2/sth2 - a2*c2;

    // Project momenta onto constraint surface for exact energy conservation
    float P = r*r + a2 - a*xi;
    float K = (xi - a)*(xi - a) + eta;
    float Delta = r*r - r + a2;
    float R_val = P*P - max(Delta, 0.0)*K;
    pr = (rdot < 0.0) ? -sqrt(max(R_val, 0.0)) : sqrt(max(R_val, 0.0));

    float Theta_tilde = eta*(1.0-c2) + a2*c2*(1.0-c2) - xi*xi*c2;
    pcth = (cdot < 0.0) ? -sqrt(max(Theta_tilde, 0.0)) : sqrt(max(Theta_tilde, 0.0));
}

// Single Kerr integration step (Leapfrog or RK4).
void integrate_kerr_step(inout float r, inout float cth, inout float phi,
        inout float pr, inout float pcth,
        float step, float a, float xi, float eta) {
    {{#rk4_integration}}
    // ── RK4 for Kerr ───────────────────────────────────────────────
    float h = step;
    float k1_r = pr, k1_c = pcth;
    float k1_pr  = kerr_r_accel(r, a, xi, eta);
    float k1_pc  = kerr_cth_accel(cth, a, xi, eta);
    float k1_phi = kerr_phi_dot(r, cth, a, xi);

    float r2  = r + 0.5*h*k1_r;
    float c2_ = clamp(cth + 0.5*h*k1_c, -0.9999, 0.9999);
    float k2_r = pr + 0.5*h*k1_pr;
    float k2_c = pcth + 0.5*h*k1_pc;
    float k2_pr  = kerr_r_accel(r2, a, xi, eta);
    float k2_pc  = kerr_cth_accel(c2_, a, xi, eta);
    float k2_phi = kerr_phi_dot(r2, c2_, a, xi);

    float r3  = r + 0.5*h*k2_r;
    float c3_ = clamp(cth + 0.5*h*k2_c, -0.9999, 0.9999);
    float k3_r = pr + 0.5*h*k2_pr;
    float k3_c = pcth + 0.5*h*k2_pc;
    float k3_pr  = kerr_r_accel(r3, a, xi, eta);
    float k3_pc  = kerr_cth_accel(c3_, a, xi, eta);
    float k3_phi = kerr_phi_dot(r3, c3_, a, xi);

    float r4  = r + h*k3_r;
    float c4_ = clamp(cth + h*k3_c, -0.9999, 0.9999);
    float k4_r = pr + h*k3_pr;
    float k4_c = pcth + h*k3_pc;
    float k4_pr  = kerr_r_accel(r4, a, xi, eta);
    float k4_pc  = kerr_cth_accel(c4_, a, xi, eta);
    float k4_phi = kerr_phi_dot(r4, c4_, a, xi);

    r    += (h/6.0)*(k1_r   + 2.0*k2_r   + 2.0*k3_r   + k4_r);
    cth  += (h/6.0)*(k1_c   + 2.0*k2_c   + 2.0*k3_c   + k4_c);
    pr   += (h/6.0)*(k1_pr  + 2.0*k2_pr  + 2.0*k3_pr  + k4_pr);
    pcth += (h/6.0)*(k1_pc  + 2.0*k2_pc  + 2.0*k3_pc  + k4_pc);
    phi  += (h/6.0)*(k1_phi + 2.0*k2_phi + 2.0*k3_phi + k4_phi);
    cth   = clamp(cth, -0.9999, 0.9999);
    {{/rk4_integration}}
    {{^rk4_integration}}
    // ── Leapfrog / Störmer-Verlet (symplectic 2nd-order) ───────────
    pr   += 0.5 * step * kerr_r_accel(r, a, xi, eta);
    pcth += 0.5 * step * kerr_cth_accel(cth, a, xi, eta);
    r   += step * pr;
    cth += step * pcth;
    cth  = clamp(cth, -0.9999, 0.9999);
    phi += step * kerr_phi_dot(r, cth, a, xi);
    pr   += 0.5 * step * kerr_r_accel(r, a, xi, eta);
    pcth += 0.5 * step * kerr_cth_accel(cth, a, xi, eta);
    {{/rk4_integration}}
}
