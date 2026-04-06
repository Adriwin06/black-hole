// Role: Relativistic jet emission models. Two modes are provided:
//       - Simple: smooth parabolic jet with spine/limb synchrotron approximation.
//       - Physical: more detailed GRMHD-inspired model with magnetization
//         parameter, spine/sheath structure, reconfinement knots,
//         jet-corona connection, and disk occultation.

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
