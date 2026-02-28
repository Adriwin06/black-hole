// Role: Render quality preset library — centralizes quality tier parameter
//       values and selector labels. Applied by setupGUI() in gui.js so quality
//       behavior is data-driven and separated from UI wiring.

/*global QUALITY_PRESETS:true, QUALITY_PRESET_LABELS:true, KERR_MODE_LABELS:true */
var QUALITY_PRESETS = {
    mobile: {
        standard: { n_steps: 28, sample_count: 1, max_revolutions: 1.4, rk4_integration: false },
        kerr: { n_steps: 120, sample_count: 1, max_revolutions: 2.0, rk4_integration: false },
        cinematic_tonemap: true,
        resolution_scale: 0.55,
        taa_enabled: true,
        taa: {
            history_weight: 0.82,
            clip_box: 0.08,
            motion_rejection: 10.0,
            max_camera_delta: 0.07,
            motion_clip_scale: 0.8
        },
        hide_planet_controls: true
    },
    fast: {
        standard: { n_steps: 40, sample_count: 1, max_revolutions: 1.5, rk4_integration: false },
        kerr: { n_steps: 200, sample_count: 2, max_revolutions: 2.5, rk4_integration: true },
        cinematic_tonemap: true,
        resolution_scale: 1.0,
        taa_enabled: false,
        taa: {
            history_weight: 0.88,
            clip_box: 0.06,
            motion_rejection: 8.0,
            max_camera_delta: 0.08,
            motion_clip_scale: 0.6
        },
        hide_planet_controls: true
    },
    medium: {
        standard: { n_steps: 100, sample_count: 1, max_revolutions: 2.0, rk4_integration: false },
        kerr: { n_steps: 400, sample_count: 3, max_revolutions: 3.0, rk4_integration: true },
        cinematic_tonemap: true,
        resolution_scale: 1.0,
        taa_enabled: false,
        taa: {
            history_weight: 0.88,
            clip_box: 0.06,
            motion_rejection: 8.0,
            max_camera_delta: 0.08,
            motion_clip_scale: 0.6
        },
        hide_planet_controls: false
    },
    high: {
        standard: { n_steps: 320, sample_count: 4, max_revolutions: 3.2, rk4_integration: true },
        kerr: { n_steps: 520, sample_count: 4, max_revolutions: 3.5, rk4_integration: true },
        cinematic_tonemap: true,
        resolution_scale: 1.0,
        taa_enabled: false,
        taa: {
            history_weight: 0.88,
            clip_box: 0.06,
            motion_rejection: 8.0,
            max_camera_delta: 0.08,
            motion_clip_scale: 0.6
        },
        hide_planet_controls: false
    },
    optimal: {
        standard: { n_steps: 100, sample_count: 1, max_revolutions: 2.0, rk4_integration: false },
        kerr: { n_steps: 400, sample_count: 1, max_revolutions: 3.5, rk4_integration: true },
        cinematic_tonemap: true,
        resolution_scale: 0.8,
        taa_enabled: true,
        taa: {
            history_weight: 0.82,
            clip_box: 0.08,
            motion_rejection: 10.0,
            max_camera_delta: 0.07,
            motion_clip_scale: 0.8
        },
        hide_planet_controls: false
    }
};

var QUALITY_PRESET_LABELS = {
    'Custom': 'custom',
    'Mobile (low power)': 'mobile',
    'Fast (preview)': 'fast',
    'Medium': 'medium',
    'High': 'high',
    'Optimal': 'optimal'
};

var KERR_MODE_LABELS = {
    'Fast (approximate lensing)': 'fast',
    'Realtime Kerr core': 'realtime_full_kerr_core'
};
