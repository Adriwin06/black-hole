// Role: Render quality preset library — centralizes quality tier parameter
//       values and selector labels. Applied by setupGUI() in gui.js so quality
//       behavior is data-driven and separated from UI wiring.

/*global QUALITY_PRESETS:true, QUALITY_PRESET_LABELS:true, KERR_MODE_LABELS:true, applyQualityPresetValues:true */
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

function applyQualityPresetValues(parameters, presetName) {
    if (!parameters) return null;

    var preset = QUALITY_PRESETS[presetName];
    if (!preset) return null;

    var isKerr = (parameters.kerr_mode === 'realtime_full_kerr_core');
    var modeValues = isKerr ? preset.kerr : preset.standard;
    if (!modeValues) return null;

    parameters.quality = presetName;
    parameters.n_steps = modeValues.n_steps;
    parameters.sample_count = modeValues.sample_count;
    parameters.max_revolutions = modeValues.max_revolutions;
    parameters.rk4_integration = modeValues.rk4_integration;
    parameters.cinematic_tonemap = preset.cinematic_tonemap;
    parameters.resolution_scale = preset.resolution_scale;
    parameters.taa_enabled = preset.taa_enabled;
    parameters.taa.history_weight = preset.taa.history_weight;
    parameters.taa.clip_box = preset.taa.clip_box;
    parameters.taa.motion_rejection = preset.taa.motion_rejection;
    parameters.taa.max_camera_delta = preset.taa.max_camera_delta;
    parameters.taa.motion_clip_scale = preset.taa.motion_clip_scale;

    return preset;
}

var QUALITY_PRESET_LABELS = {
    'Custom': 'custom',
    'Mobile (fastest)': 'mobile',
    'Optimal (recommended)': 'optimal',
    'High (full resolution)': 'medium',
    'Ultra (max quality)': 'high'
};

var KERR_MODE_LABELS = {
    'Fast (approximate lensing)': 'fast',
    'Realtime Kerr core': 'realtime_full_kerr_core'
};
