export function setupRenderingFolder(options) {
    options = options || {};

    var gui = options.gui;
    var parameters = options.parameters;
    var addControl = options.addControl;
    var qualityPresetLabels = options.qualityPresetLabels;
    var kerrModeLabels = options.kerrModeLabels;
    var applyBlackHolePreset = options.applyBlackHolePreset;
    var applyQualityPreset = options.applyQualityPreset;
    var applyKerrMode = options.applyKerrMode;
    var updateShader = options.updateShader;
    var onResolutionScaleChange = options.onResolutionScaleChange || function() {};
    var onTaaEnabledChange = options.onTaaEnabledChange || function() {};
    var onTaaSettingChange = options.onTaaSettingChange || function() {};

    var presetObj = { preset: 'Default' };
    var presetsFolder = gui.addFolder('Black Hole Presets');
    var blackHolePresetController = addControl(presetsFolder, presetObj, 'preset', {
        name: 'preset',
        options: ['Custom', 'Default', 'M87*', 'Sgr A*', 'Cygnus X-1', 'GRS 1915+105', 'Gargantua (Interstellar visuals)', 'Schwarzschild'],
        onChange: function(value) { if (applyBlackHolePreset) applyBlackHolePreset(value); },
        help: 'Load a literature-inspired black-hole preset. These are illustrative starting points, not definitive measurements.'
    });
    presetsFolder.open();

    var renderFolder = gui.addFolder('Rendering');
    var qualityPresetController = addControl(renderFolder, parameters, 'quality', {
        options: qualityPresetLabels,
        name: 'quality preset',
        onChange: applyQualityPreset,
        help: 'Global render preset. Higher presets use more ray-marching steps and samples.'
    });
    addControl(renderFolder, parameters, 'kerr_mode', {
        options: kerrModeLabels,
        name: 'solver mode',
        onChange: applyKerrMode,
        help: 'Fast = Schwarzschild/Binet photon lensing with a perturbative frame-drag term. Kerr-inspired disk velocities = same photon lensing, but emitting matter uses Kerr angular velocity in a simplified local-speed model.'
    });
    addControl(renderFolder, parameters, 'n_steps', {
        min: 20,
        max: 1400,
        step: 1,
        name: 'ray steps',
        trackQualityPreset: true,
        onChange: updateShader,
        help: 'More steps improve thin features and strong lensing, but reduce FPS.'
    });
    addControl(renderFolder, parameters, 'sample_count', {
        min: 1,
        max: 12,
        step: 1,
        name: 'samples / pixel',
        trackQualityPreset: true,
        onChange: updateShader,
        help: 'Supersampling for anti-aliasing and smoother edges. Higher values are slower.'
    });
    addControl(renderFolder, parameters, 'max_revolutions', {
        min: 1.0,
        max: 8.0,
        step: 0.1,
        name: 'max orbit turns',
        trackQualityPreset: true,
        onChange: updateShader,
        help: 'How many wrapped photon turns are traced before escape/capture cutoff.'
    });
    addControl(renderFolder, parameters, 'resolution_scale', {
        min: 0.35,
        max: 2.0,
        step: 0.05,
        name: 'resolution scale',
        trackQualityPreset: true,
        onChange: onResolutionScaleChange,
        help: 'Internal rendering resolution multiplier. <1 lowers cost, >1 supersamples for sharper output.'
    });
    var taaEnabledCtrl = addControl(renderFolder, parameters, 'taa_enabled', {
        name: 'TAA (stable)',
        trackQualityPreset: true,
        onChange: onTaaEnabledChange,
        help: 'Temporal anti-aliasing with aggressive history clamping to avoid ghosting/details loss.'
    });
    var taaHistoryWeightCtrl = addControl(renderFolder, parameters.taa, 'history_weight', {
        min: 0.0,
        max: 0.98,
        step: 0.01,
        name: 'taa history weight',
        trackQualityPreset: true,
        onChange: onTaaSettingChange,
        help: 'Higher = steadier image but more trailing risk. Lower = cleaner motion.'
    });
    var taaClipBoxCtrl = addControl(renderFolder, parameters.taa, 'clip_box', {
        min: 0.01,
        max: 0.5,
        step: 0.01,
        name: 'taa clip box',
        trackQualityPreset: true,
        onChange: onTaaSettingChange,
        help: 'History clamp size. Lower values reject stale history more aggressively.'
    });
    var taaMotionRejectCtrl = addControl(renderFolder, parameters.taa, 'motion_rejection', {
        min: 0.0,
        max: 20.0,
        step: 0.1,
        name: 'taa motion reject',
        trackQualityPreset: true,
        onChange: onTaaSettingChange,
        help: 'How quickly history influence drops as camera motion increases.'
    });
    var taaMaxDeltaCtrl = addControl(renderFolder, parameters.taa, 'max_camera_delta', {
        min: 0.005,
        max: 0.5,
        step: 0.005,
        name: 'taa max cam delta',
        trackQualityPreset: true,
        onChange: onTaaSettingChange,
        help: 'Hard camera-motion cutoff where TAA history is fully discarded.'
    });
    var taaMotionClipCtrl = addControl(renderFolder, parameters.taa, 'motion_clip_scale', {
        min: 0.0,
        max: 2.0,
        step: 0.01,
        name: 'taa motion clip',
        trackQualityPreset: true,
        onChange: onTaaSettingChange,
        help: 'Extra clamp expansion with motion. Higher values reduce ghosting further in movement.'
    });
    addControl(renderFolder, parameters, 'rk4_integration', {
        name: 'RK4 integration',
        trackQualityPreset: true,
        onChange: updateShader,
        help: 'Higher-order integration for better stability in curved trajectories.'
    });
    renderFolder.open();

    return {
        presetObj: presetObj,
        blackHolePresetController: blackHolePresetController,
        qualityPresetController: qualityPresetController,
        taaEnabledCtrl: taaEnabledCtrl,
        taaHistoryWeightCtrl: taaHistoryWeightCtrl,
        taaClipBoxCtrl: taaClipBoxCtrl,
        taaMotionRejectCtrl: taaMotionRejectCtrl,
        taaMaxDeltaCtrl: taaMaxDeltaCtrl,
        taaMotionClipCtrl: taaMotionClipCtrl
    };
}
