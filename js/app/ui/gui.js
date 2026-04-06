// Role: dat.GUI configuration â€” builds every control panel and folder, wires up
//       all parameter bindings, manages conditional visibility for jet / torus /
//       planet controls, and defines applyBlackHolePreset() to bulk-apply preset
//       data from presets.js and quality-presets.js. Must be called after init()
//       has created the scene.

import { $, dat } from '../vendor.js';
import {
    DISK_TEMPERATURE_MIN,
    DISK_TEMPERATURE_MAX,
    shader,
    scene,
    camera,
    cameraControls,
    cameraPan,
    observer,
    updateUniforms,
    setDistanceController,
    setRefreshAllControllersGlobal
} from '../core/runtime/runtime-state.js';
import {
    OBSERVER_ORBIT_MIN,
    OBSERVER_DISTANCE_MIN,
    OBSERVER_DISTANCE_MAX,
    PLANET_ORBIT_MIN,
    clampObserverDistance
} from '../core/observer.js';
import { updateCamera } from '../scene/camera.js';
import {
    applyQualityPresetValues,
    QUALITY_PRESETS,
    QUALITY_PRESET_LABELS,
    KERR_MODE_LABELS
} from './quality-presets.js';
import { BH_PRESETS } from './presets.js';
import { diveState, hoverState } from '../core/scenarios/scenario-state.js';
import { startDive, resetDive, seekDive } from '../core/scenarios/dive.js';
import { startHover, resetHover, seekHover } from '../core/scenarios/hover.js';
import {
    toggleAnimationTimelineCapture,
    setAnimationTimelineCaptureCameraSmoothingEnabled,
    updateAnimationTimelineCaptureUi
} from '../core/scenarios/animation-capture.js';
import { buildTimelinePanel } from '../presentation/editor/timeline-panel.js';
import { getTimelinePanelBinding } from '../core/runtime/ui-bindings.js';
import { registerBlackHoleUiBinding } from '../core/runtime/runtime-registry.js';
import { setupAnimationsPanel } from './panels/animations-panel.js';
import { setupObserverOrbitWidget } from './widgets/observer-orbit-widget.js';
import { setupControlsPanel } from './panels/controls-panel.js';
import { setupRenderingFolder } from './panels/rendering-folder.js';
import { setupAccretionLookFolders } from './panels/accretion-look-folders.js';
import { setupSimulationFolders } from './panels/simulation-folders.js';

export function setupGUI() {

    var hint = $('#hint-text');
    var p = shader.parameters;

    function updateShader() {
        hint.hide();
        scene.updateShader();
    }

    function updateUniformsLive() {
        updateUniforms();
        shader.needsUpdate = true;
    }

    function observerControlHelp(prefix) {
        var lead = prefix ? (prefix + ' ') : '';
        return lead + 'Left drag: orbit, Right drag: pan, Left+Right drag: roll.';
    }

    function showObserverControlHint(prefix) {
        hint.text(observerControlHelp(prefix));
        hint.stop(true, true).fadeIn(120).delay(1400).fadeOut(360);
    }

    function resetObserverCamera() {
        if (cameraControls) {
            cameraControls.reset();
            cameraControls.target.set(0, 0, 0);
        }
        if (cameraPan) {
            cameraPan.set(0, 0);
        }
        updateCamera();
        shader.needsUpdate = true;
        hint.text('Camera reset. ' + observerControlHelp(''));
        hint.stop(true, true).fadeIn(120).delay(900).fadeOut(360);
    }

    var gui = new dat.GUI({ width: 360 });

    setupControlsPanel(gui);

    var syncObserverWidgetControls = null;

    // Recursively update all dat.GUI controllers to reflect programmatic changes.
    // dat.GUI does NOT auto-sync the displayed value when the bound property
    // is changed from code â€” you must call updateDisplay() on each controller.
    function refreshAllControllers() {
        function recurse(folder) {
            for (var i = 0; i < folder.__controllers.length; i++) {
                folder.__controllers[i].updateDisplay();
            }
            for (var key in folder.__folders) {
                recurse(folder.__folders[key]);
            }
        }
        recurse(gui);
        if (syncObserverWidgetControls) syncObserverWidgetControls();
    }
    setRefreshAllControllersGlobal(refreshAllControllers);
    registerBlackHoleUiBinding('refreshControllers', refreshAllControllers);

    function setGuiRowClass(guiEl, klass) {
        $(guiEl.domElement).parent().parent().addClass(klass);
    }

    function addHelpText(controller, text) {
        if (!text) return;
        // Use the parent .cr row for the tooltip
        var row = $(controller.domElement).closest('.cr')[0];
        if (row) row.title = text;
    }

    var blackHolePresetController = null;
    var qualityPresetController = null;
    var blackHoleCustomState = null;
    var qualityCustomState = null;

    function deepClonePlain(value) {
        return JSON.parse(JSON.stringify(value));
    }

    function captureQualityPresetState() {
        return {
            n_steps: p.n_steps,
            sample_count: p.sample_count,
            max_revolutions: p.max_revolutions,
            rk4_integration: p.rk4_integration,
            cinematic_tonemap: p.cinematic_tonemap,
            resolution_scale: p.resolution_scale,
            taa_enabled: p.taa_enabled,
            taa: deepClonePlain(p.taa)
        };
    }

    function restoreQualityPresetState(state) {
        if (!state) return;
        p.n_steps = state.n_steps;
        p.sample_count = state.sample_count;
        p.max_revolutions = state.max_revolutions;
        p.rk4_integration = state.rk4_integration;
        p.cinematic_tonemap = state.cinematic_tonemap;
        p.resolution_scale = state.resolution_scale;
        p.taa_enabled = state.taa_enabled;
        p.taa.history_weight = state.taa.history_weight;
        p.taa.clip_box = state.taa.clip_box;
        p.taa.motion_rejection = state.taa.motion_rejection;
        p.taa.max_camera_delta = state.taa.max_camera_delta;
        p.taa.motion_clip_scale = state.taa.motion_clip_scale;
    }

    function captureBlackHolePresetState() {
        return {
            black_hole: deepClonePlain(p.black_hole),
            accretion_disk: p.accretion_disk,
            accretion_mode: p.accretion_mode,
            disk_self_irradiation: p.disk_self_irradiation,
            disk_temperature: p.disk_temperature,
            torus: deepClonePlain(p.torus),
            slim: deepClonePlain(p.slim),
            jet: deepClonePlain(p.jet),
            grmhd: deepClonePlain(p.grmhd),
            beaming: p.beaming,
            physical_beaming: p.physical_beaming,
            doppler_shift: p.doppler_shift,
            look: {
                disk_gain: p.look.disk_gain,
                glow: p.look.glow,
                tonemap_mode: p.look.tonemap_mode
            },
            bloom: deepClonePlain(p.bloom)
        };
    }

    function restoreBlackHolePresetState(state) {
        if (!state) return;
        p.black_hole.spin_enabled = state.black_hole.spin_enabled;
        p.black_hole.spin = state.black_hole.spin;
        p.black_hole.spin_strength = state.black_hole.spin_strength;
        p.accretion_disk = state.accretion_disk;
        p.accretion_mode = state.accretion_mode;
        if (state.disk_self_irradiation !== undefined) p.disk_self_irradiation = state.disk_self_irradiation;
        p.disk_temperature = state.disk_temperature;
        p.torus.r0 = state.torus.r0;
        p.torus.h_ratio = state.torus.h_ratio;
        p.torus.radial_falloff = state.torus.radial_falloff;
        p.torus.opacity = state.torus.opacity;
        p.torus.outer_radius = state.torus.outer_radius;
        p.slim.h_ratio = state.slim.h_ratio;
        p.slim.opacity = state.slim.opacity;
        p.slim.puff_factor = state.slim.puff_factor;
        p.jet.enabled = state.jet.enabled;
        p.jet.mode = state.jet.mode;
        p.jet.half_angle = state.jet.half_angle;
        p.jet.lorentz_factor = state.jet.lorentz_factor;
        p.jet.brightness = state.jet.brightness;
        p.jet.length = state.jet.length;
        p.jet.magnetization = state.jet.magnetization;
        p.jet.knot_spacing = state.jet.knot_spacing;
        p.jet.corona_brightness = state.jet.corona_brightness;
        p.jet.base_width = state.jet.base_width;
        p.jet.corona_extent = state.jet.corona_extent;
        if (state.grmhd) {
            p.grmhd.enabled = state.grmhd.enabled;
            p.grmhd.r_high = state.grmhd.r_high;
            p.grmhd.magnetic_beta = state.grmhd.magnetic_beta;
            p.grmhd.mad_flux = state.grmhd.mad_flux;
            p.grmhd.density_scale = state.grmhd.density_scale;
            p.grmhd.turbulence_amp = state.grmhd.turbulence_amp;
            p.grmhd.electron_kappa = state.grmhd.electron_kappa;
            p.grmhd.magnetic_field_str = state.grmhd.magnetic_field_str;
        }
        p.beaming = state.beaming;
        p.physical_beaming = state.physical_beaming;
        p.doppler_shift = state.doppler_shift;
        p.look.disk_gain = state.look.disk_gain;
        p.look.glow = state.look.glow;
        p.look.tonemap_mode = state.look.tonemap_mode;
        p.bloom.enabled = state.bloom.enabled;
        p.bloom.strength = state.bloom.strength;
        p.bloom.threshold = state.bloom.threshold;
        p.bloom.radius = state.bloom.radius;
    }

    function markBlackHolePresetCustom() {
        blackHoleCustomState = captureBlackHolePresetState();
        if (typeof presetObj !== 'undefined' && presetObj.preset !== 'Custom') {
            presetObj.preset = 'Custom';
            if (blackHolePresetController) blackHolePresetController.updateDisplay();
        }
    }

    function markQualityPresetCustom() {
        qualityCustomState = captureQualityPresetState();
        if (p.quality !== 'custom') {
            p.quality = 'custom';
            if (qualityPresetController) qualityPresetController.updateDisplay();
        }
    }

    function withBlackHolePresetTracking(handler) {
        return function() {
            markBlackHolePresetCustom();
            if (handler) return handler.apply(this, arguments);
        };
    }

    function withQualityPresetTracking(handler) {
        return function() {
            markQualityPresetCustom();
            if (handler) return handler.apply(this, arguments);
        };
    }

    function addControl(folder, obj, key, cfg) {
        cfg = cfg || {};
        var c;
        if (cfg.options) c = folder.add(obj, key, cfg.options);
        else c = folder.add(obj, key);
        // dat.GUI may replace numeric controllers when min/max are set.
        // Keep the returned controller so later name/onChange/visibility
        // always target the live row in the DOM.
        if (cfg.min !== undefined) c = c.min(cfg.min) || c;
        if (cfg.max !== undefined) c = c.max(cfg.max) || c;
        if (cfg.step !== undefined) c = c.step(cfg.step) || c;
        if (cfg.name) c = c.name(cfg.name) || c;
        var onChange = cfg.onChange;
        if (cfg.trackBlackHolePreset) onChange = withBlackHolePresetTracking(onChange);
        if (cfg.trackQualityPreset) onChange = withQualityPresetTracking(onChange);
        if (onChange) c = c.onChange(onChange) || c;
        if (cfg.help) addHelpText(c, cfg.help);
        if (cfg.className) setGuiRowClass(c, cfg.className);
        return c;
    }

    function setControlVisible(controller, isVisible) {
        if (!controller) return;
        // dat.GUI stores the row <li> on __li; this is the most reliable target.
        var row = controller && controller.__li ? controller.__li :
            $(controller.domElement).closest('li')[0];
        if (!row) return;
        row.style.display = isVisible ? '' : 'none';
    }

    function setControlsVisible(controllers, isVisible) {
        controllers.forEach(function(ctrl) {
            setControlVisible(ctrl, isVisible);
        });
    }

    function getOptionLabel(options, value) {
        for (var label in options) {
            if (options.hasOwnProperty(label) && options[label] === value) {
                return label;
            }
        }
        return String(value);
    }

    var updateDependentVisibility = function() {};

    function applyKerrMode(mode) {
        // Re-apply quality preset with new mode (preset values depend on kerr_mode)
        applyQualityPresetInternal(p.quality);
        hint.text('Solver mode: ' + getOptionLabel(KERR_MODE_LABELS, mode));
        hint.stop(true, true).fadeIn(120).delay(900).fadeOut(350);
    }

    function applyQualityPresetInternal(value) {
        $('.planet-controls').show();

        if (value === 'custom') {
            restoreQualityPresetState(qualityCustomState);
        } else {
            var preset = null;
            if (typeof applyQualityPresetValues === 'function') {
                preset = applyQualityPresetValues(p, value);
            } else {
                preset = QUALITY_PRESETS[value];
                if (preset) {
                    var isKerr = (
                        p.kerr_mode === 'kerr_inspired_disk_velocity' ||
                        p.kerr_mode === 'realtime_full_kerr_core'
                    );
                    var modeValues = isKerr ? preset.kerr : preset.standard;
                    if (modeValues) {
                        p.quality = value;
                        p.n_steps = modeValues.n_steps;
                        p.sample_count = modeValues.sample_count;
                        p.max_revolutions = modeValues.max_revolutions;
                        p.rk4_integration = modeValues.rk4_integration;
                        p.cinematic_tonemap = preset.cinematic_tonemap;
                        p.resolution_scale = preset.resolution_scale;
                        p.taa_enabled = preset.taa_enabled;
                        p.taa.history_weight = preset.taa.history_weight;
                        p.taa.clip_box = preset.taa.clip_box;
                        p.taa.motion_rejection = preset.taa.motion_rejection;
                        p.taa.max_camera_delta = preset.taa.max_camera_delta;
                        p.taa.motion_clip_scale = preset.taa.motion_clip_scale;
                    } else {
                        preset = null;
                    }
                }
            }
            if (!preset) return;

            if (preset.hide_planet_controls) {
                $('.planet-controls').hide();
            }
        }

        if (typeof applyRenderScaleFromSettings === 'function') {
            applyRenderScaleFromSettings();
        }
        if (typeof resetTemporalAAHistory === 'function') {
            resetTemporalAAHistory();
        }
        updateDependentVisibility();
        updateShader();
        refreshAllControllers();
    }

    function applyQualityPreset(value) {
        applyQualityPresetInternal(value);
        if (value === 'custom') {
            qualityCustomState = captureQualityPresetState();
        }
    }

    // applyBlackHolePreset is assigned below, after visibility helpers are defined.
    var presetObj = { preset: 'Default' };
    var applyBlackHolePreset;
    var renderingFolder = setupRenderingFolder({
        gui: gui,
        parameters: p,
        addControl: addControl,
        qualityPresetLabels: QUALITY_PRESET_LABELS,
        kerrModeLabels: KERR_MODE_LABELS,
        applyBlackHolePreset: function(name) {
            if (applyBlackHolePreset) applyBlackHolePreset(name);
        },
        applyQualityPreset: applyQualityPreset,
        applyKerrMode: applyKerrMode,
        updateShader: updateShader,
        onResolutionScaleChange: function() {
            if (typeof applyRenderScaleFromSettings === 'function') {
                applyRenderScaleFromSettings();
            }
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
        },
        onTaaEnabledChange: function() {
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
            updateDependentVisibility();
            shader.needsUpdate = true;
        },
        onTaaSettingChange: function() {
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
            shader.needsUpdate = true;
        }
    });
    presetObj = renderingFolder.presetObj;
    blackHolePresetController = renderingFolder.blackHolePresetController;
    qualityPresetController = renderingFolder.qualityPresetController;
    var taaEnabledCtrl = renderingFolder.taaEnabledCtrl;
    var taaHistoryWeightCtrl = renderingFolder.taaHistoryWeightCtrl;
    var taaClipBoxCtrl = renderingFolder.taaClipBoxCtrl;
    var taaMotionRejectCtrl = renderingFolder.taaMotionRejectCtrl;
    var taaMaxDeltaCtrl = renderingFolder.taaMaxDeltaCtrl;
    var taaMotionClipCtrl = renderingFolder.taaMotionClipCtrl;

    applyQualityPreset(p.quality);
    blackHoleCustomState = captureBlackHolePresetState();
    qualityCustomState = captureQualityPresetState();

    // Custom scroll handler to control observer distance
    // Placed here so we have access to refreshAllControllers and distanceController
    renderer.domElement.addEventListener( 'wheel', function(e) {
        e.preventDefault();
        var delta = e.deltaY > 0 ? 1.15 : 0.87; // zoom out / zoom in
        var newDist = p.observer.distance * delta;
        newDist = clampObserverDistance(newDist, p.observer.motion);
        p.observer.distance = newDist;
        updateCamera();
        shader.needsUpdate = true;
        refreshAllControllers();
    }, { passive: false } );

    var accretionLookFolders = setupAccretionLookFolders({
        gui: gui,
        parameters: p,
        observer: observer,
        addControl: addControl,
        setControlVisible: setControlVisible,
        setControlsVisible: setControlsVisible,
        getUpdateDependentVisibility: function() { return updateDependentVisibility; },
        updateShader: updateShader,
        updateUniformsLive: updateUniformsLive,
        markShaderDirty: function() { shader.needsUpdate = true; },
        diskTemperatureMin: DISK_TEMPERATURE_MIN,
        diskTemperatureMax: DISK_TEMPERATURE_MAX
    });
    var accretionModeCtrl = accretionLookFolders.accretionModeCtrl;
    var diskSelfIrradiationCtrl = accretionLookFolders.diskSelfIrradiationCtrl;
    var diskTempCtrl = accretionLookFolders.diskTempCtrl;
    var torusRows = accretionLookFolders.torusRows;
    var slimRows = accretionLookFolders.slimRows;
    var grmhdEnabledCtrl = accretionLookFolders.grmhdEnabledCtrl;
    var grmhdRows = accretionLookFolders.grmhdRows;
    var jetCommonCtrls = accretionLookFolders.jetCommonCtrls;
    var jetPhysicalCtrls = accretionLookFolders.jetPhysicalCtrls;
    var spinCtrl = accretionLookFolders.spinCtrl;
    var spinStrengthCtrl = accretionLookFolders.spinStrengthCtrl;
    var lookDiskGainCtrl = accretionLookFolders.lookDiskGainCtrl;
    var lookGlowCtrl = accretionLookFolders.lookGlowCtrl;
    var lookDopplerBoostCtrl = accretionLookFolders.lookDopplerBoostCtrl;
    var lookAberrationStrengthCtrl = accretionLookFolders.lookAberrationStrengthCtrl;
    var bloomStrengthCtrl = accretionLookFolders.bloomStrengthCtrl;
    var bloomThresholdCtrl = accretionLookFolders.bloomThresholdCtrl;
    var bloomRadiusCtrl = accretionLookFolders.bloomRadiusCtrl;

    applyBlackHolePreset = function(name) {
        if (name === 'Custom') {
            restoreBlackHolePresetState(blackHoleCustomState);
            updateDependentVisibility();
            updateCamera();
            updateShader();
            refreshAllControllers();
            blackHoleCustomState = captureBlackHolePresetState();
            return;
        }
        var preset = BH_PRESETS[name];
        if (!preset) return;

        p.black_hole.spin_enabled    = preset.spin_enabled;
        p.black_hole.spin            = preset.spin;
        p.black_hole.spin_strength   = preset.spin_strength;
        p.accretion_disk             = preset.accretion_disk;
        p.accretion_mode             = preset.accretion_mode;
        p.disk_temperature           = preset.disk_temperature;
        p.torus.r0                   = preset.torus.r0;
        p.torus.h_ratio              = preset.torus.h_ratio;
        p.torus.radial_falloff       = preset.torus.radial_falloff;
        p.torus.opacity              = preset.torus.opacity;
        p.torus.outer_radius         = preset.torus.outer_radius;
        p.slim.h_ratio               = preset.slim.h_ratio;
        p.slim.opacity               = preset.slim.opacity;
        p.slim.puff_factor           = preset.slim.puff_factor;
        p.jet.enabled                = preset.jet.enabled;
        p.jet.mode                   = preset.jet.mode;
        p.jet.half_angle             = preset.jet.half_angle;
        p.jet.lorentz_factor         = preset.jet.lorentz_factor;
        p.jet.brightness             = preset.jet.brightness;
        p.jet.length                 = preset.jet.length;
        p.jet.magnetization          = preset.jet.magnetization;
        p.jet.knot_spacing           = preset.jet.knot_spacing;
        p.jet.corona_brightness      = preset.jet.corona_brightness;
        p.jet.base_width             = preset.jet.base_width;
        p.jet.corona_extent          = preset.jet.corona_extent;
        if (preset.grmhd) {
            p.grmhd.enabled              = preset.grmhd.enabled;
            p.grmhd.r_high               = preset.grmhd.r_high;
            p.grmhd.magnetic_beta        = preset.grmhd.magnetic_beta;
            p.grmhd.mad_flux             = preset.grmhd.mad_flux;
            p.grmhd.density_scale        = preset.grmhd.density_scale;
            p.grmhd.turbulence_amp       = preset.grmhd.turbulence_amp;
            p.grmhd.electron_kappa       = preset.grmhd.electron_kappa;
            p.grmhd.magnetic_field_str   = preset.grmhd.magnetic_field_str;
        } else {
            // Default GRMHD off when preset doesn't specify it
            p.grmhd.enabled = false;
        }
        p.beaming                    = preset.beaming;
        p.physical_beaming           = preset.physical_beaming;
        p.doppler_shift              = preset.doppler_shift;
        p.look.disk_gain             = preset.disk_gain;
        p.look.glow                  = preset.glow;
        p.look.tonemap_mode          = preset.tonemap_mode;

        if (preset.bloom) {
            p.bloom.enabled   = preset.bloom.enabled;
            p.bloom.strength  = preset.bloom.strength;
            p.bloom.threshold = preset.bloom.threshold;
            p.bloom.radius    = preset.bloom.radius;
        }

        updateDependentVisibility();
        updateCamera();
        updateShader();
        refreshAllControllers();
    };

    var simulationFolders = setupSimulationFolders({
        gui: gui,
        parameters: p,
        observer: observer,
        addControl: addControl,
        updateUniforms: updateUniforms,
        updateShader: updateShader,
        markShaderDirty: function() { shader.needsUpdate = true; },
        planetOrbitMin: PLANET_ORBIT_MIN,
        getUpdateDependentVisibility: function() { return updateDependentVisibility; }
    });
    var planetDistanceCtrl = simulationFolders.planetDistanceCtrl;
    var planetRadiusCtrl = simulationFolders.planetRadiusCtrl;
    var physicalBeamingCtrl = simulationFolders.physicalBeamingCtrl;
    var lorentzContractionCtrl = simulationFolders.lorentzContractionCtrl;


    updateDependentVisibility = function() {
        var diskEnabled = !!p.accretion_disk;
        var thinDiskEnabled = diskEnabled && p.accretion_mode === 'thin_disk';
        var thickTorusEnabled = diskEnabled && p.accretion_mode === 'thick_torus';
        var slimDiskEnabled = diskEnabled && p.accretion_mode === 'slim_disk';

        setControlVisible(accretionModeCtrl, diskEnabled);
        setControlVisible(diskSelfIrradiationCtrl, diskEnabled);
        setControlVisible(diskTempCtrl, diskEnabled);
        setControlsVisible(torusRows, thickTorusEnabled);
        setControlsVisible(slimRows, slimDiskEnabled);

        // GRMHD controls: toggle visible when disk is enabled, parameters visible when GRMHD mode is on
        var grmhdEnabled = diskEnabled && !!p.grmhd.enabled;
        setControlVisible(grmhdEnabledCtrl, diskEnabled);
        setControlsVisible(grmhdRows, grmhdEnabled);

        var jetEnabled = !!p.jet.enabled;
        setControlsVisible(jetCommonCtrls, jetEnabled);
        setControlsVisible(jetPhysicalCtrls, jetEnabled && p.jet.mode === 'physical');

        var spinEnabled = !!p.black_hole.spin_enabled;
        setControlVisible(spinCtrl, spinEnabled);
        setControlVisible(spinStrengthCtrl, spinEnabled);

        var bloomEnabled = !!p.bloom.enabled;
        setControlsVisible([bloomStrengthCtrl, bloomThresholdCtrl, bloomRadiusCtrl], bloomEnabled);

        var beamingEnabled = !!p.beaming;
        var physicalBeamingEnabled = beamingEnabled && !!p.physical_beaming;
        setControlVisible(physicalBeamingCtrl, beamingEnabled);
        setControlVisible(lookDopplerBoostCtrl, beamingEnabled && !physicalBeamingEnabled);

        setControlVisible(lookAberrationStrengthCtrl, !!p.aberration);
        setControlVisible(lookGlowCtrl, thinDiskEnabled);
        setControlVisible(lookDiskGainCtrl, diskEnabled || jetEnabled);

        var planetEnabled = !!p.planet.enabled;
        setControlsVisible([planetDistanceCtrl, planetRadiusCtrl], planetEnabled);
        setControlVisible(lorentzContractionCtrl, planetEnabled);

        var taaEnabled = !!p.taa_enabled;
        setControlVisible(taaEnabledCtrl, true);
        setControlsVisible(
            [taaHistoryWeightCtrl, taaClipBoxCtrl, taaMotionRejectCtrl, taaMaxDeltaCtrl, taaMotionClipCtrl],
            taaEnabled
        );
    };
    updateDependentVisibility();

    // Side panels extracted into focused helpers to keep the main GUI wiring readable.
    setupAnimationsPanel({
        initialObserverDistance: p.observer.distance,
        diveState: diveState,
        hoverState: hoverState,
        startDive: startDive,
        resetDive: resetDive,
        seekDive: seekDive,
        startHover: startHover,
        resetHover: resetHover,
        seekHover: seekHover,
        toggleAnimationTimelineCapture: toggleAnimationTimelineCapture,
        setAnimationTimelineCaptureCameraSmoothingEnabled:
            setAnimationTimelineCaptureCameraSmoothingEnabled,
        updateAnimationTimelineCaptureUi: updateAnimationTimelineCaptureUi,
        buildTimelinePanel: buildTimelinePanel,
        getTimelinePanelBinding: getTimelinePanelBinding
    });

    var observerOrbitWidget = setupObserverOrbitWidget({
        parameters: p,
        observerOrbitMin: OBSERVER_ORBIT_MIN,
        observerDistanceMin: OBSERVER_DISTANCE_MIN,
        observerDistanceMax: OBSERVER_DISTANCE_MAX,
        clampObserverDistance: clampObserverDistance,
        updateCamera: updateCamera,
        markShaderDirty: function() {
            shader.needsUpdate = true;
        },
        updateShader: updateShader,
        resetObserverCamera: resetObserverCamera,
        showObserverControlHint: showObserverControlHint,
        registerBlackHoleUiBinding: registerBlackHoleUiBinding,
        setDistanceController: setDistanceController
    });
    syncObserverWidgetControls = observerOrbitWidget.syncControls;
}


