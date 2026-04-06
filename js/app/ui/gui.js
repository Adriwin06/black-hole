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

    // â”€â”€ Controls Side Panel â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    (function setupControlsPanel() {
        var STORED_WIDTH_KEY = 'black-hole.controls-panel.width';
        var DEFAULT_WIDTH = 370;
        var MIN_WIDTH = 300;
        var MAX_WIDTH = 600;

        var ctrlPanel = document.createElement('div');
        ctrlPanel.id = 'controls-panel';
        ctrlPanel.className = 'sp-panel sp-panel--right sp-panel--collapsed';
        ctrlPanel.innerHTML =
            '<div class="sp-resize sp-resize--left"></div>' +
            '<div class="sp-header">' +
                '<span class="sp-title">CONTROLS</span>' +
                '<button id="controls-close-btn" class="sp-close" type="button">&times;</button>' +
            '</div>' +
            '<div id="controls-content" class="sp-content"></div>';
        document.body.appendChild(ctrlPanel);

        // Move dat.GUI into the panel
        var dgAc = gui.domElement.parentElement;
        document.getElementById('controls-content').appendChild(dgAc);

        // Create open button
        var openBtn = document.createElement('button');
        openBtn.id = 'controls-open-btn';
        openBtn.className = 'sp-open-btn sp-open-btn--right';
        openBtn.type = 'button';
        openBtn.innerHTML = 'CONTROLS &#9654;';
        openBtn.title = 'Open Controls panel';
        document.body.appendChild(openBtn);

        // Panel state
        var panelOpen = false;
        function setPanelOpen(open) {
            panelOpen = open;
            ctrlPanel.classList.toggle('sp-panel--collapsed', !open);
            openBtn.classList.toggle('sp-open-btn--hidden', open);
            document.body.classList.toggle('has-controls-panel', open);
        }

        openBtn.addEventListener('click', function() { setPanelOpen(true); });
        document.getElementById('controls-close-btn')
            .addEventListener('click', function() { setPanelOpen(false); });

        // Resize (left edge â€” drag left to widen)
        var resizer = ctrlPanel.querySelector('.sp-resize');
        var resizing = false;
        var startX = 0;
        var startW = DEFAULT_WIDTH;

        function readWidth() {
            try { var w = parseFloat(localStorage.getItem(STORED_WIDTH_KEY)); return isFinite(w) ? w : null; } catch(e) { return null; }
        }
        function saveWidth(w) {
            try { localStorage.setItem(STORED_WIDTH_KEY, String(Math.round(w))); } catch(e) {}
        }
        function applyWidth(w, persist) {
            w = Math.round(Math.max(MIN_WIDTH, Math.min(MAX_WIDTH, w)));
            ctrlPanel.style.width = w + 'px';
            if (persist) saveWidth(w);
        }

        applyWidth(readWidth() || DEFAULT_WIDTH, false);

        resizer.addEventListener('pointerdown', function(e) {
            if (e.button !== 0) return;
            resizing = true;
            startX = e.clientX;
            startW = ctrlPanel.getBoundingClientRect().width || DEFAULT_WIDTH;
            document.body.classList.add('sp-resizing');
            if (resizer.setPointerCapture) resizer.setPointerCapture(e.pointerId);
            document.addEventListener('pointermove', onMove);
            document.addEventListener('pointerup', onEnd);
            e.preventDefault();
        });
        function onMove(e) {
            if (!resizing) return;
            applyWidth(startW + (startX - e.clientX), false);
            e.preventDefault();
        }
        function onEnd() {
            if (!resizing) return;
            resizing = false;
            document.body.classList.remove('sp-resizing');
            applyWidth(parseFloat(ctrlPanel.style.width) || DEFAULT_WIDTH, true);
            document.removeEventListener('pointermove', onMove);
            document.removeEventListener('pointerup', onEnd);
        }

        // Auto-open on desktop
        if (!(window.matchMedia && window.matchMedia('(max-width: 960px)').matches)) {
            setPanelOpen(true);
        }
    })();

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

    // â”€â”€â”€ Real Black Hole Presets â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    // applyBlackHolePreset is assigned below, after visibility helpers are defined.
    var presetObj = { preset: 'Default' };
    var applyBlackHolePreset;
    var presetsFolder = gui.addFolder('Black Hole Presets');
    blackHolePresetController = addControl(presetsFolder, presetObj, 'preset', {
        name: 'preset',
        options: ['Custom', 'Default', 'M87*', 'Sgr A*', 'Cygnus X-1', 'GRS 1915+105', 'Gargantua (Interstellar visuals)', 'Schwarzschild'],
        onChange: function(val) { if (applyBlackHolePreset) applyBlackHolePreset(val); },
        help: 'Load a literature-inspired black-hole preset. These are illustrative starting points, not definitive measurements.'
    });
    presetsFolder.open();
    // â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

    var renderFolder = gui.addFolder('Rendering');
    qualityPresetController = addControl(renderFolder, p, 'quality', {
        options: QUALITY_PRESET_LABELS,
        name: 'quality preset',
        onChange: applyQualityPreset,
        help: 'Global render preset. Higher presets use more ray-marching steps and samples.'
    });
    addControl(renderFolder, p, 'kerr_mode', {
        options: KERR_MODE_LABELS,
        name: 'solver mode',
        onChange: applyKerrMode,
        help: 'Fast = Schwarzschild/Binet photon lensing with a perturbative frame-drag term. Kerr-inspired disk velocities = same photon lensing, but emitting matter uses Kerr angular velocity in a simplified local-speed model.'
    });

    addControl(renderFolder, p, 'n_steps', {
        min: 20,
        max: 1400,
        step: 1,
        name: 'ray steps',
        trackQualityPreset: true,
        onChange: updateShader,
        help: 'More steps improve thin features and strong lensing, but reduce FPS.'
    });
    addControl(renderFolder, p, 'sample_count', {
        min: 1,
        max: 12,
        step: 1,
        name: 'samples / pixel',
        trackQualityPreset: true,
        onChange: updateShader,
        help: 'Supersampling for anti-aliasing and smoother edges. Higher values are slower.'
    });
    addControl(renderFolder, p, 'max_revolutions', {
        min: 1.0,
        max: 8.0,
        step: 0.1,
        name: 'max orbit turns',
        trackQualityPreset: true,
        onChange: updateShader,
        help: 'How many wrapped photon turns are traced before escape/capture cutoff.'
    });
    addControl(renderFolder, p, 'resolution_scale', {
        min: 0.35,
        max: 2.0,
        step: 0.05,
        name: 'resolution scale',
        trackQualityPreset: true,
        onChange: function() {
            if (typeof applyRenderScaleFromSettings === 'function') {
                applyRenderScaleFromSettings();
            }
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
        },
        help: 'Internal rendering resolution multiplier. <1 lowers cost, >1 supersamples for sharper output.'
    });
    var taaEnabledCtrl = addControl(renderFolder, p, 'taa_enabled', {
        name: 'TAA (stable)',
        trackQualityPreset: true,
        onChange: function() {
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
            updateDependentVisibility();
            shader.needsUpdate = true;
        },
        help: 'Temporal anti-aliasing with aggressive history clamping to avoid ghosting/details loss.'
    });
    var taaHistoryWeightCtrl = addControl(renderFolder, p.taa, 'history_weight', {
        min: 0.0,
        max: 0.98,
        step: 0.01,
        name: 'taa history weight',
        trackQualityPreset: true,
        onChange: function() {
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
            shader.needsUpdate = true;
        },
        help: 'Higher = steadier image but more trailing risk. Lower = cleaner motion.'
    });
    var taaClipBoxCtrl = addControl(renderFolder, p.taa, 'clip_box', {
        min: 0.01,
        max: 0.5,
        step: 0.01,
        name: 'taa clip box',
        trackQualityPreset: true,
        onChange: function() {
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
            shader.needsUpdate = true;
        },
        help: 'History clamp size. Lower values reject stale history more aggressively.'
    });
    var taaMotionRejectCtrl = addControl(renderFolder, p.taa, 'motion_rejection', {
        min: 0.0,
        max: 20.0,
        step: 0.1,
        name: 'taa motion reject',
        trackQualityPreset: true,
        onChange: function() {
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
            shader.needsUpdate = true;
        },
        help: 'How quickly history influence drops as camera motion increases.'
    });
    var taaMaxDeltaCtrl = addControl(renderFolder, p.taa, 'max_camera_delta', {
        min: 0.005,
        max: 0.5,
        step: 0.005,
        name: 'taa max cam delta',
        trackQualityPreset: true,
        onChange: function() {
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
            shader.needsUpdate = true;
        },
        help: 'Hard camera-motion cutoff where TAA history is fully discarded.'
    });
    var taaMotionClipCtrl = addControl(renderFolder, p.taa, 'motion_clip_scale', {
        min: 0.0,
        max: 2.0,
        step: 0.01,
        name: 'taa motion clip',
        trackQualityPreset: true,
        onChange: function() {
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
            shader.needsUpdate = true;
        },
        help: 'Extra clamp expansion with motion. Higher values reduce ghosting further in movement.'
    });
    addControl(renderFolder, p, 'rk4_integration', {
        name: 'RK4 integration',
        trackQualityPreset: true,
        onChange: updateShader,
        help: 'Higher-order integration for better stability in curved trajectories.'
    });
    renderFolder.open();

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

    var diskFolder = gui.addFolder('Accretion disk');
    addControl(diskFolder, p, 'accretion_disk', {
        name: 'enabled',
        trackBlackHolePreset: true,
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Toggles thermal emission from the accretion disk.'
    });
    var accretionModeCtrl = addControl(diskFolder, p, 'accretion_mode', {
        name: 'accretion type',
        options: ['thin_disk', 'thick_torus', 'slim_disk'],
        trackBlackHolePreset: true,
        onChange: function() {
            updateDependentVisibility();
            updateShader();
            observer.turbulenceTimeOffset = -observer.time;
            shader.needsUpdate = true;
        },
        help: 'Thin disk: zero-torque Shakura-Sunyaev / Novikov-Thorne-style proxy (quasars/XRBs). Thick torus: ADAF/RIAF (M87*/Sgr A*). Slim disk: super-Eddington.'
    });
    var diskSelfIrradiationCtrl = addControl(diskFolder, p, 'disk_self_irradiation', {
        name: 'self-irradiation',
        trackBlackHolePreset: true,
        onChange: function() {
            updateShader();
        },
        help: 'Heuristic Cunningham-inspired inner-disk brightening. Increases local flux near the ISCO with a spin-scaled falloff; this is not recursive ray-traced returning radiation.'
    });
    var diskTempCtrl = addControl(diskFolder, p, 'disk_temperature', {
        min: DISK_TEMPERATURE_MIN,
        max: DISK_TEMPERATURE_MAX,
        step: 1,
        name: 'temperature (K)',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Rest-frame disk color temperature before relativistic shifts.'
    });
    var torusCenterCtrl = addControl(diskFolder, p.torus, 'r0', {
        min: 1.5,
        max: 10.0,
        step: 0.1,
        name: 'torus center r',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Center radius of the torus in r_s units (thick torus mode only).'
    });
    var torusHRCtrl = addControl(diskFolder, p.torus, 'h_ratio', {
        min: 0.1,
        max: 1.0,
        step: 0.01,
        name: 'torus H/R',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Height-to-radius ratio of the torus cross-section (thick torus mode only).'
    });
    var torusFalloffCtrl = addControl(diskFolder, p.torus, 'radial_falloff', {
        min: 0.5,
        max: 5.0,
        step: 0.1,
        name: 'radial falloff',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Power-law index for radial emissivity decay. Physical ADAF: 2-4. Lower = flatter profile, higher = more concentrated.'
    });
    var torusOpacityCtrl = addControl(diskFolder, p.torus, 'opacity', {
        min: 0.001,
        max: 0.15,
        step: 0.001,
        name: 'opacity',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Absorption coefficient. Low = optically thin (ADAF). Higher = more self-shielding and defined torus surface.'
    });
    var torusOuterCtrl = addControl(diskFolder, p.torus, 'outer_radius', {
        min: 1.5,
        max: 8.0,
        step: 0.1,
        name: 'outer extent',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Outer edge multiplier relative to torus center r0. Controls how far the torus extends.'
    });

    // â”€â”€â”€ Slim disk controls â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    var slimHRCtrl = addControl(diskFolder, p.slim, 'h_ratio', {
        min: 0.05,
        max: 0.5,
        step: 0.01,
        name: 'slim H/R',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Base height-to-radius ratio of the slim disk. Higher = geometrically thicker.'
    });
    var slimOpacityCtrl = addControl(diskFolder, p.slim, 'opacity', {
        min: 0.1,
        max: 2.0,
        step: 0.01,
        name: 'slim opacity',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Absorption coefficient. Higher = more opaque, surface-like. Lower = more translucent, volumetric.'
    });
    var slimPuffCtrl = addControl(diskFolder, p.slim, 'puff_factor', {
        min: 0.0,
        max: 6.0,
        step: 0.1,
        name: 'ISCO puff',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'How much the disk puffs up near the ISCO due to radiation pressure. Higher = more pronounced thickening.'
    });

    // Store GUI row elements for conditional visibility
    var torusRows = [torusCenterCtrl, torusHRCtrl, torusFalloffCtrl, torusOpacityCtrl, torusOuterCtrl];
    var slimRows = [slimHRCtrl, slimOpacityCtrl, slimPuffCtrl];

    // â”€â”€â”€ GRMHD-inspired toggle â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    var grmhdEnabledCtrl = addControl(diskFolder, p.grmhd, 'enabled', {
        name: 'GRMHD-inspired',
        trackBlackHolePreset: true,
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Enables GRMHD-inspired morphology controls: two-temperature plasma weighting, synchrotron-inspired emissivity corrections, MRI-inspired turbulence, and magnetic-field effects. These are semi-analytic visual models, not full GRMHD evolution.'
    });
    var grmhdRHighCtrl = addControl(diskFolder, p.grmhd, 'r_high', {
        min: 1.0,
        max: 160.0,
        step: 1.0,
        name: 'R_high',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Renderer-side Ti:Te / electron-heating proxy. Higher R_high suppresses high-beta disk-body emission most strongly in the GRMHD torus branch, while thin/slim disks keep the effect deliberately milder so the RGB disk remains visible.'
    });
    var grmhdBetaCtrl = addControl(diskFolder, p.grmhd, 'magnetic_beta', {
        min: 0.01,
        max: 100.0,
        step: 0.1,
        name: 'plasma Î²',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Ratio of gas-to-magnetic pressure in the disk midplane. Low Î² = strongly magnetized, high Î² = gas dominated. GRMHD simulations: Î² ~ 1-30.'
    });
    var grmhdMADCtrl = addControl(diskFolder, p.grmhd, 'mad_flux', {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        name: 'MAD flux',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Magnetically Arrested Disk saturation. 0 = SANE (standard), 1 = full MAD. MAD state creates stronger B-fields, prominent m=1 spiral arms, and more powerful jets (Tchekhovskoy+ 2011).'
    });
    var grmhdDensityCtrl = addControl(diskFolder, p.grmhd, 'density_scale', {
        min: 0.1,
        max: 5.0,
        step: 0.1,
        name: 'density scale',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Normalization of the GRMHD density profile. Higher values increase overall emissivity and absorption.'
    });
    var grmhdTurbCtrl = addControl(diskFolder, p.grmhd, 'turbulence_amp', {
        min: 0.0,
        max: 3.0,
        step: 0.1,
        name: 'MRI turbulence',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Amplitude of MRI-driven turbulent density fluctuations. 0 = smooth, 1 = typical GRMHD, 3 = strongly turbulent. Controls log-normal density PDF and spiral arm contrast.'
    });
    var grmhdKappaCtrl = addControl(diskFolder, p.grmhd, 'electron_kappa', {
        min: 2.5,
        max: 8.0,
        step: 0.1,
        name: 'Îº (non-thermal)',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Kappa-distribution index for non-thermal electrons from magnetic reconnection. Îº ~ 3.5 = strong non-thermal tail, Îº ~ 8 = nearly thermal. (Ball+ 2018, Werner+ 2018).'
    });
    var grmhdBFieldCtrl = addControl(diskFolder, p.grmhd, 'magnetic_field_str', {
        min: 0.1,
        max: 5.0,
        step: 0.1,
        name: 'B-field strength',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Overall magnetic field strength scaling. Affects synchrotron emissivity (j âˆ BÂ²), self-absorption, and electron temperature.'
    });

    var grmhdRows = [grmhdRHighCtrl, grmhdBetaCtrl, grmhdMADCtrl, grmhdDensityCtrl,
                     grmhdTurbCtrl, grmhdKappaCtrl, grmhdBFieldCtrl];

    // Apply initial accretion visibility immediately, even before the global
    // dependency updater is assigned at the end of setupGUI().
    setControlVisible(accretionModeCtrl, !!p.accretion_disk);
    setControlVisible(diskSelfIrradiationCtrl, !!p.accretion_disk);
    setControlVisible(diskTempCtrl, !!p.accretion_disk);
    setControlsVisible(torusRows, !!p.accretion_disk && p.accretion_mode === 'thick_torus');
    setControlsVisible(slimRows, !!p.accretion_disk && p.accretion_mode === 'slim_disk');
    setControlVisible(grmhdEnabledCtrl, !!p.accretion_disk);
    setControlsVisible(grmhdRows, !!p.accretion_disk && !!p.grmhd.enabled);

    diskFolder.open();

    // â”€â”€â”€ Relativistic Jets â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    var jetFolder = gui.addFolder('Relativistic jets');
    addControl(jetFolder, p.jet, 'enabled', {
        name: 'enabled',
        trackBlackHolePreset: true,
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Bipolar analytic jets aligned with the spin axis. The detailed mode adds GRMHD-inspired structure and Blandford-Znajek-like power scaling.'
    });
    var jetModeCtrl = addControl(jetFolder, p.jet, 'mode', {
        name: 'jet model',
        options: ['simple', 'physical'],
        trackBlackHolePreset: true,
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Simple: smooth analytic jet. Physical: more detailed GRMHD-inspired jet model with spine/sheath, reconfinement knots, corona base, and disk occultation; RGB colour still uses an effective-temperature proxy rather than a literal synchrotron spectrum.'
    });
    var jetAngleCtrl = addControl(jetFolder, p.jet, 'half_angle', {
        min: 1.0,
        max: 25.0,
        step: 0.5,
        name: 'half-angle (Â°)',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Opening half-angle of the jet. Typical AGN jets: 2-7Â°. Parabolic collimation narrows the beam with distance.'
    });
    var jetLorentzCtrl = addControl(jetFolder, p.jet, 'lorentz_factor', {
        min: 1.1,
        max: 20.0,
        step: 0.1,
        name: 'Lorentz Î“',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Bulk Lorentz factor. Î“~2-5 for AGN jets, Î“~100+ for GRBs. Controls relativistic beaming.'
    });
    var jetBrightCtrl = addControl(jetFolder, p.jet, 'brightness', {
        min: 0.05,
        max: 3.0,
        step: 0.01,
        name: 'brightness',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Overall jet synchrotron emission strength.'
    });
    var jetLengthCtrl = addControl(jetFolder, p.jet, 'length', {
        min: 5.0,
        max: 60.0,
        step: 1.0,
        name: 'length (r_s)',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Visible jet length in Schwarzschild radii.'
    });
    var jetMagCtrl = addControl(jetFolder, p.jet, 'magnetization', {
        min: 1.0,
        max: 50.0,
        step: 0.5,
        name: 'Ïƒ (magnetization)',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Plasma magnetization at jet base (Ïƒ = BÂ²/4Ï€ÏcÂ²). Higher Ïƒ = more magnetically dominated, narrower spine. MAD jets: Ïƒ ~ 10-30.'
    });
    var jetKnotCtrl = addControl(jetFolder, p.jet, 'knot_spacing', {
        min: 2.0,
        max: 15.0,
        step: 0.5,
        name: 'knot spacing',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Spacing of reconfinement shock knots (in r_s). Observed in M87 (HST-1), 3C 273. Set high to minimize knots.'
    });
    var jetCoronaCtrl = addControl(jetFolder, p.jet, 'corona_brightness', {
        min: 0.0,
        max: 5.0,
        step: 0.1,
        name: 'corona glow',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Brightness of the jet-corona connection at the base. Hot plasma where the funnel meets the inner accretion flow.'
    });
    var jetBaseWidthCtrl = addControl(jetFolder, p.jet, 'base_width', {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        name: 'base width',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Controls how wide the jet funnel is at the base. 0 = narrow collimated base, 1 = wide split-monopole funnel.'
    });
    var jetCoronaExtentCtrl = addControl(jetFolder, p.jet, 'corona_extent', {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        name: 'corona extent',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Radial spread of the corona base emission. 0 = wide spread, 1 = tightly concentrated near the jet axis.'
    });

    // Store jet row groups for conditional visibility
    var jetCommonCtrls = [jetModeCtrl, jetAngleCtrl, jetLorentzCtrl, jetBrightCtrl, jetLengthCtrl];
    var jetPhysicalCtrls = [jetMagCtrl, jetKnotCtrl, jetCoronaCtrl, jetBaseWidthCtrl, jetCoronaExtentCtrl];

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

    var spinFolder = gui.addFolder('Black hole');
    addControl(spinFolder, p.black_hole, 'spin_enabled', {
        name: 'rotation enabled',
        trackBlackHolePreset: true,
        onChange: function() {
            updateDependentVisibility();
            updateUniformsLive();
        },
        help: 'Enables Kerr-like rotation effects. Disable for a Schwarzschild-style shadow.'
    });
    var spinCtrl = addControl(spinFolder, p.black_hole, 'spin', {
        min: -0.99,
        max: 0.99,
        step: 0.01,
        name: 'a/M',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Signed Kerr spin. Positive and negative values mirror the spin direction; the current UI does not expose a separate retrograde-disk toggle.'
    });
    var spinStrengthCtrl = addControl(spinFolder, p.black_hole, 'spin_strength', {
        min: 0.0,
        max: 1.4,
        step: 0.01,
        name: 'shadow squeeze',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Visual strength multiplier for rotation-induced asymmetry in the shadow.'
    });
    spinFolder.open();

    var lookFolder = gui.addFolder('Look');
    addControl(lookFolder, p.look, 'tonemap_mode', {
        options: { 'ACES Filmic': 0, 'AgX': 1, 'Scientific (Log)': 2 },
        name: 'tonemapper',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'ACES: classic cinematic. AgX: better saturation handling for extreme HDR. Scientific: logarithmic false-color like EHT papers.'
    });
    addControl(lookFolder, p.look, 'exposure', {
        min: 0.6,
        max: 2.5,
        step: 0.01,
        name: 'exposure',
        onChange: updateUniformsLive,
        help: 'Overall tone-mapped brightness.'
    });
    var lookDiskGainCtrl = addControl(lookFolder, p.look, 'disk_gain', {
        min: 0.4,
        max: 4.0,
        step: 0.01,
        name: 'disk intensity',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Brightness multiplier applied to the accretion disk emission.'
    });
    var lookGlowCtrl = addControl(lookFolder, p.look, 'glow', {
        min: 0.0,
        max: 2.0,
        step: 0.01,
        name: 'inner glow',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Extra bloom-like emphasis near the hotter inner disk region.'
    });
    var lookDopplerBoostCtrl = addControl(lookFolder, p.look, 'doppler_boost', {
        min: 0.0,
        max: 2.5,
        step: 0.01,
        name: 'doppler boost',
        onChange: updateUniformsLive,
        help: 'Controls visual strength of relativistic beaming contrast.'
    });
    var lookAberrationStrengthCtrl = addControl(lookFolder, p.look, 'aberration_strength', {
        min: 0.0,
        max: 3.0,
        step: 0.01,
        name: 'aberration strength',
        onChange: updateUniformsLive,
        help: 'Scales apparent directional warping from observer motion.'
    });
    addControl(lookFolder, p.look, 'star_gain', {
        min: 0.0,
        max: 2.5,
        step: 0.01,
        name: 'star gain',
        onChange: updateUniformsLive,
        help: 'Brightness multiplier for background stars.'
    });
    addControl(lookFolder, p.look, 'galaxy_gain', {
        min: 0.0,
        max: 2.5,
        step: 0.01,
        name: 'galaxy gain',
        onChange: updateUniformsLive,
        help: 'Brightness multiplier for the background galaxy map.'
    });
    lookFolder.open();

    // â”€â”€â”€ Post-processing (bloom) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    var ppFolder = gui.addFolder('Post-processing');
    addControl(ppFolder, p.bloom, 'enabled', {
        name: 'bloom',
        trackBlackHolePreset: true,
        onChange: function() {
            updateDependentVisibility();
            shader.needsUpdate = true;
        },
        help: 'Screen-space Gaussian bloom approximating lens glow and optical spill. Bright regions of the accretion disk bleed light into surrounding pixels.'
    });
    var bloomStrengthCtrl = addControl(ppFolder, p.bloom, 'strength', {
        min: 0.0,
        max: 2.0,
        step: 0.01,
        name: 'bloom strength',
        trackBlackHolePreset: true,
        onChange: function() { shader.needsUpdate = true; },
        help: 'Overall intensity of the bloom glow added to the image.'
    });
    var bloomThresholdCtrl = addControl(ppFolder, p.bloom, 'threshold', {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        name: 'bloom threshold',
        trackBlackHolePreset: true,
        onChange: function() { shader.needsUpdate = true; },
        help: 'Minimum pixel brightness that contributes to bloom. Lower values bloom more of the disk.'
    });
    var bloomRadiusCtrl = addControl(ppFolder, p.bloom, 'radius', {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        name: 'bloom radius',
        trackBlackHolePreset: true,
        onChange: function() { shader.needsUpdate = true; },
        help: 'Controls the width of the bloom halo. Higher values give wider, more diffuse glow matching real optical PSFs.'
    });
    ppFolder.open();

    var folder = gui.addFolder('Planet');
    addControl(folder, p.planet, 'enabled', {
        name: 'enabled',
        help: 'Adds an orbiting planet used as a reference object.',
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        }
    });
    var planetDistanceCtrl = addControl(folder, p.planet, 'distance', {
        min: PLANET_ORBIT_MIN,
        step: 0.1,
        name: 'distance',
        onChange: updateUniforms,
        help: 'Orbital radius of the planet. Public controls clamp this to r >= 3 r_s so the reference body stays in the stable timelike-orbit regime.'
    });
    var planetRadiusCtrl = addControl(folder, p.planet, 'radius', {
        min: 0.01,
        max: 2.0,
        step: 0.01,
        name: 'radius',
        onChange: updateUniforms,
        help: 'Planet size in simulation units.'
    });
    $(folder.domElement).addClass('planet-controls');
    //folder.open();

    folder = gui.addFolder('Relativistic effects');
    addControl(folder, p, 'aberration', {
        name: 'aberration (ray dir)',
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Changes apparent incoming ray direction due to observer velocity.'
    });
    addControl(folder, p, 'beaming', {
        name: 'beaming (intensity)',
        trackBlackHolePreset: true,
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Applies relativistic intensity boosting/dimming to thermal emitters and the background-sky proxy. Jet branches keep their own synchrotron beaming model.'
    });
    var physicalBeamingCtrl = addControl(folder, p, 'physical_beaming', {
        name: 'physical (DÂ³ Liouville)',
        trackBlackHolePreset: true,
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Uses physically motivated Liouville transfer scaling for thermal emitters and the background-sky proxy instead of the softened cinematic curve.'
    });
    addControl(folder, p, 'doppler_shift', {
        name: 'doppler shift (color)',
        trackBlackHolePreset: true,
        onChange: updateShader,
        help: 'Shifts thermal emitters and the background-sky proxy by red/blue shift factors. Physical jet mode keeps its own synchrotron temperature-proxy shift.'
    });
    addControl(folder, p, 'gravitational_time_dilation', {
        name: 'time dilation',
        onChange: updateShader,
        className: 'planet-controls indirect-planet-controls',
        help: 'Accounts for rate differences between local and distant observer clocks.'
    });
    var lorentzContractionCtrl = addControl(folder, p, 'lorentz_contraction', {
        name: 'lorentz contraction',
        onChange: updateShader,
        className: 'planet-controls indirect-planet-controls',
        help: 'Applies length contraction effects to moving scene elements.'
    });

    folder.open();

    folder = gui.addFolder('Time');
    addControl(folder, p, 'light_travel_time', {
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Enables retarded-time rendering, where events are seen after light delay.'
    });
    addControl(folder, p, 'time_scale', {
        min: 0,
        max: 6,
        step: 0.01,
        name: 'time scale',
        help: 'Simulation clock multiplier.'
    });
    var resetSimObj = { 'reset simulation': function() {
        observer.time = 0.0;
        shader.needsUpdate = true;
    }};
    folder.add(resetSimObj, 'reset simulation').name('â†º reset simulation');
    //folder.open();


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


