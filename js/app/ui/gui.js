// Role: dat.GUI configuration — builds every control panel and folder, wires up
//       all parameter bindings, manages conditional visibility for jet / torus /
//       planet controls, and defines applyBlackHolePreset() to bulk-apply preset
//       data from presets.js. Must be called after init() has created the scene.

function setupGUI() {

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

    var gui = new dat.GUI({ width: 360 });

    // Recursively update all dat.GUI controllers to reflect programmatic changes.
    // dat.GUI does NOT auto-sync the displayed value when the bound property
    // is changed from code — you must call updateDisplay() on each controller.
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
    }
    refreshAllControllersGlobal = refreshAllControllers;

    function setGuiRowClass(guiEl, klass) {
        $(guiEl.domElement).parent().parent().addClass(klass);
    }

    function addHelpText(controller, text) {
        if (!text) return;
        // Use the parent .cr row for the tooltip
        var row = $(controller.domElement).closest('.cr')[0];
        if (row) row.title = text;
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
        if (cfg.onChange) c = c.onChange(cfg.onChange) || c;
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

    var updateDependentVisibility = function() {};

    function applyKerrMode(mode) {
        // Re-apply quality preset with new mode (preset values depend on kerr_mode)
        applyQualityPresetInternal(p.quality);
        hint.text('Solver mode: ' + mode.replace(/_/g, ' '));
        hint.stop(true, true).fadeIn(120).delay(900).fadeOut(350);
    }

    function applyQualityPresetInternal(value) {
        var isKerr = (p.kerr_mode === 'realtime_full_kerr_core');
        $('.planet-controls').show();

        switch(value) {
        case 'fast':
            if (isKerr) {
                // Kerr mode needs more steps for accuracy, but still faster preset
                p.n_steps = 200;
                p.sample_count = 2;
                p.max_revolutions = 2.5;
                p.rk4_integration = true;
            } else {
                p.n_steps = 40;
                p.sample_count = 1;
                p.max_revolutions = 1.5;
                p.rk4_integration = false;
            }
            p.cinematic_tonemap = true;
            $('.planet-controls').hide();
            p.resolution_scale = 1.0;
            p.taa_enabled = false;
            p.taa.history_weight = 0.88;
            p.taa.clip_box = 0.06;
            p.taa.motion_rejection = 8.0;
            p.taa.max_camera_delta = 0.08;
            p.taa.motion_clip_scale = 0.6;
            break;
        case 'medium':
            if (isKerr) {
                p.n_steps = 400;
                p.sample_count = 3;
                p.max_revolutions = 3.0;
                p.rk4_integration = true;
            } else {
                p.n_steps = 100;
                p.sample_count = 1;
                p.max_revolutions = 2.0;
                p.rk4_integration = false;
            }
            p.cinematic_tonemap = true;
            p.resolution_scale = 1.0;
            p.taa_enabled = false;
            p.taa.history_weight = 0.88;
            p.taa.clip_box = 0.06;
            p.taa.motion_rejection = 8.0;
            p.taa.max_camera_delta = 0.08;
            p.taa.motion_clip_scale = 0.6;
            break;
        case 'high':
            if (isKerr) {
                p.n_steps = 520;
                p.sample_count = 4;
                p.max_revolutions = 3.5;
                p.rk4_integration = true;
            } else {
                p.n_steps = 320;
                p.sample_count = 4;
                p.max_revolutions = 3.2;
                p.rk4_integration = true;
            }
            p.cinematic_tonemap = true;
            p.resolution_scale = 1.0;
            p.taa_enabled = false;
            p.taa.history_weight = 0.88;
            p.taa.clip_box = 0.06;
            p.taa.motion_rejection = 8.0;
            p.taa.max_camera_delta = 0.08;
            p.taa.motion_clip_scale = 0.6;
            break;
        case 'mobile':
            if (isKerr) {
                p.n_steps = 120;
                p.sample_count = 1;
                p.max_revolutions = 2.0;
                p.rk4_integration = false;
            } else {
                p.n_steps = 28;
                p.sample_count = 1;
                p.max_revolutions = 1.4;
                p.rk4_integration = false;
            }
            p.cinematic_tonemap = true;
            p.resolution_scale = 0.55;
            p.taa_enabled = true;
            p.taa.history_weight = 0.82;
            p.taa.clip_box = 0.08;
            p.taa.motion_rejection = 10.0;
            p.taa.max_camera_delta = 0.07;
            p.taa.motion_clip_scale = 0.8;
            $('.planet-controls').hide();
            break;
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
    }

    var qualityLabels = {
        'Fast (preview)': 'fast',
        'Medium': 'medium',
        'High': 'high',
        'Mobile (low power)': 'mobile'
    };

    var kerrModeLabels = {
        'Fast (approximate lensing)': 'fast',
        'Realtime Kerr core': 'realtime_full_kerr_core'
    };

    // ─── Real Black Hole Presets ─────────────────────────────────────────────────
    // applyBlackHolePreset is assigned below, after visibility helpers are defined.
    var presetObj = { preset: 'Default' };
    var applyBlackHolePreset;
    var presetsFolder = gui.addFolder('Black Hole Presets');
    addControl(presetsFolder, presetObj, 'preset', {
        name: 'preset',
        options: ['Custom', 'Default', 'M87*', 'Sgr A*', 'Cygnus X-1', 'GRS 1915+105', 'Gargantua (Interstellar visuals)', 'Schwarzschild'],
        onChange: function(val) { if (applyBlackHolePreset) applyBlackHolePreset(val); },
        help: 'Load physically motivated parameters for a real black hole. Overrides all settings below. Select Custom to tweak freely.'
    });
    presetsFolder.open();
    // ─────────────────────────────────────────────────────────────────────────────

    var renderFolder = gui.addFolder('Rendering');
    addControl(renderFolder, p, 'quality', {
        options: qualityLabels,
        name: 'quality preset',
        onChange: applyQualityPreset,
        help: 'Global render preset. High preset uses more ray-marching steps and samples.'
    });
    addControl(renderFolder, p, 'kerr_mode', {
        options: kerrModeLabels,
        name: 'solver mode',
        onChange: applyKerrMode,
        help: 'Fast = fastest approximate. Realtime Kerr core = accurate full GR with good performance.'
    });

    addControl(renderFolder, p, 'n_steps', {
        min: 20,
        max: 1400,
        step: 1,
        name: 'ray steps',
        onChange: updateShader,
        help: 'More steps improve thin features and strong lensing, but reduce FPS.'
    });
    addControl(renderFolder, p, 'sample_count', {
        min: 1,
        max: 12,
        step: 1,
        name: 'samples / pixel',
        onChange: updateShader,
        help: 'Supersampling for anti-aliasing and smoother edges. Higher values are slower.'
    });
    addControl(renderFolder, p, 'max_revolutions', {
        min: 1.0,
        max: 8.0,
        step: 0.1,
        name: 'max orbit turns',
        onChange: updateShader,
        help: 'How many wrapped photon turns are traced before escape/capture cutoff.'
    });
    addControl(renderFolder, p, 'resolution_scale', {
        min: 0.35,
        max: 2.0,
        step: 0.05,
        name: 'resolution scale',
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
        onChange: function() {
            if (typeof resetTemporalAAHistory === 'function') {
                resetTemporalAAHistory();
            }
            updateDependentVisibility();
            shader.needsUpdate = true;
        },
        help: 'Temporal anti-aliasing with aggressive history clamping to avoid ghosting.'
    });
    var taaHistoryWeightCtrl = addControl(renderFolder, p.taa, 'history_weight', {
        min: 0.0,
        max: 0.98,
        step: 0.01,
        name: 'taa history weight',
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
        onChange: updateShader,
        help: 'Higher-order integration for better stability in curved trajectories.'
    });
    renderFolder.open();

    applyQualityPreset(p.quality);

    // Custom scroll handler to control observer distance
    // Placed here so we have access to refreshAllControllers and distanceController
    renderer.domElement.addEventListener( 'wheel', function(e) {
        e.preventDefault();
        var delta = e.deltaY > 0 ? 1.15 : 0.87; // zoom out / zoom in
        var newDist = p.observer.distance * delta;
        newDist = Math.max(1.5, Math.min(30, newDist));
        p.observer.distance = newDist;
        updateCamera();
        shader.needsUpdate = true;
        refreshAllControllers();
    }, { passive: false } );

    var diskFolder = gui.addFolder('Accretion disk');
    addControl(diskFolder, p, 'accretion_disk', {
        name: 'enabled',
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Toggles thermal emission from the accretion disk.'
    });
    var accretionModeCtrl = addControl(diskFolder, p, 'accretion_mode', {
        name: 'accretion type',
        options: ['thin_disk', 'thick_torus', 'slim_disk'],
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Thin disk: Novikov-Thorne (quasars/XRBs). Thick torus: ADAF/RIAF (M87*/Sgr A*). Slim disk: super-Eddington.'
    });
    var diskTempCtrl = addControl(diskFolder, p, 'disk_temperature', {
        min: DISK_TEMPERATURE_MIN,
        max: DISK_TEMPERATURE_MAX,
        step: 1,
        name: 'temperature (K)',
        onChange: updateUniformsLive,
        help: 'Rest-frame disk color temperature before relativistic shifts.'
    });
    var torusCenterCtrl = addControl(diskFolder, p.torus, 'r0', {
        min: 1.5,
        max: 10.0,
        step: 0.1,
        name: 'torus center r',
        onChange: updateUniformsLive,
        help: 'Center radius of the torus in r_s units (thick torus mode only).'
    });
    var torusHRCtrl = addControl(diskFolder, p.torus, 'h_ratio', {
        min: 0.1,
        max: 1.0,
        step: 0.01,
        name: 'torus H/R',
        onChange: updateUniformsLive,
        help: 'Height-to-radius ratio of the torus cross-section (thick torus mode only).'
    });
    var torusFalloffCtrl = addControl(diskFolder, p.torus, 'radial_falloff', {
        min: 0.5,
        max: 5.0,
        step: 0.1,
        name: 'radial falloff',
        onChange: updateUniformsLive,
        help: 'Power-law index for radial emissivity decay. Physical ADAF: 2-4. Lower = flatter profile, higher = more concentrated.'
    });
    var torusOpacityCtrl = addControl(diskFolder, p.torus, 'opacity', {
        min: 0.001,
        max: 0.15,
        step: 0.001,
        name: 'opacity',
        onChange: updateUniformsLive,
        help: 'Absorption coefficient. Low = optically thin (ADAF). Higher = more self-shielding and defined torus surface.'
    });
    var torusOuterCtrl = addControl(diskFolder, p.torus, 'outer_radius', {
        min: 1.5,
        max: 8.0,
        step: 0.1,
        name: 'outer extent',
        onChange: updateUniformsLive,
        help: 'Outer edge multiplier relative to torus center r0. Controls how far the torus extends.'
    });

    // ─── Slim disk controls ───────────────────────────────
    var slimHRCtrl = addControl(diskFolder, p.slim, 'h_ratio', {
        min: 0.05,
        max: 0.5,
        step: 0.01,
        name: 'slim H/R',
        onChange: updateUniformsLive,
        help: 'Base height-to-radius ratio of the slim disk. Higher = geometrically thicker.'
    });
    var slimOpacityCtrl = addControl(diskFolder, p.slim, 'opacity', {
        min: 0.1,
        max: 2.0,
        step: 0.01,
        name: 'slim opacity',
        onChange: updateUniformsLive,
        help: 'Absorption coefficient. Higher = more opaque, surface-like. Lower = more translucent, volumetric.'
    });
    var slimPuffCtrl = addControl(diskFolder, p.slim, 'puff_factor', {
        min: 0.0,
        max: 6.0,
        step: 0.1,
        name: 'ISCO puff',
        onChange: updateUniformsLive,
        help: 'How much the disk puffs up near the ISCO due to radiation pressure. Higher = more pronounced thickening.'
    });

    // Store GUI row elements for conditional visibility
    var torusRows = [torusCenterCtrl, torusHRCtrl, torusFalloffCtrl, torusOpacityCtrl, torusOuterCtrl];
    var slimRows = [slimHRCtrl, slimOpacityCtrl, slimPuffCtrl];

    // Apply initial accretion visibility immediately, even before the global
    // dependency updater is assigned at the end of setupGUI().
    setControlVisible(accretionModeCtrl, !!p.accretion_disk);
    setControlVisible(diskTempCtrl, !!p.accretion_disk);
    setControlsVisible(torusRows, !!p.accretion_disk && p.accretion_mode === 'thick_torus');
    setControlsVisible(slimRows, !!p.accretion_disk && p.accretion_mode === 'slim_disk');

    diskFolder.open();

    // ─── Relativistic Jets ───────────────────────────────
    var jetFolder = gui.addFolder('Relativistic jets');
    addControl(jetFolder, p.jet, 'enabled', {
        name: 'enabled',
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Bipolar relativistic jets along the spin axis (Blandford-Znajek mechanism).'
    });
    var jetModeCtrl = addControl(jetFolder, p.jet, 'mode', {
        name: 'jet model',
        options: ['simple', 'physical'],
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Simple: smooth parabolic jet. Physical: GRMHD-calibrated model with spine/sheath, reconfinement knots, corona base, disk occultation.'
    });
    var jetAngleCtrl = addControl(jetFolder, p.jet, 'half_angle', {
        min: 1.0,
        max: 25.0,
        step: 0.5,
        name: 'half-angle (°)',
        onChange: updateUniformsLive,
        help: 'Opening half-angle of the jet. Typical AGN jets: 2-7°. Parabolic collimation narrows the beam with distance.'
    });
    var jetLorentzCtrl = addControl(jetFolder, p.jet, 'lorentz_factor', {
        min: 1.1,
        max: 20.0,
        step: 0.1,
        name: 'Lorentz Γ',
        onChange: updateUniformsLive,
        help: 'Bulk Lorentz factor. Γ~2-5 for AGN jets, Γ~100+ for GRBs. Controls relativistic beaming.'
    });
    var jetBrightCtrl = addControl(jetFolder, p.jet, 'brightness', {
        min: 0.05,
        max: 3.0,
        step: 0.01,
        name: 'brightness',
        onChange: updateUniformsLive,
        help: 'Overall jet synchrotron emission strength.'
    });
    var jetLengthCtrl = addControl(jetFolder, p.jet, 'length', {
        min: 5.0,
        max: 60.0,
        step: 1.0,
        name: 'length (r_s)',
        onChange: updateUniformsLive,
        help: 'Visible jet length in Schwarzschild radii.'
    });
    var jetMagCtrl = addControl(jetFolder, p.jet, 'magnetization', {
        min: 1.0,
        max: 50.0,
        step: 0.5,
        name: 'σ (magnetization)',
        onChange: updateUniformsLive,
        help: 'Plasma magnetization at jet base (σ = B²/4πρc²). Higher σ = more magnetically dominated, narrower spine. MAD jets: σ ~ 10-30.'
    });
    var jetKnotCtrl = addControl(jetFolder, p.jet, 'knot_spacing', {
        min: 2.0,
        max: 15.0,
        step: 0.5,
        name: 'knot spacing',
        onChange: updateUniformsLive,
        help: 'Spacing of reconfinement shock knots (in r_s). Observed in M87 (HST-1), 3C 273. Set high to minimize knots.'
    });
    var jetCoronaCtrl = addControl(jetFolder, p.jet, 'corona_brightness', {
        min: 0.0,
        max: 5.0,
        step: 0.1,
        name: 'corona glow',
        onChange: updateUniformsLive,
        help: 'Brightness of the jet-corona connection at the base. Hot plasma where the funnel meets the inner accretion flow.'
    });
    var jetBaseWidthCtrl = addControl(jetFolder, p.jet, 'base_width', {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        name: 'base width',
        onChange: updateUniformsLive,
        help: 'Controls how wide the jet funnel is at the base. 0 = narrow collimated base, 1 = wide split-monopole funnel.'
    });
    var jetCoronaExtentCtrl = addControl(jetFolder, p.jet, 'corona_extent', {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        name: 'corona extent',
        onChange: updateUniformsLive,
        help: 'Radial spread of the corona base emission. 0 = wide spread, 1 = tightly concentrated near the jet axis.'
    });

    // Store jet row groups for conditional visibility
    var jetCommonCtrls = [jetModeCtrl, jetAngleCtrl, jetLorentzCtrl, jetBrightCtrl, jetLengthCtrl];
    var jetPhysicalCtrls = [jetMagCtrl, jetKnotCtrl, jetCoronaCtrl, jetBaseWidthCtrl, jetCoronaExtentCtrl];

    applyBlackHolePreset = function(name) {
        if (name === 'Custom') return;
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
        p.observer.distance          = preset.observer.distance;
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
        onChange: updateUniformsLive,
        help: 'Dimensionless spin. Positive = prograde disk, negative = retrograde.'
    });
    var spinStrengthCtrl = addControl(spinFolder, p.black_hole, 'spin_strength', {
        min: 0.0,
        max: 1.4,
        step: 0.01,
        name: 'shadow squeeze',
        onChange: updateUniformsLive,
        help: 'Visual strength multiplier for rotation-induced asymmetry in the shadow.'
    });
    spinFolder.open();

    var lookFolder = gui.addFolder('Look');
    addControl(lookFolder, p.look, 'tonemap_mode', {
        options: { 'ACES Filmic': 0, 'AGX (Black Hole)': 1, 'Scientific (Log)': 2 },
        name: 'tonemapper',
        onChange: updateUniformsLive,
        help: 'ACES: classic cinematic. AGX: better saturation handling for extreme HDR. Scientific: logarithmic false-color like EHT papers.'
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
        onChange: updateUniformsLive,
        help: 'Brightness multiplier applied to the accretion disk emission.'
    });
    var lookGlowCtrl = addControl(lookFolder, p.look, 'glow', {
        min: 0.0,
        max: 2.0,
        step: 0.01,
        name: 'inner glow',
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

    // ─── Post-processing (bloom) ────────────────────────────────────────────
    var ppFolder = gui.addFolder('Post-processing');
    addControl(ppFolder, p.bloom, 'enabled', {
        name: 'bloom',
        onChange: function() {
            updateDependentVisibility();
            shader.needsUpdate = true;
        },
        help: 'Physically-based bloom simulating optical diffraction and lens glow. Bright regions of the accretion disk bleed light into surrounding pixels.'
    });
    var bloomStrengthCtrl = addControl(ppFolder, p.bloom, 'strength', {
        min: 0.0,
        max: 2.0,
        step: 0.01,
        name: 'bloom strength',
        onChange: function() { shader.needsUpdate = true; },
        help: 'Overall intensity of the bloom glow added to the image.'
    });
    var bloomThresholdCtrl = addControl(ppFolder, p.bloom, 'threshold', {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        name: 'bloom threshold',
        onChange: function() { shader.needsUpdate = true; },
        help: 'Minimum pixel brightness that contributes to bloom. Lower values bloom more of the disk.'
    });
    var bloomRadiusCtrl = addControl(ppFolder, p.bloom, 'radius', {
        min: 0.0,
        max: 1.0,
        step: 0.01,
        name: 'bloom radius',
        onChange: function() { shader.needsUpdate = true; },
        help: 'Controls the width of the bloom halo. Higher values give wider, more diffuse glow matching real optical PSFs.'
    });
    ppFolder.open();

    var folder = gui.addFolder('Observer');
    addControl(folder, p.observer, 'motion', {
        name: 'orbital motion',
        help: 'When enabled, the observer follows a circular orbit around the black hole.',
        onChange: function(motion) {
        updateCamera();
        updateShader();
        if (motion) {
            hint.text('Moving observer; drag to rotate camera');
        } else {
            hint.text('Stationary observer; drag to orbit around');
        }
        hint.fadeIn();
        }
    });
    distanceController = addControl(folder, p.observer, 'distance', {
        min: 1.5,
        max: 30,
        step: 0.1,
        name: 'distance',
        onChange: function() {
            updateCamera();
            shader.needsUpdate = true;
        },
        help: 'Observer distance from the black hole center in Schwarzschild-radius units.'
    });
    var observerActions = {
        reset_camera: function() {
            if (cameraControls) {
                cameraControls.reset();
                cameraControls.target.set(0, 0, 0);
            }
            if (cameraPan) {
                cameraPan.set(0, 0);
            }
            updateCamera();
            shader.needsUpdate = true;
            hint.text('Camera reset');
            hint.stop(true, true).fadeIn(120).delay(700).fadeOut(320);
        }
    };
    addControl(folder, observerActions, 'reset_camera', {
        name: 'reset camera',
        help: 'Reset camera pan, tilt/roll, and orbit orientation.'
    });
    folder.open();

    folder = gui.addFolder('Planet');
    addControl(folder, p.planet, 'enabled', {
        name: 'enabled',
        help: 'Adds an orbiting planet used as a reference object.',
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        }
    });
    var planetDistanceCtrl = addControl(folder, p.planet, 'distance', {
        min: 1.5,
        step: 0.1,
        name: 'distance',
        onChange: updateUniforms,
        help: 'Orbital radius of the planet.'
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
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Applies relativistic intensity boosting/dimming.'
    });
    var physicalBeamingCtrl = addControl(folder, p, 'physical_beaming', {
        name: 'physical (D³ Liouville)',
        onChange: function() {
            updateDependentVisibility();
            updateShader();
        },
        help: 'Uses physically motivated Liouville transfer scaling instead of cinematic curve.'
    });
    addControl(folder, p, 'doppler_shift', {
        name: 'doppler shift (color)',
        onChange: updateShader,
        help: 'Shifts observed spectrum by red/blue shift factors.'
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
    //folder.open();

    updateDependentVisibility = function() {
        var diskEnabled = !!p.accretion_disk;
        var thinDiskEnabled = diskEnabled && p.accretion_mode === 'thin_disk';
        var thickTorusEnabled = diskEnabled && p.accretion_mode === 'thick_torus';
        var slimDiskEnabled = diskEnabled && p.accretion_mode === 'slim_disk';

        setControlVisible(accretionModeCtrl, diskEnabled);
        setControlVisible(diskTempCtrl, diskEnabled);
        setControlsVisible(torusRows, thickTorusEnabled);
        setControlsVisible(slimRows, slimDiskEnabled);

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

    // ─── Freefall Dive Panel ────────────────────────────────────────────────────
    // Creates a separate animation panel positioned to the left of the dat.GUI
    // controls.  The observer plunges radially from rest at infinity along a
    // Schwarzschild geodesic.  Inside the event horizon the shader switches to
    // interior coordinate mode: the Binet equation is integrated past u = 1
    // (r < r_s) so backward-traced photons can exit the horizon and reveal the
    // external universe as a shrinking, violently lensed window.
    (function setupDivePanel() {
        var panel = document.createElement('div');
        panel.id = 'dive-panel';
        panel.innerHTML =
            '<div class="dive-title">FREEFALL DIVE</div>' +
            '<div class="dive-desc">Radial plunge into the black hole interior. ' +
            'Rays are traced through the event horizon with the full Schwarzschild ' +
            'geodesic equation &mdash; no approximations.</div>' +
            '<button id="dive-start-btn" class="dive-btn dive-btn-start">' +
                '\u25b6 START DIVE</button>' +
            '<div class="dive-control-row">' +
                '<label>Fall speed</label>' +
                '<input type="range" id="dive-speed" min="0.01" max="5.0" ' +
                    'step="0.01" value="1.0">' +
                '<span id="dive-speed-val">1.0\u00d7</span>' +
            '</div>' +
            '<div class="dive-control-row dive-cinematic-row">' +
                '<label for="dive-cinematic">Auto-speed</label>' +
                '<input type="checkbox" id="dive-cinematic">' +
                '<span class="dive-cinematic-hint">Slow near photon sphere &amp; horizon</span>' +
            '</div>' +
            '<div id="dive-horizon-track" class="dive-horizon-track">' +
                '<div id="dive-horizon-bar" class="dive-horizon-fill outside">' +
                '</div>' +
                '<div class="dive-horizon-label">Event Horizon</div>' +
            '</div>' +
            '<div class="dive-readout">' +
                '<div id="dive-radius" class="dive-metric">' +
                    'r = ' + p.observer.distance.toFixed(2) +
                    ' r<sub>s</sub></div>' +
                '<div id="dive-velocity" class="dive-metric">v = 0.000 c</div>' +
                '<div id="dive-status" class="dive-status ready">Ready</div>' +
            '</div>' +
            '<button id="dive-reset-btn" class="dive-btn dive-btn-reset" ' +
                'disabled>\u21ba RESET</button>';
        document.body.appendChild(panel);

        document.getElementById('dive-start-btn').addEventListener('click',
            function() { startDive(); });
        document.getElementById('dive-reset-btn').addEventListener('click',
            function() { resetDive(); });
        document.getElementById('dive-speed').addEventListener('input',
            function() {
                diveState.speed = parseFloat(this.value);
                var v = diveState.speed;
                document.getElementById('dive-speed-val').textContent =
                    (v < 0.1 ? v.toFixed(2) : v.toFixed(1)) + '\u00d7';
            });

        document.getElementById('dive-cinematic').addEventListener('change',
            function() { diveState.cinematic = this.checked; });

        // ── Clickable / draggable progress bar ──────────────────────────
        var track = document.getElementById('dive-horizon-track');
        var dragging = false;
        function handleTrackSeek(e) {
            if (!diveState.active && !diveState.reachedSingularity) return;
            var rect = track.getBoundingClientRect();
            var clientX = e.touches ? e.touches[0].clientX : e.clientX;
            var x = Math.max(0, Math.min(clientX - rect.left, rect.width));
            var progress = x / rect.width;
            var startR = Math.max(diveState.prevDistance, 1);
            var targetR = startR * (1.0 - progress);
            seekDive(targetR);
        }
        track.addEventListener('mousedown', function(e) {
            dragging = true; handleTrackSeek(e); e.preventDefault();
        });
        document.addEventListener('mousemove', function(e) {
            if (dragging) handleTrackSeek(e);
        });
        document.addEventListener('mouseup', function() { dragging = false; });
        track.addEventListener('touchstart', function(e) {
            handleTrackSeek(e); e.preventDefault();
        });
        track.addEventListener('touchmove', function(e) {
            handleTrackSeek(e); e.preventDefault();
        });

        // ── 3-D axes orientation gizmo ──────────────────────────────────
        var gizmo = document.createElement('div');
        gizmo.id = 'axes-gizmo-container';
        gizmo.innerHTML = '<canvas id="axes-gizmo" width="80" height="80"></canvas>';
        document.body.appendChild(gizmo);
    })();
    // ─────────────────────────────────────────────────────────────────────────────

}
