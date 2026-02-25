"use strict"
/*global THREE, SHADER_LOADER, Mustache, Stats, Detector, $, dat:false */
/*global document, window, setTimeout, requestAnimationFrame:false */
/*global ProceduralTextures:false */

if ( ! Detector.webgl ) Detector.addGetWebGLMessage();

var DISK_TEMPERATURE_MIN = 4500.0;
var DISK_TEMPERATURE_MAX = 30000.0;

function Observer() {
    this.position = new THREE.Vector3(10,0,0);
    this.velocity = new THREE.Vector3(0,1,0);
    this.orientation = new THREE.Matrix3();
    this.time = 0.0;
}

function formatThousands(value) {
    return Math.round(value).toString().replace(/\B(?=(\d{3})+(?!\d))/g, " ");
}

Observer.prototype.orbitalFrame = function() {

    //var orbital_y = observer.velocity.clone().normalize();
    var orbital_y = (new THREE.Vector3())
        .subVectors(observer.velocity.clone().normalize().multiplyScalar(4.0),
            observer.position).normalize();

    var orbital_z = (new THREE.Vector3())
        .crossVectors(observer.position, orbital_y).normalize();
    var orbital_x = (new THREE.Vector3()).crossVectors(orbital_y, orbital_z);


    return (new THREE.Matrix4()).makeBasis(
        orbital_x,
        orbital_y,
        orbital_z
    ).linearPart();
};

Observer.prototype.move = function(dt) {

    dt *= shader.parameters.time_scale;

    var r;
    var v = 0;

    // motion on a pre-defined cirular orbit
    if (shader.parameters.observer.motion) {

        r = shader.parameters.observer.distance;
        v =  1.0 / Math.sqrt(2.0*(r-1.0));
        // Convert local velocity to coordinate angular velocity
        // Ω = v·sqrt(1-r_s/r)/r = 1/sqrt(2r³) for circular Schwarzschild orbit
        var ang_vel = v * Math.sqrt(1.0 - 1.0/r) / r;
        var angle = this.time * ang_vel;

        var s = Math.sin(angle), c = Math.cos(angle);

        this.position.set(c*r, s*r, 0);
        this.velocity.set(-s*v, c*v, 0);

        var alpha = degToRad(shader.parameters.observer.orbital_inclination);
        var orbit_coords = (new THREE.Matrix4()).makeRotationY(alpha);

        this.position.applyMatrix4(orbit_coords);
        this.velocity.applyMatrix4(orbit_coords);
    }
    else {
        r = this.position.length();
    }

    if (shader.parameters.gravitational_time_dilation) {
        if (v > 0) {
            // Circular orbit: combined gravitational + kinematic time dilation
            // dτ/dt = sqrt(1 - 3M/r) = sqrt(1 - 3/(2r)) for r_s = 1, M = 0.5
            dt = dt / Math.sqrt(Math.max(1.0 - 1.5/r, 0.001));
        } else {
            // Stationary observer: gravitational time dilation only
            // dτ/dt = sqrt(1 - r_s/r) = sqrt(1 - 1/r)
            dt = dt / Math.sqrt(Math.max(1.0 - 1.0/r, 0.001));
        }
    }

    this.time += dt;
};

var container, stats;
var camera, scene, renderer, cameraControls, shader = null;
var observer = new Observer();
var distanceController = null;
var refreshAllControllersGlobal = null; // Will be set in setupGUI
var effectLabels = {
    spin: null,
    temperature: null
};

function updateEffectLabels() {
    if (!shader || !effectLabels.spin || !effectLabels.temperature) return;

    var spinPercent = shader.parameters.black_hole.spin * 100.0;
    effectLabels.spin.textContent = 'a/M = ' + spinPercent.toFixed(1) + '%';
    effectLabels.temperature.textContent =
        'T = ' + formatThousands(shader.parameters.disk_temperature) + ' K';

    if (shader.parameters.black_hole.spin_enabled) {
        effectLabels.spin.classList.remove('is-disabled');
    } else {
        effectLabels.spin.classList.add('is-disabled');
    }
}

function Shader(mustacheTemplate) {
    // Compile-time shader parameters
    this.parameters = {
        n_steps: 100,
        sample_count: 1,
        max_revolutions: 2.0,
        rk4_integration: false,
        cinematic_tonemap: true,
        quality: 'high',
        kerr_mode: 'realtime_full_kerr_core',
        accretion_disk: true,
        accretion_mode: 'thin_disk',
        disk_temperature: 10000.0,
        torus: {
            r0: 4.0,
            h_ratio: 0.45
        },
        jet: {
            enabled: false,
            mode: 'simple',
            half_angle: 5.0,
            lorentz_factor: 3.0,
            brightness: 1.2,
            length: 30.0,
            magnetization: 10.0,
            knot_spacing: 6.0,
            corona_brightness: 1.5
        },
        black_hole: {
            spin_enabled: true,
            spin: 0.90,
            spin_strength: 1.0
        },
        look: {
            tonemap_mode: 1,
            exposure: 1.0,
            disk_gain: 1.0,
            glow: 0.0,
            doppler_boost: 1.0,
            aberration_strength: 1.0,
            star_gain: 0.5,
            galaxy_gain: 0.5
        },
        planet: {
            enabled: true,
            distance: 14.0,
            radius: 0.4
        },
        lorentz_contraction: true,
        gravitational_time_dilation: true,
        aberration: true,
        beaming: true,
        physical_beaming: true,
        doppler_shift: true,
        light_travel_time: true,
        time_scale: 1.0,
        observer: {
            motion: true,
            distance: 11.0,
            orbital_inclination: -10
        },

        planetEnabled: function() {
            return this.planet.enabled && this.quality !== 'fast';
        },

        observerMotion: function() {
            return this.observer.motion;
        }
    };
    var that = this;
    this.needsUpdate = false;

    this.hasMovingParts = function() {
        return this.parameters.planet.enabled || this.parameters.observer.motion;
    };

    this.compile = function() {
        that.parameters.kerr_fast_mode = (that.parameters.kerr_mode === 'fast');
        that.parameters.kerr_full_core = (that.parameters.kerr_mode === 'realtime_full_kerr_core');

        var accMode = that.parameters.accretion_mode;
        var diskOn = that.parameters.accretion_disk;
        that.parameters.accretion_thin_disk = diskOn && (accMode === 'thin_disk');
        that.parameters.accretion_thick_torus = diskOn && (accMode === 'thick_torus');
        that.parameters.accretion_slim_disk = diskOn && (accMode === 'slim_disk');
        that.parameters.jet_enabled = that.parameters.jet.enabled;
        that.parameters.jet_simple = that.parameters.jet.enabled && (that.parameters.jet.mode === 'simple');
        that.parameters.jet_physical = that.parameters.jet.enabled && (that.parameters.jet.mode === 'physical');

        return Mustache.render(mustacheTemplate, that.parameters);
    };
}

function degToRad(a) { return Math.PI * a / 180.0; }

(function(){
    var textures = {};

    function whenLoaded() {
        init(textures);
        $('#loader').hide();
        $('.initially-hidden').removeClass('initially-hidden');
        animate();
    }

    function checkLoaded() {
        if (shader === null) return;
        for (var key in textures) if (textures[key] === null) return;
        whenLoaded();
    }

    SHADER_LOADER.load(function(shaders) {
        shader = new Shader(shaders.raytracer.fragment);
        checkLoaded();
    });

    var texLoader = new THREE.TextureLoader();
    function loadTexture(symbol, filename, interpolation) {
        textures[symbol] = null;
        texLoader.load(filename, function(tex) {
            tex.magFilter = interpolation;
            tex.minFilter = interpolation;
            textures[symbol] = tex;
            checkLoaded();
        });
    }

    loadTexture('galaxy', 'img/milkyway.jpg', THREE.NearestFilter);
    loadTexture('spectra', 'img/spectra.png', THREE.LinearFilter);
    loadTexture('moon', 'img/beach-ball.png', THREE.LinearFilter);
    loadTexture('stars', 'img/stars.png', THREE.LinearFilter);
})();

var updateUniforms;

function init(textures) {

    container = document.createElement( 'div' );
    document.body.appendChild( container );

    scene = new THREE.Scene();

    var geometry = new THREE.PlaneBufferGeometry( 2, 2 );

    var uniforms = {
        time: { type: "f", value: 0 },
        resolution: { type: "v2", value: new THREE.Vector2() },
        cam_pos: { type: "v3", value: new THREE.Vector3() },
        cam_x: { type: "v3", value: new THREE.Vector3() },
        cam_y: { type: "v3", value: new THREE.Vector3() },
        cam_z: { type: "v3", value: new THREE.Vector3() },
        cam_vel: { type: "v3", value: new THREE.Vector3() },

        planet_distance: { type: "f" },
        planet_radius: { type: "f" },
        disk_temperature: { type: "f", value: 10000.0 },
        accretion_inner_r: { type: "f", value: 3.0 },
        bh_spin: { type: "f", value: 0.90 },
        bh_spin_strength: { type: "f", value: 1.0 },
        bh_rotation_enabled: { type: "f", value: 1.0 },
        look_exposure: { type: "f", value: 1.0 },
        look_disk_gain: { type: "f", value: 1.0 },
        look_glow: { type: "f", value: 0.0 },
        look_doppler_boost: { type: "f", value: 1.0 },
        look_aberration_strength: { type: "f", value: 1.0 },
        look_star_gain: { type: "f", value: 1.0 },
        look_galaxy_gain: { type: "f", value: 1.0 },
        look_tonemap_mode: { type: "f", value: 0.0 },

        torus_r0: { type: "f", value: 4.0 },
        torus_h_ratio: { type: "f", value: 0.45 },

        jet_half_angle: { type: "f", value: 5.0 },
        jet_lorentz: { type: "f", value: 3.0 },
        jet_brightness: { type: "f", value: 1.2 },
        jet_length: { type: "f", value: 30.0 },
        jet_magnetization: { type: "f", value: 10.0 },
        jet_knot_spacing: { type: "f", value: 6.0 },
        jet_corona_brightness: { type: "f", value: 1.5 },

        star_texture: { type: "t", value: textures.stars },
        galaxy_texture: { type: "t", value: textures.galaxy },
        planet_texture: { type: "t", value: textures.moon },
        spectrum_texture: { type: "t", value: textures.spectra }
    };

    // Calculate ISCO radius using Bardeen-Press-Teukolsky formula
    // chi is dimensionless spin parameter (-1 to 1)
    // Returns ISCO in units of Schwarzschild radius (r_s = 1)
    function calculateISCO(chi, isPrograde) {
        var chi2 = chi * chi;
        var cbrt_1_minus_chi2 = Math.pow(Math.max(1 - chi2, 0), 1/3);
        var cbrt_1_plus_chi = Math.pow(1 + Math.abs(chi), 1/3);
        var cbrt_1_minus_chi = Math.pow(Math.max(1 - Math.abs(chi), 0), 1/3);
        
        var Z1 = 1 + cbrt_1_minus_chi2 * (cbrt_1_plus_chi + cbrt_1_minus_chi);
        var Z2 = Math.sqrt(3 * chi2 + Z1 * Z1);
        
        // Prograde orbits (co-rotating with black hole) have smaller ISCO
        // Retrograde orbits have larger ISCO
        var sign = isPrograde ? -1 : 1;
        var isco_rg = 3 + Z2 + sign * Math.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2));
        
        // Convert from gravitational radii (r_g = GM/c^2) to Schwarzschild radii (r_s = 2*r_g)
        // Since our units have r_s = 1, we have r_g = 0.5
        return isco_rg * 0.5;
    }

    updateUniforms = function() {
        uniforms.planet_distance.value = shader.parameters.planet.distance;
        uniforms.planet_radius.value = shader.parameters.planet.radius;
        uniforms.disk_temperature.value = shader.parameters.disk_temperature;
        
        // Calculate ISCO based on spin (prograde disk assumed)
        var spin = shader.parameters.black_hole.spin;
        var spinEnabled = shader.parameters.black_hole.spin_enabled;
        uniforms.accretion_inner_r.value = spinEnabled ? calculateISCO(spin, true) : 3.0;
        
        uniforms.bh_spin.value = shader.parameters.black_hole.spin;
        uniforms.bh_spin_strength.value = shader.parameters.black_hole.spin_strength;
        uniforms.bh_rotation_enabled.value = shader.parameters.black_hole.spin_enabled ? 1.0 : 0.0;
        uniforms.look_exposure.value = shader.parameters.look.exposure;
        uniforms.look_disk_gain.value = shader.parameters.look.disk_gain;
        uniforms.look_glow.value = shader.parameters.look.glow;
        uniforms.look_doppler_boost.value = shader.parameters.look.doppler_boost;
        uniforms.look_aberration_strength.value = shader.parameters.look.aberration_strength;
        uniforms.look_star_gain.value = shader.parameters.look.star_gain;
        uniforms.look_galaxy_gain.value = shader.parameters.look.galaxy_gain;
        uniforms.look_tonemap_mode.value = parseFloat(shader.parameters.look.tonemap_mode);

        uniforms.torus_r0.value = shader.parameters.torus.r0;
        uniforms.torus_h_ratio.value = shader.parameters.torus.h_ratio;

        uniforms.jet_half_angle.value = shader.parameters.jet.half_angle;
        uniforms.jet_lorentz.value = shader.parameters.jet.lorentz_factor;
        uniforms.jet_brightness.value = shader.parameters.jet.brightness;
        uniforms.jet_length.value = shader.parameters.jet.length;
        uniforms.jet_magnetization.value = shader.parameters.jet.magnetization;
        uniforms.jet_knot_spacing.value = shader.parameters.jet.knot_spacing;
        uniforms.jet_corona_brightness.value = shader.parameters.jet.corona_brightness;

        uniforms.resolution.value.x = renderer.domElement.width;
        uniforms.resolution.value.y = renderer.domElement.height;

        uniforms.time.value = observer.time;
        uniforms.cam_pos.value = observer.position;

        var e = observer.orientation.elements;

        uniforms.cam_x.value.set(e[0], e[1], e[2]);
        uniforms.cam_y.value.set(e[3], e[4], e[5]);
        uniforms.cam_z.value.set(e[6], e[7], e[8]);

        function setVec(target, value) {
            uniforms[target].value.set(value.x, value.y, value.z);
        }

        setVec('cam_pos', observer.position);
        setVec('cam_vel', observer.velocity);

        updateEffectLabels();
    };

    var material = new THREE.ShaderMaterial( {
        uniforms: uniforms,
        vertexShader: $('#vertex-shader').text(),
    });

    scene.updateShader = function() {
        material.fragmentShader = shader.compile();
        material.needsUpdate = true;
        shader.needsUpdate = true;
    };

    scene.updateShader();

    var mesh = new THREE.Mesh( geometry, material );
    scene.add( mesh );

    renderer = new THREE.WebGLRenderer({
        antialias: true,
        powerPreference: 'high-performance'
    });
    renderer.setPixelRatio( window.devicePixelRatio );
    container.appendChild( renderer.domElement );

    stats = new Stats();
    stats.domElement.style.position = 'absolute';
    stats.domElement.style.top = '0px';
    container.appendChild( stats.domElement );
    $(stats.domElement).addClass('hidden-phone');

    effectLabels.spin = document.getElementById('spin-label');
    effectLabels.temperature = document.getElementById('temperature-label');
    updateEffectLabels();

    // Orbit camera from three.js
    camera = new THREE.PerspectiveCamera( 45, window.innerWidth / window.innerHeight, 1, 80000 );
    initializeCamera(camera);

    cameraControls = new THREE.OrbitControls( camera, renderer.domElement );
    cameraControls.target.set( 0, 0, 0 );
    cameraControls.enableZoom = false; // We handle zoom manually for distance sync
    cameraControls.addEventListener( 'change', updateCamera );
    updateCamera();

    onWindowResize();

    window.addEventListener( 'resize', onWindowResize, false );

    setupGUI();
}

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
        if (cfg.min !== undefined) c.min(cfg.min);
        if (cfg.max !== undefined) c.max(cfg.max);
        if (cfg.step !== undefined) c.step(cfg.step);
        if (cfg.name) c.name(cfg.name);
        if (cfg.onChange) c.onChange(cfg.onChange);
        if (cfg.help) addHelpText(c, cfg.help);
        if (cfg.className) setGuiRowClass(c, cfg.className);
        return c;
    }

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
            break;
        }

        updateShader();
        refreshAllControllers();
    }

    function applyQualityPreset(value) {
        applyQualityPresetInternal(value);
    }

    var qualityLabels = {
        'Fast (preview)': 'fast',
        'Medium': 'medium',
        'High': 'high'
    };

    var kerrModeLabels = {
        'Fast (approximate lensing)': 'fast',
        'Realtime Kerr core': 'realtime_full_kerr_core'
    };

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
        onChange: updateShader,
        help: 'Toggles thermal emission from the accretion disk.'
    });
    addControl(diskFolder, p, 'accretion_mode', {
        name: 'accretion type',
        options: ['thin_disk', 'thick_torus', 'slim_disk'],
        onChange: function(mode) {
            updateAccretionModeVisibility(mode);
            updateShader();
        },
        help: 'Thin disk: Novikov-Thorne (quasars/XRBs). Thick torus: ADAF/RIAF (M87*/Sgr A*). Slim disk: super-Eddington.'
    });
    addControl(diskFolder, p, 'disk_temperature', {
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

    // Store GUI row elements for conditional visibility
    var torusRows = [torusCenterCtrl, torusHRCtrl];

    function updateAccretionModeVisibility(mode) {
        var showTorus = (mode === 'thick_torus');
        torusRows.forEach(function(ctrl) {
            $(ctrl.domElement).closest('li').css('display', showTorus ? '' : 'none');
        });
    }
    // Initialize visibility
    updateAccretionModeVisibility(p.accretion_mode);

    diskFolder.open();

    // ─── Relativistic Jets ───────────────────────────────
    var jetFolder = gui.addFolder('Relativistic jets');
    addControl(jetFolder, p.jet, 'enabled', {
        name: 'enabled',
        onChange: updateShader,
        help: 'Bipolar relativistic jets along the spin axis (Blandford-Znajek mechanism).'
    });
    var jetModeCtrl = addControl(jetFolder, p.jet, 'mode', {
        name: 'jet model',
        options: ['simple', 'physical'],
        onChange: function(mode) {
            updateJetModeVisibility(mode, p.jet.enabled);
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

    // Show/hide jet parameter controls based on jet.enabled and jet.mode
    var jetCommonCtrls = [jetModeCtrl, jetAngleCtrl, jetLorentzCtrl, jetBrightCtrl, jetLengthCtrl];
    var jetPhysicalCtrls = [jetMagCtrl, jetKnotCtrl, jetCoronaCtrl];
    function updateJetModeVisibility(mode, enabled) {
        jetCommonCtrls.forEach(function(ctrl) {
            $(ctrl.domElement).closest('li').css('display', enabled ? '' : 'none');
        });
        var showPhysical = enabled && (mode === 'physical');
        jetPhysicalCtrls.forEach(function(ctrl) {
            $(ctrl.domElement).closest('li').css('display', showPhysical ? '' : 'none');
        });
    }
    updateJetModeVisibility(p.jet.mode, p.jet.enabled);
    // Patch the enabled onChange to also update visibility
    jetFolder.__controllers[0].onChange(function(val) {
        updateJetModeVisibility(p.jet.mode, val);
        updateShader();
    });

    var spinFolder = gui.addFolder('Black hole');
    addControl(spinFolder, p.black_hole, 'spin_enabled', {
        name: 'rotation enabled',
        onChange: updateUniformsLive,
        help: 'Enables Kerr-like rotation effects. Disable for a Schwarzschild-style shadow.'
    });
    addControl(spinFolder, p.black_hole, 'spin', {
        min: -0.99,
        max: 0.99,
        step: 0.01,
        name: 'a/M',
        onChange: updateUniformsLive,
        help: 'Dimensionless spin. Positive = prograde disk, negative = retrograde.'
    });
    addControl(spinFolder, p.black_hole, 'spin_strength', {
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
    addControl(lookFolder, p.look, 'disk_gain', {
        min: 0.4,
        max: 4.0,
        step: 0.01,
        name: 'disk intensity',
        onChange: updateUniformsLive,
        help: 'Brightness multiplier applied to the accretion disk emission.'
    });
    addControl(lookFolder, p.look, 'glow', {
        min: 0.0,
        max: 2.0,
        step: 0.01,
        name: 'inner glow',
        onChange: updateUniformsLive,
        help: 'Extra bloom-like emphasis near the hotter inner disk region.'
    });
    addControl(lookFolder, p.look, 'doppler_boost', {
        min: 0.0,
        max: 2.5,
        step: 0.01,
        name: 'doppler boost',
        onChange: updateUniformsLive,
        help: 'Controls visual strength of relativistic beaming contrast.'
    });
    addControl(lookFolder, p.look, 'aberration_strength', {
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
    folder.open();

    folder = gui.addFolder('Planet');
    addControl(folder, p.planet, 'enabled', {
        name: 'enabled',
        help: 'Adds an orbiting planet used as a reference object.',
        onChange: function(enabled) {
        updateShader();
        var controls = $('.indirect-planet-controls').show();
        if (enabled) controls.show();
        else controls.hide();
        }
    });
    addControl(folder, p.planet, 'distance', {
        min: 1.5,
        step: 0.1,
        name: 'distance',
        onChange: updateUniforms,
        help: 'Orbital radius of the planet.'
    });
    addControl(folder, p.planet, 'radius', {
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
        onChange: updateShader,
        help: 'Changes apparent incoming ray direction due to observer velocity.'
    });
    addControl(folder, p, 'beaming', {
        name: 'beaming (intensity)',
        onChange: updateShader,
        help: 'Applies relativistic intensity boosting/dimming.'
    });
    addControl(folder, p, 'physical_beaming', {
        name: 'physical (D³ Liouville)',
        onChange: updateShader,
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
    addControl(folder, p, 'lorentz_contraction', {
        name: 'lorentz contraction',
        onChange: updateShader,
        className: 'planet-controls indirect-planet-controls',
        help: 'Applies length contraction effects to moving scene elements.'
    });

    folder.open();

    folder = gui.addFolder('Time');
    addControl(folder, p, 'light_travel_time', {
        onChange: updateShader,
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

}

function onWindowResize( event ) {
    renderer.setSize( window.innerWidth, window.innerHeight );
    updateUniforms();
}

function initializeCamera(camera) {

    var pitchAngle = 3.0, yawAngle = 0.0;

    // there are nicely named methods such as "lookAt" in the camera object
    // but there do not do a thing to the projection matrix due to an internal
    // representation of the camera coordinates using a quaternion (nice)
    camera.matrixWorldInverse.makeRotationX(degToRad(-pitchAngle));
    camera.matrixWorldInverse.multiply(new THREE.Matrix4().makeRotationY(degToRad(-yawAngle)));

    var m = camera.matrixWorldInverse.elements;

    camera.position.set(m[2], m[6], m[10]);
}

function updateCamera( event ) {

    var zoom_dist = camera.position.length();
    var m = camera.matrixWorldInverse.elements;
    var camera_matrix;

    if (shader.parameters.observer.motion) {
        camera_matrix = new THREE.Matrix3();
    }
    else {
        camera_matrix = observer.orientation;
    }

    camera_matrix.set(
        // row-major, not the same as .elements (nice)
        // y and z swapped for a nicer coordinate system
        m[0], m[1], m[2],
        m[8], m[9], m[10],
        m[4], m[5], m[6]
    );

    if (shader.parameters.observer.motion) {

        observer.orientation = observer.orbitalFrame().multiply(camera_matrix);

    } else {

        var p = new THREE.Vector3(
            camera_matrix.elements[6],
            camera_matrix.elements[7],
            camera_matrix.elements[8]);

        var dist = shader.parameters.observer.distance;
        observer.position.set(-p.x*dist, -p.y*dist, -p.z*dist);
        observer.velocity.set(0,0,0);
    }
}

function frobeniusDistance(matrix1, matrix2) {
    var sum = 0.0;
    for (var i in matrix1.elements) {
        var diff = matrix1.elements[i] - matrix2.elements[i];
        sum += diff*diff;
    }
    return Math.sqrt(sum);
}

function animate() {
    requestAnimationFrame( animate );

    camera.updateMatrixWorld();
    camera.matrixWorldInverse.getInverse( camera.matrixWorld );

    if (shader.needsUpdate || shader.hasMovingParts() ||
        frobeniusDistance(camera.matrixWorldInverse, lastCameraMat) > 1e-10) {

        shader.needsUpdate = false;
        render();
        lastCameraMat = camera.matrixWorldInverse.clone();
    }
    stats.update();
}

var lastCameraMat = new THREE.Matrix4().identity();

var getFrameDuration = (function() {
    var lastTimestamp = new Date().getTime();
    return function() {
        var timestamp = new Date().getTime();
        var diff = (timestamp - lastTimestamp) / 1000.0;
        lastTimestamp = timestamp;
        return diff;
    };
})();

function render() {
    observer.move(getFrameDuration());
    if (shader.parameters.observer.motion) updateCamera();
    updateUniforms();
    renderer.render( scene, camera );
}
