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
        var ang_vel = v / r;
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
        dt = Math.sqrt((dt*dt * (1.0 - v*v)) / (1-1.0/r));
    }

    this.time += dt;
};

var container, stats;
var camera, scene, renderer, cameraControls, shader = null;
var observer = new Observer();
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
        cinematic_tonemap: false,
        quality: 'high',
        kerr_mode: 'realtime_full_kerr_core',
        accretion_disk: true,
        disk_temperature: 8500.0,
        black_hole: {
            spin_enabled: true,
            spin: 0.60,
            spin_strength: 0.55
        },
        look: {
            exposure: 0.96,
            disk_gain: 1.75,
            glow: 0.42,
            doppler_boost: 1.0,
            aberration_strength: 1.4,
            star_gain: 0.55,
            galaxy_gain: 0.55
        },
        planet: {
            enabled: true,
            distance: 7.0,
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
        that.parameters.kerr_full_core = (that.parameters.kerr_mode === 'realtime_full_kerr_core' ||
            that.parameters.kerr_mode === 'offline_accurate');
        that.parameters.kerr_offline = (that.parameters.kerr_mode === 'offline_accurate');
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
        disk_temperature: { type: "f", value: 8500.0 },
        accretion_inner_r: { type: "f", value: 3.0 },
        bh_spin: { type: "f", value: 0.60 },
        bh_spin_strength: { type: "f", value: 0.55 },
        bh_rotation_enabled: { type: "f", value: 1.0 },
        look_exposure: { type: "f", value: 0.96 },
        look_disk_gain: { type: "f", value: 1.75 },
        look_glow: { type: "f", value: 0.42 },
        look_doppler_boost: { type: "f", value: 1.0 },
        look_aberration_strength: { type: "f", value: 1.4 },
        look_star_gain: { type: "f", value: 0.55 },
        look_galaxy_gain: { type: "f", value: 0.55 },

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

    var gui = new dat.GUI();

    function applyKerrMode(mode) {
        // Preserve selected quality in fast mode, but enforce physically heavier settings
        // for Kerr geodesic integration modes.
        if (mode === 'realtime_full_kerr_core') {
            p.rk4_integration = true;
            p.max_revolutions = Math.max(p.max_revolutions, 3.0);
            p.n_steps = Math.max(p.n_steps, 520);
            p.sample_count = Math.max(p.sample_count, 4);
            p.cinematic_tonemap = true;
        } else if (mode === 'offline_accurate') {
            p.rk4_integration = true;
            p.max_revolutions = 4.0;
            p.n_steps = 960;
            p.sample_count = 8;
            p.cinematic_tonemap = true;
        }

        updateShader();
    }

    function applyQualityPreset(value) {
        $('.planet-controls').show();

        switch(value) {
        case 'fast':
            p.n_steps = 40;
            p.sample_count = 1;
            p.max_revolutions = 1.5;
            p.rk4_integration = false;
            p.cinematic_tonemap = true;
            $('.planet-controls').hide();
            break;
        case 'medium':
            p.n_steps = 100;
            p.sample_count = 1;
            p.max_revolutions = 2.0;
            p.rk4_integration = false;
            p.cinematic_tonemap = true;
            break;
        case 'high':
            p.n_steps = 320;
            p.sample_count = 4;
            p.max_revolutions = 3.2;
            p.rk4_integration = true;
            p.cinematic_tonemap = true;
            break;
        }

        applyKerrMode(p.kerr_mode);
    }

    gui.add(p, 'quality', ['fast', 'medium', 'high']).onChange(applyQualityPreset);
    gui.add(p, 'kerr_mode', [
        'fast',
        'realtime_full_kerr_core',
        'offline_accurate'
    ]).name('kerr solver mode').onChange(applyKerrMode);
    applyQualityPreset(p.quality);
    applyKerrMode(p.kerr_mode);
    var diskFolder = gui.addFolder('Accretion disk');
    diskFolder.add(p, 'accretion_disk').onChange(updateShader);
    diskFolder.add(p, 'disk_temperature')
        .min(DISK_TEMPERATURE_MIN)
        .max(DISK_TEMPERATURE_MAX)
        .step(1)
        .name('temperature (K)')
        .onChange(updateUniformsLive);
    diskFolder.open();

    var spinFolder = gui.addFolder('Black hole');
    spinFolder.add(p.black_hole, 'spin_enabled')
        .name('rotating shadow')
        .onChange(updateUniformsLive);
    spinFolder.add(p.black_hole, 'spin')
        .min(-0.99)
        .max(0.99)
        .step(0.01)
        .name('a/M')
        .onChange(updateUniformsLive);
    spinFolder.add(p.black_hole, 'spin_strength')
        .min(0.0)
        .max(1.4)
        .step(0.01)
        .name('shadow squeeze')
        .onChange(updateUniformsLive);
    spinFolder.open();

    var lookFolder = gui.addFolder('Look');
    lookFolder.add(p.look, 'exposure')
        .min(0.6)
        .max(2.5)
        .step(0.01)
        .name('exposure')
        .onChange(updateUniformsLive);
    lookFolder.add(p.look, 'disk_gain')
        .min(0.4)
        .max(4.0)
        .step(0.01)
        .name('disk intensity')
        .onChange(updateUniformsLive);
    lookFolder.add(p.look, 'glow')
        .min(0.0)
        .max(2.0)
        .step(0.01)
        .name('inner glow')
        .onChange(updateUniformsLive);
    lookFolder.add(p.look, 'doppler_boost')
        .min(0.0)
        .max(2.5)
        .step(0.01)
        .name('doppler boost')
        .onChange(updateUniformsLive);
    lookFolder.add(p.look, 'aberration_strength')
        .min(0.0)
        .max(3.0)
        .step(0.01)
        .name('aberration strength')
        .onChange(updateUniformsLive);
    lookFolder.add(p.look, 'star_gain')
        .min(0.0)
        .max(2.5)
        .step(0.01)
        .name('star gain')
        .onChange(updateUniformsLive);
    lookFolder.add(p.look, 'galaxy_gain')
        .min(0.0)
        .max(2.5)
        .step(0.01)
        .name('galaxy gain')
        .onChange(updateUniformsLive);
    lookFolder.open();

    var folder = gui.addFolder('Observer');
    folder.add(p.observer, 'motion').onChange(function(motion) {
        updateCamera();
        updateShader();
        if (motion) {
            hint.text('Moving observer; drag to rotate camera');
        } else {
            hint.text('Stationary observer; drag to orbit around');
        }
        hint.fadeIn();
    });
    folder.add(p.observer, 'distance').min(1.5).max(30).onChange(updateCamera);
    folder.open();

    folder = gui.addFolder('Planet');
    folder.add(p.planet, 'enabled').onChange(function(enabled) {
        updateShader();
        var controls = $('.indirect-planet-controls').show();
        if (enabled) controls.show();
        else controls.hide();
    });
    folder.add(p.planet, 'distance').min(1.5).onChange(updateUniforms);
    folder.add(p.planet, 'radius').min(0.01).max(2.0).onChange(updateUniforms);
    $(folder.domElement).addClass('planet-controls');
    //folder.open();

    function setGuiRowClass(guiEl, klass) {
        $(guiEl.domElement).parent().parent().addClass(klass);
    }

    folder = gui.addFolder('Relativistic effects');
    folder.add(p, 'aberration').name('aberration (ray dir)').onChange(updateShader);
    folder.add(p, 'beaming').name('beaming (intensity)').onChange(updateShader);
    folder.add(p, 'physical_beaming').name('physical (D³ Liouville)').onChange(updateShader);
    folder.add(p, 'doppler_shift').name('doppler shift (color)').onChange(updateShader);
    setGuiRowClass(
        folder.add(p, 'gravitational_time_dilation').onChange(updateShader),
        'planet-controls indirect-planet-controls');
    setGuiRowClass(
        folder.add(p, 'lorentz_contraction').onChange(updateShader),
        'planet-controls indirect-planet-controls');

    folder.open();

    folder = gui.addFolder('Time');
    folder.add(p, 'light_travel_time').onChange(updateShader);
    folder.add(p, 'time_scale').min(0);
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
