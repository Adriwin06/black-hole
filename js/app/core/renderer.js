// Role: Core renderer — declares all shared globals, builds the Three.js scene,
//       wires uniforms, initialises bloom, camera and GUI, then drives the
//       animate/render loop. init() is called by bootstrap.js once all GLSL
//       shards and textures have been fetched and are ready.

"use strict";
/*global THREE, Mustache, Stats, Detector, $, dat:false */
/*global document, window, setTimeout, requestAnimationFrame:false */

if ( ! Detector.webgl ) Detector.addGetWebGLMessage();

var DISK_TEMPERATURE_MIN = 4500.0;
var DISK_TEMPERATURE_MAX = 30000.0;

var container, stats;
var camera, scene, renderer, cameraControls, shader = null;
var observer = new Observer();
var distanceController = null;
var refreshAllControllersGlobal = null; // Will be set in setupGUI
var bloomPass = null;
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

var updateUniforms;

function init(glslSource, textures) {

    shader = new Shader(glslSource);

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
        torus_radial_falloff: { type: "f", value: 2.5 },
        torus_opacity: { type: "f", value: 0.015 },
        torus_outer_radius: { type: "f", value: 3.5 },

        slim_h_ratio: { type: "f", value: 0.15 },
        slim_opacity: { type: "f", value: 0.6 },
        slim_puff_factor: { type: "f", value: 2.5 },

        jet_half_angle: { type: "f", value: 5.0 },
        jet_lorentz: { type: "f", value: 3.0 },
        jet_brightness: { type: "f", value: 1.2 },
        jet_length: { type: "f", value: 30.0 },
        jet_magnetization: { type: "f", value: 10.0 },
        jet_knot_spacing: { type: "f", value: 6.0 },
        jet_corona_brightness: { type: "f", value: 1.5 },
        jet_base_width: { type: "f", value: 0.4 },
        jet_corona_extent: { type: "f", value: 0.5 },

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
        uniforms.torus_radial_falloff.value = shader.parameters.torus.radial_falloff;
        uniforms.torus_opacity.value = shader.parameters.torus.opacity;
        uniforms.torus_outer_radius.value = shader.parameters.torus.outer_radius;

        uniforms.slim_h_ratio.value = shader.parameters.slim.h_ratio;
        uniforms.slim_opacity.value = shader.parameters.slim.opacity;
        uniforms.slim_puff_factor.value = shader.parameters.slim.puff_factor;

        uniforms.jet_half_angle.value = shader.parameters.jet.half_angle;
        uniforms.jet_lorentz.value = shader.parameters.jet.lorentz_factor;
        uniforms.jet_brightness.value = shader.parameters.jet.brightness;
        uniforms.jet_length.value = shader.parameters.jet.length;
        uniforms.jet_magnetization.value = shader.parameters.jet.magnetization;
        uniforms.jet_knot_spacing.value = shader.parameters.jet.knot_spacing;
        uniforms.jet_corona_brightness.value = shader.parameters.jet.corona_brightness;
        uniforms.jet_base_width.value = shader.parameters.jet.base_width;
        uniforms.jet_corona_extent.value = shader.parameters.jet.corona_extent;

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

    // ============== BLOOM POST-PROCESSING ==============
    bloomPass = setupBloom();
    // ============== END BLOOM ==============

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

function onWindowResize( event ) {
    renderer.setSize( window.innerWidth, window.innerHeight );
    if (bloomPass) {
        bloomPass.resize(renderer.domElement.width, renderer.domElement.height);
    }
    updateUniforms();
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

    if (shader.parameters.bloom.enabled && bloomPass) {
        bloomPass.render(renderer, scene, camera, shader.parameters.bloom);
    } else {
        renderer.render( scene, camera );
    }
}
