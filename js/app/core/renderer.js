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
var cameraPan = new THREE.Vector2(0, 0);
var distanceController = null;
var refreshAllControllersGlobal = null; // Will be set in setupGUI
var bloomPass = null;

// ─── Freefall Dive State ──────────────────────────────────────────────────────
// Tracks the dive animation that plunges the observer through the event horizon
// into the black hole interior with physically accurate geodesic ray tracing.
var diveState = {
    active: false,
    paused: false,
    speed: 1.0,
    cinematic: false,   // auto-vary speed for maximum visual drama
    autoOrient: true,
    currentR: 11.0,
    direction: new THREE.Vector3(1, 0, 0),
    startPosition: new THREE.Vector3(10, 0, 0),
    startVelocity: new THREE.Vector3(0, 1, 0),
    prevMotionState: true,
    prevDistance: 11.0,
    reachedSingularity: false
};
// ─────────────────────────────────────────────────────────────────────────────
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
        cam_pan: { type: "v2", value: new THREE.Vector2() },

        interior_mode: { type: "f", value: 0.0 },

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
        uniforms.cam_pan.value.set(cameraPan.x, cameraPan.y);

        // Interior mode: enable when observer is inside the event horizon.
        // The Binet equation is valid at all r; interior_mode tells the shader
        // to trace past u = 1 and use the analytical escape classification.
        // Must be exactly at the horizon (r_s = 1), not a padded threshold,
        // because the escape classifier assumes u0 > 1.
        var obsR = observer.position.length();
        uniforms.interior_mode.value = (obsR < 1.0) ? 1.0 : 0.0;

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
    cameraControls.panCallback = function(deltaX, deltaY, width, height) {
        var panSpeed = 0.75;
        var maxPan = 0.45;
        cameraPan.x -= 2.0 * deltaX / width * panSpeed;
        cameraPan.y += 2.0 * deltaY / height * panSpeed;
        cameraPan.x = Math.max(-maxPan, Math.min(maxPan, cameraPan.x));
        cameraPan.y = Math.max(-maxPan, Math.min(maxPan, cameraPan.y));
        shader.needsUpdate = true;
        return true;
    };
    cameraControls.addEventListener( 'change', updateCamera );
    updateCamera();

    onWindowResize();

    window.addEventListener( 'resize', onWindowResize, false );

    setupGUI();
}

// ─── Dive Animation Functions ───────────────────────────────────────────────
// Physics: Free-fall from rest at infinity in Schwarzschild geometry.
// Proper-time equation of motion: dr/dτ = -sqrt(r_s/r) = -sqrt(1/r)
// in units where r_s = 1.  Observer velocity (locally measured three-velocity):
// v = sqrt(1/r), capped at 0.998 c to avoid numerical divergence at horizon.
//
// The Binet equation d²u/dφ² = -u + (3/2)u² is valid for all r including
// inside the horizon.  Rays traced backward from an interior observer:
//   - Rays with impact parameter b < b_crit escape outward through the horizon
//     and show the external universe (background stars, accretion disk).
//   - Rays with b >= b_crit fall to the singularity → rendered black.
// This naturally produces the shrinking "window to the universe" effect as
// the observer approaches the singularity.
// ─────────────────────────────────────────────────────────────────────────────

function startDive() {
    if (diveState.active && !diveState.paused) {
        diveState.paused = true;
        updateDiveUI();
        return;
    }
    if (diveState.paused) {
        diveState.paused = false;
        updateDiveUI();
        return;
    }

    // Save current observer state for reset
    diveState.prevMotionState = shader.parameters.observer.motion;
    diveState.prevDistance = shader.parameters.observer.distance;
    diveState.startPosition = observer.position.clone();
    diveState.startVelocity = observer.velocity.clone();

    // Disable orbital motion — dive controls the observer now
    shader.parameters.observer.motion = false;

    // Dive direction = radially inward from current position
    diveState.direction = observer.position.clone().normalize();
    diveState.currentR = observer.position.length();
    diveState.active = true;
    diveState.paused = false;
    diveState.reachedSingularity = false;

    // Boost ray steps for interior — need more integration steps to trace
    // rays that cross the horizon boundary twice.
    if (shader.parameters.n_steps < 400) {
        shader.parameters.n_steps = 400;
    }
    if (shader.parameters.max_revolutions < 3.0) {
        shader.parameters.max_revolutions = 3.0;
    }
    shader.parameters.rk4_integration = true;

    scene.updateShader();
    updateCamera();
    updateDiveUI();
    if (refreshAllControllersGlobal) refreshAllControllersGlobal();
}

function resetDive() {
    diveState.active = false;
    diveState.paused = false;
    diveState.reachedSingularity = false;

    // Restore pre-dive observer state
    shader.parameters.observer.motion = diveState.prevMotionState;
    shader.parameters.observer.distance = diveState.prevDistance;
    diveState.currentR = diveState.prevDistance;

    observer.position.copy(diveState.startPosition);
    observer.velocity.copy(diveState.startVelocity);

    scene.updateShader();
    updateCamera();
    shader.needsUpdate = true;
    updateDiveUI();
    updateDiveFade();
    if (refreshAllControllersGlobal) refreshAllControllersGlobal();
}

// Cinematic speed envelope: scales the fall speed so visually rich regions
// (photon sphere r≈1.5, event horizon r≈1.0) play out slowly while the
// uneventful far-field approach is fast-forwarded.
//   r > 3  : up to 3× faster than base speed
//   r ≈ 1.5: ~0.25× (photon-sphere lensing slowdown)
//   r ≈ 1.0: ~0.15× (horizon-crossing slowdown)
//   r < 0.8: ~0.5× (inside — watch the escape cone shrink)
function cinematicFactor(r) {
    var farBoost   = 2.0 * Math.max(r - 3.0, 0.0) / 7.0;  // speeds up distant approach
    var photonSlow = 3.0 * Math.exp(-Math.pow((r - 1.5) / 0.30, 2));
    var horizonSlow= 5.0 * Math.exp(-Math.pow((r - 1.0) / 0.22, 2));
    return (1.0 + farBoost) / (1.0 + photonSlow + horizonSlow);
}

function seekDive(targetR) {
    if (!diveState.active && !diveState.reachedSingularity) return;
    targetR = Math.max(0.08, Math.min(diveState.prevDistance, targetR));

    // Allow scrubbing back from singularity
    if (diveState.reachedSingularity && targetR > 0.12) {
        diveState.reachedSingularity = false;
        diveState.active = true;
    }

    diveState.currentR = targetR;
    if (targetR <= 0.09) {
        diveState.reachedSingularity = true;
    }
    diveState.paused = true;

    // Update observer position and velocity for the new radius
    observer.position.copy(diveState.direction.clone().multiplyScalar(targetR));
    var v = Math.min(Math.sqrt(1.0 / targetR), 0.998);
    observer.velocity.copy(diveState.direction.clone().multiplyScalar(-v));
    shader.parameters.observer.distance = targetR;

    shader.needsUpdate = true;
    updateCamera();
    updateDiveUI();
}

function updateDive(dt) {
    if (!diveState.active || diveState.paused || diveState.reachedSingularity) return;

    var r = diveState.currentR;
    if (r < 0.08) {
        diveState.reachedSingularity = true;
        updateDiveUI();
        return;
    }

    // Free-fall from rest at infinity: dr/dτ = -sqrt(r_s/r) = -sqrt(1/r)
    // RK2 (midpoint method) for stable integration
    var effectiveSpeed = diveState.cinematic
        ? diveState.speed * cinematicFactor(r)
        : diveState.speed;
    var fallDt = dt * effectiveSpeed * shader.parameters.time_scale;
    var k1 = -Math.sqrt(1.0 / r) * fallDt;
    var rMid = Math.max(r + k1 * 0.5, 0.01);
    var k2 = -Math.sqrt(1.0 / rMid) * fallDt;
    var newR = Math.max(r + k2, 0.08);

    diveState.currentR = newR;

    // Update observer position along dive direction
    observer.position.copy(diveState.direction.clone().multiplyScalar(newR));

    // Free-fall three-velocity (capped < c to avoid numerical singularity)
    var v = Math.min(Math.sqrt(1.0 / newR), 0.998);
    observer.velocity.copy(diveState.direction.clone().multiplyScalar(-v));

    // Sync shader distance parameter
    shader.parameters.observer.distance = newR;

    // Advance observer time (proper time of the free-falling observer)
    observer.time += dt * effectiveSpeed * shader.parameters.time_scale;

    // Trigger shader recompile when crossing the horizon (interior mode transition)
    if (newR < 1.0 && r >= 1.0) {
        scene.updateShader();
    }

    shader.needsUpdate = true;
    updateDiveUI();
}

function updateDiveUI() {
    var radiusEl = document.getElementById('dive-radius');
    var velocityEl = document.getElementById('dive-velocity');
    var statusEl = document.getElementById('dive-status');
    var btnEl = document.getElementById('dive-start-btn');
    var resetBtn = document.getElementById('dive-reset-btn');
    var horizonBar = document.getElementById('dive-horizon-bar');

    if (!radiusEl) return;

    var r = diveState.currentR;
    var v = r > 0.01 ? Math.min(Math.sqrt(1.0 / r), 0.999) : 0.999;

    radiusEl.innerHTML = 'r = ' + r.toFixed(3) + ' r<sub>s</sub>';
    velocityEl.textContent = 'v = ' + v.toFixed(3) + ' c';

    // Show effective speed (with cinematic multiplier if active)
    var speedEl = document.getElementById('dive-speed-val');
    if (speedEl && diveState.active) {
        var effSpd = diveState.cinematic
            ? diveState.speed * cinematicFactor(r) : diveState.speed;
        speedEl.textContent = (effSpd < 0.1
            ? effSpd.toFixed(2) : effSpd.toFixed(1)) + '×';
    }

    // Update horizon proximity bar
    if (horizonBar) {
        // Map r from 0..startR to bar progress (100% = at singularity)
        var progress = Math.max(0, Math.min(100, (1.0 - r / Math.max(diveState.prevDistance, 1)) * 100));
        horizonBar.style.width = progress + '%';
        if (r < 1.0) {
            horizonBar.className = 'dive-horizon-fill inside';
        } else {
            horizonBar.className = 'dive-horizon-fill outside';
        }
    }

    if (diveState.reachedSingularity) {
        statusEl.textContent = '\u26a0 Singularity reached';
        statusEl.className = 'dive-status singularity';
        btnEl.textContent = '\u25b6 START DIVE';
        btnEl.disabled = true;
    } else if (!diveState.active) {
        statusEl.textContent = 'Ready';
        statusEl.className = 'dive-status ready';
        btnEl.textContent = '\u25b6 START DIVE';
        btnEl.disabled = false;
    } else if (diveState.paused) {
        statusEl.textContent = '\u23f8 Paused at r = ' + r.toFixed(2);
        statusEl.className = 'dive-status paused';
        btnEl.textContent = '\u25b6 RESUME';
    } else if (r > 1.5) {
        statusEl.textContent = '\u2193 Approaching horizon';
        statusEl.className = 'dive-status approaching';
        btnEl.textContent = '\u23f8 PAUSE';
    } else if (r > 1.0) {
        statusEl.textContent = '\u26a1 Near event horizon!';
        statusEl.className = 'dive-status near-horizon';
        btnEl.textContent = '\u23f8 PAUSE';
    } else if (r > 0.3) {
        statusEl.textContent = '\u26a1 INSIDE event horizon';
        statusEl.className = 'dive-status inside';
        btnEl.textContent = '\u23f8 PAUSE';
    } else {
        statusEl.textContent = '\ud83c\udf00 Approaching singularity';
        statusEl.className = 'dive-status deep';
        btnEl.textContent = '\u23f8 PAUSE';
    }

    if (resetBtn) {
        resetBtn.disabled = !diveState.active && !diveState.reachedSingularity;
    }
    updateDiveFade();
}

function updateDiveFade() {
    // No artificial overlay — the physics naturally darkens the view
    // as the escape window shrinks toward the singularity.
}

function updateAxesGizmo() {
    var canvas = document.getElementById('axes-gizmo');
    if (!canvas) return;
    var ctx = canvas.getContext('2d');
    if (!ctx) return;
    var w = canvas.width, h = canvas.height;
    var cx = w * 0.5, cy = h * 0.5;
    var len = 26;

    ctx.clearRect(0, 0, w, h);

    var e = observer.orientation.elements;
    // Orientation cols: cam_x=col0, cam_y=col1, cam_z=col2
    // Project world axis V to screen: sx = dot(V, cam_x), sy = -dot(V, cam_y)
    // (canvas Y is downward hence negation)
    var axes = [
        { name: 'X', color: '#ff4444', sx: e[0], sy: -e[3], depth: e[6] },
        { name: 'Y', color: '#44ff44', sx: e[1], sy: -e[4], depth: e[7] },
        { name: 'Z', color: '#4488ff', sx: e[2], sy: -e[5], depth: e[8] }
    ];

    // Draw back-to-front (most-toward-camera drawn last)
    axes.sort(function(a, b) { return a.depth - b.depth; });

    for (var i = 0; i < axes.length; i++) {
        var ax = axes[i];
        var ex = cx + ax.sx * len;
        var ey = cy + ax.sy * len;
        ctx.globalAlpha = ax.depth > 0 ? 1.0 : 0.3;

        // Shaft
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(ex, ey);
        ctx.strokeStyle = ax.color;
        ctx.lineWidth = 2;
        ctx.stroke();

        // Arrow tip
        ctx.beginPath();
        ctx.arc(ex, ey, 3, 0, 2 * Math.PI);
        ctx.fillStyle = ax.color;
        ctx.fill();

        // Label
        ctx.font = 'bold 11px monospace';
        ctx.fillStyle = ax.color;
        ctx.fillText(ax.name, ex + 5, ey + 4);
        ctx.globalAlpha = 1.0;
    }
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
    var MAX_FRAME_DT = 0.1; // Cap at 100 ms to prevent time jumps after tab switch
    var lastTimestamp = new Date().getTime();

    // When the page becomes visible again after being hidden, reset the
    // timestamp so the first frame doesn't get a huge accumulated delta.
    document.addEventListener('visibilitychange', function() {
        if (!document.hidden) {
            lastTimestamp = new Date().getTime();
        }
    });

    return function() {
        var timestamp = new Date().getTime();
        var diff = (timestamp - lastTimestamp) / 1000.0;
        lastTimestamp = timestamp;
        return Math.min(diff, MAX_FRAME_DT);
    };
})();

function render() {
    var dt = getFrameDuration();

    if (diveState.active && !diveState.paused && !diveState.reachedSingularity) {
        updateDive(dt);
        updateCamera(); // Keep camera orientation synced
    } else {
        observer.move(dt);
        if (shader.parameters.observer.motion) updateCamera();
    }

    updateUniforms();

    if (shader.parameters.bloom.enabled && bloomPass) {
        bloomPass.render(renderer, scene, camera, shader.parameters.bloom);
    } else {
        renderer.render( scene, camera );
    }
    updateAxesGizmo();
}
