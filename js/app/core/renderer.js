"use strict";

import { THREE, Stats, Detector, $ } from '../vendor.js';
import { Shader } from './shader.js';
import {
    formatThousands,
    clampObserverDistance,
    clampPlanetOrbitDistance
} from './observer.js';
import { initializeCamera, updateCamera } from '../scene/camera.js';
import { setupBloom } from '../graphics/bloom.js';
import { setupTemporalAA } from './temporal-aa.js';
import { setupGUI } from '../ui/gui.js';
import {
    applyQualityPresetValues
} from '../ui/quality-presets.js';
import {
    beginQualityBenchmarkIfNeeded,
    isLikelyMobileDevice,
    readStoredQualityPreset,
    clampResolutionScale
} from './renderer-quality.js';
import { updateObserverDistanceBinding } from './ui-bindings.js';
import {
    DISK_TEMPERATURE_MIN,
    DISK_TEMPERATURE_MAX,
    container,
    stats,
    camera,
    scene,
    renderer,
    cameraControls,
    shader,
    observer,
    cameraPan,
    bloomPass,
    taaPass,
    shaderUniforms,
    baseDevicePixelRatio,
    lastTaaCameraMat,
    updateUniforms,
    setContainer,
    setStats,
    setCamera,
    setScene,
    setRenderer,
    setCameraControls,
    setShader,
    setBloomPass,
    setTaaPass,
    setShaderUniforms,
    setIsMobileClient,
    setRendererContextLost,
    setUpdateUniforms,
    setApplyRenderScaleFromSettings,
    setResetTemporalAAHistory
} from './runtime-state.js';

if (!Detector.webgl) Detector.addGetWebGLMessage();

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

function resetTemporalAAHistoryImpl() {
    if (taaPass) taaPass.reset();
    if (shaderUniforms && shaderUniforms.taa_jitter) {
        shaderUniforms.taa_jitter.value.set(0, 0);
    }
    lastTaaCameraMat.identity();
}

export function resizeRendererAndPasses() {
    if (!renderer || !shader) return;

    var scale = clampResolutionScale(shader.parameters.resolution_scale);
    shader.parameters.resolution_scale = scale;

    renderer.setPixelRatio(baseDevicePixelRatio * scale);
    renderer.setSize(window.innerWidth, window.innerHeight);

    var w = renderer.domElement.width;
    var h = renderer.domElement.height;

    if (bloomPass) bloomPass.resize(w, h);
    if (taaPass) taaPass.resize(w, h);
}

function applyRenderScaleFromSettingsImpl() {
    resizeRendererAndPasses();
    resetTemporalAAHistoryImpl();
    if (updateUniforms) updateUniforms();
    if (shader) shader.needsUpdate = true;
}

export function onWindowResize() {
    resizeRendererAndPasses();
    resetTemporalAAHistoryImpl();
    if (updateUniforms) updateUniforms();
}

setApplyRenderScaleFromSettings(applyRenderScaleFromSettingsImpl);
setResetTemporalAAHistory(resetTemporalAAHistoryImpl);

export function init(glslSource, textures) {
    var mobileClient = isLikelyMobileDevice();
    setShader(new Shader(glslSource));
    setIsMobileClient(mobileClient);
    if (mobileClient) {
        document.body.classList.add('mobile-ui');
    }

    var storedQualityPreset = readStoredQualityPreset();
    var initialQualityPreset = storedQualityPreset || 'optimal';
    if (typeof applyQualityPresetValues === 'function') {
        applyQualityPresetValues(shader.parameters, initialQualityPreset);
    } else {
        shader.parameters.quality = initialQualityPreset;
    }

    var nextContainer = document.createElement('div');
    document.body.appendChild(nextContainer);
    setContainer(nextContainer);

    setScene(new THREE.Scene());

    var geometry = new THREE.PlaneBufferGeometry(2, 2);
    var uniforms = {
        time: { type: 'f', value: 0 },
        turbulence_time_offset: { type: 'f', value: 0.0 },
        turbulence_loop_enabled: { type: 'f', value: 0.0 },
        turbulence_loop_seconds: { type: 'f', value: 20.0 },
        resolution: { type: 'v2', value: new THREE.Vector2() },
        cam_pos: { type: 'v3', value: new THREE.Vector3() },
        cam_x: { type: 'v3', value: new THREE.Vector3() },
        cam_y: { type: 'v3', value: new THREE.Vector3() },
        cam_z: { type: 'v3', value: new THREE.Vector3() },
        cam_vel: { type: 'v3', value: new THREE.Vector3() },
        cam_pan: { type: 'v2', value: new THREE.Vector2() },
        taa_jitter: { type: 'v2', value: new THREE.Vector2() },
        interior_mode: { type: 'f', value: 0.0 },
        planet_distance: { type: 'f' },
        planet_radius: { type: 'f' },
        disk_temperature: { type: 'f', value: 10000.0 },
        accretion_inner_r: { type: 'f', value: 3.0 },
        bh_spin: { type: 'f', value: 0.90 },
        bh_spin_strength: { type: 'f', value: 1.0 },
        bh_rotation_enabled: { type: 'f', value: 1.0 },
        photon_spin_lensing_scale: { type: 'f', value: 1.0 },
        look_exposure: { type: 'f', value: 1.0 },
        look_disk_gain: { type: 'f', value: 1.0 },
        look_glow: { type: 'f', value: 0.0 },
        look_doppler_boost: { type: 'f', value: 1.0 },
        look_aberration_strength: { type: 'f', value: 1.0 },
        look_star_gain: { type: 'f', value: 1.0 },
        look_galaxy_gain: { type: 'f', value: 1.0 },
        look_tonemap_mode: { type: 'f', value: 0.0 },
        torus_r0: { type: 'f', value: 4.0 },
        torus_h_ratio: { type: 'f', value: 0.45 },
        torus_radial_falloff: { type: 'f', value: 2.5 },
        torus_opacity: { type: 'f', value: 0.015 },
        torus_outer_radius: { type: 'f', value: 3.5 },
        slim_h_ratio: { type: 'f', value: 0.15 },
        slim_opacity: { type: 'f', value: 0.6 },
        slim_puff_factor: { type: 'f', value: 2.5 },
        jet_half_angle: { type: 'f', value: 5.0 },
        jet_lorentz: { type: 'f', value: 3.0 },
        jet_brightness: { type: 'f', value: 1.2 },
        jet_length: { type: 'f', value: 30.0 },
        jet_magnetization: { type: 'f', value: 10.0 },
        jet_knot_spacing: { type: 'f', value: 6.0 },
        jet_corona_brightness: { type: 'f', value: 1.5 },
        jet_base_width: { type: 'f', value: 0.4 },
        jet_corona_extent: { type: 'f', value: 0.5 },
        grmhd_r_high: { type: 'f', value: 40.0 },
        grmhd_magnetic_beta: { type: 'f', value: 10.0 },
        grmhd_mad_flux: { type: 'f', value: 0.0 },
        grmhd_density_scale: { type: 'f', value: 1.0 },
        grmhd_turbulence_amp: { type: 'f', value: 1.0 },
        grmhd_electron_kappa: { type: 'f', value: 5.0 },
        grmhd_magnetic_field_str: { type: 'f', value: 1.0 },
        grav_blueshift_factor: { type: 'f', value: 1.0 },
        star_texture: { type: 't', value: textures.stars },
        galaxy_texture: { type: 't', value: textures.galaxy },
        planet_texture: { type: 't', value: textures.moon },
        spectrum_texture: { type: 't', value: textures.spectra }
    };
    setShaderUniforms(uniforms);

    function calculateISCO(chi, isPrograde) {
        var chi2 = chi * chi;
        var cbrt1MinusChi2 = Math.pow(Math.max(1 - chi2, 0), 1 / 3);
        var cbrt1PlusChi = Math.pow(1 + Math.abs(chi), 1 / 3);
        var cbrt1MinusChi = Math.pow(Math.max(1 - Math.abs(chi), 0), 1 / 3);
        var z1 = 1 + cbrt1MinusChi2 * (cbrt1PlusChi + cbrt1MinusChi);
        var z2 = Math.sqrt(3 * chi2 + z1 * z1);
        var sign = isPrograde ? -1 : 1;
        var iscoRg = 3 + z2 + sign * Math.sqrt((3 - z1) * (3 + z1 + 2 * z2));
        return iscoRg * 0.5;
    }

    function updateUniformsImpl() {
        shader.parameters.disk_temperature = Math.max(
            DISK_TEMPERATURE_MIN,
            Math.min(DISK_TEMPERATURE_MAX, shader.parameters.disk_temperature)
        );
        shader.parameters.planet.distance =
            clampPlanetOrbitDistance(shader.parameters.planet.distance);
        uniforms.planet_distance.value = shader.parameters.planet.distance;
        uniforms.planet_radius.value = shader.parameters.planet.radius;
        uniforms.disk_temperature.value = shader.parameters.disk_temperature;
        shader.parameters.observer.distance = clampObserverDistance(
            shader.parameters.observer.distance,
            shader.parameters.observer.motion
        );

        var spin = shader.parameters.black_hole.spin;
        var spinEnabled = shader.parameters.black_hole.spin_enabled;
        var spinMagnitude = Math.abs(spin);
        uniforms.accretion_inner_r.value = spinEnabled
            ? calculateISCO(spinMagnitude, true)
            : 3.0;

        uniforms.bh_spin.value = shader.parameters.black_hole.spin;
        uniforms.bh_spin_strength.value = shader.parameters.black_hole.spin_strength;
        uniforms.bh_rotation_enabled.value = spinEnabled ? 1.0 : 0.0;
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

        uniforms.grmhd_r_high.value = shader.parameters.grmhd.r_high;
        uniforms.grmhd_magnetic_beta.value = shader.parameters.grmhd.magnetic_beta;
        uniforms.grmhd_mad_flux.value = shader.parameters.grmhd.mad_flux;
        uniforms.grmhd_density_scale.value = shader.parameters.grmhd.density_scale;
        uniforms.grmhd_turbulence_amp.value = shader.parameters.grmhd.turbulence_amp;
        uniforms.grmhd_electron_kappa.value = shader.parameters.grmhd.electron_kappa;
        uniforms.grmhd_magnetic_field_str.value = shader.parameters.grmhd.magnetic_field_str;
        uniforms.turbulence_loop_enabled.value =
            shader.parameters.turbulence_loop_enabled ? 1.0 : 0.0;
        uniforms.turbulence_loop_seconds.value = Math.max(
            1e-4,
            parseFloat(shader.parameters.turbulence_loop_seconds) || 20.0
        );

        uniforms.resolution.value.x = renderer.domElement.width;
        uniforms.resolution.value.y = renderer.domElement.height;

        uniforms.time.value = observer.time;
        uniforms.turbulence_time_offset.value = observer.turbulenceTimeOffset || 0.0;

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

        var obsR = observer.position.length();
        uniforms.interior_mode.value = (obsR < 1.0) ? 1.0 : 0.0;

        var observerSpeedSq = observer.velocity.lengthSq();
        var photonSpinLensingScale = 1.0;
        if (observerSpeedSq < 1e-8) {
            var fade = (obsR - 1.5) / 0.2;
            fade = Math.max(0.0, Math.min(1.0, fade));
            photonSpinLensingScale = fade * fade * (3.0 - 2.0 * fade);
        }
        uniforms.photon_spin_lensing_scale.value = photonSpinLensingScale;

        if (obsR > 1.0) {
            uniforms.grav_blueshift_factor.value =
                Math.sqrt(Math.max(1.0 - 1.0 / obsR, 0.001));
        } else {
            uniforms.grav_blueshift_factor.value = 1.0;
        }

        updateEffectLabels();
    }

    setUpdateUniforms(updateUniformsImpl);

    var material = new THREE.ShaderMaterial({
        uniforms: uniforms,
        vertexShader: $('#vertex-shader').text()
    });

    scene.updateShader = function() {
        material.fragmentShader = shader.compile();
        material.needsUpdate = true;
        shader.needsUpdate = true;
        resetTemporalAAHistoryImpl();
    };

    scene.updateShader();
    scene.add(new THREE.Mesh(geometry, material));

    var nextRenderer = new THREE.WebGLRenderer({
        antialias: true,
        powerPreference: 'high-performance'
    });
    nextRenderer.domElement.style.touchAction = 'none';
    container.appendChild(nextRenderer.domElement);
    setRenderer(nextRenderer);

    renderer.domElement.addEventListener('webglcontextlost', function(e) {
        e.preventDefault();
        setRendererContextLost(true);
        console.warn('WebGL context lost - GPU driver may have reset (TDR).');
    }, false);
    renderer.domElement.addEventListener('webglcontextrestored', function() {
        setRendererContextLost(false);
        console.info('WebGL context restored.');
        if (shader) shader.needsUpdate = true;
        resetTemporalAAHistoryImpl();
    }, false);

    setBloomPass(setupBloom());
    setTaaPass(setupTemporalAA());

    var nextStats = new Stats();
    nextStats.domElement.style.position = 'fixed';
    nextStats.domElement.style.top = '0px';
    nextStats.domElement.style.left = '0px';
    nextStats.domElement.style.zIndex = '1000';
    container.appendChild(nextStats.domElement);
    setStats(nextStats);
    $(stats.domElement).addClass('hidden-phone');

    effectLabels.spin = document.getElementById('spin-label');
    effectLabels.temperature = document.getElementById('temperature-label');
    updateEffectLabels();

    var nextCamera = new THREE.PerspectiveCamera(
        45,
        window.innerWidth / window.innerHeight,
        1,
        80000
    );
    setCamera(nextCamera);
    initializeCamera(camera);

    var controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.target.set(0, 0, 0);
    controls.enableZoom = false;
    controls.panCallback = function(deltaX, deltaY, width, height) {
        var panSpeed = 0.75;
        var maxPan = 0.45;
        cameraPan.x -= 2.0 * deltaX / width * panSpeed;
        cameraPan.y += 2.0 * deltaY / height * panSpeed;
        cameraPan.x = Math.max(-maxPan, Math.min(maxPan, cameraPan.x));
        cameraPan.y = Math.max(-maxPan, Math.min(maxPan, cameraPan.y));
        shader.needsUpdate = true;
        return true;
    };
    controls.zoomCallback = function(dollyDeltaY) {
        var zoomBase = 1.08;
        var zoomFactor = dollyDeltaY > 0 ? (1.0 / zoomBase) : zoomBase;
        var newDist = shader.parameters.observer.distance * zoomFactor;
        newDist = clampObserverDistance(newDist, shader.parameters.observer.motion);
        shader.parameters.observer.distance = newDist;
        updateCamera();
        shader.needsUpdate = true;
        updateObserverDistanceBinding();
        return true;
    };
    controls.addEventListener('change', updateCamera);
    setCameraControls(controls);
    updateCamera();

    applyRenderScaleFromSettingsImpl();
    window.addEventListener('resize', onWindowResize, false);

    setupGUI();
    beginQualityBenchmarkIfNeeded();
}

export function updateAxesGizmo() {
    var canvas = document.getElementById('axes-gizmo');
    if (!canvas) return;
    var ctx = canvas.getContext('2d');
    if (!ctx) return;
    var w = canvas.width;
    var h = canvas.height;
    var cx = w * 0.5;
    var cy = h * 0.5;
    var len = 26;

    ctx.clearRect(0, 0, w, h);

    var e = observer.orientation.elements;
    var axes = [
        { name: 'X', color: '#ff4444', sx: e[0], sy: -e[3], depth: e[6] },
        { name: 'Y', color: '#44ff44', sx: e[1], sy: -e[4], depth: e[7] },
        { name: 'Z', color: '#4488ff', sx: e[2], sy: -e[5], depth: e[8] }
    ];

    axes.sort(function(a, b) { return a.depth - b.depth; });

    for (var i = 0; i < axes.length; i++) {
        var ax = axes[i];
        var ex = cx + ax.sx * len;
        var ey = cy + ax.sy * len;
        ctx.globalAlpha = ax.depth > 0 ? 1.0 : 0.3;

        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(ex, ey);
        ctx.strokeStyle = ax.color;
        ctx.lineWidth = 2;
        ctx.stroke();

        ctx.beginPath();
        ctx.arc(ex, ey, 3, 0, 2 * Math.PI);
        ctx.fillStyle = ax.color;
        ctx.fill();

        ctx.font = 'bold 11px monospace';
        ctx.fillStyle = ax.color;
        ctx.fillText(ax.name, ex + 5, ey + 4);
        ctx.globalAlpha = 1.0;
    }
}
