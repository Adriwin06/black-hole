// Role: Core renderer — declares all shared globals, builds the Three.js scene,
//       wires uniforms, initialises bloom, camera and GUI, then drives the
//       animate/render loop. init() is called by bootstrap.js once all GLSL
//       shards and textures have been fetched and are ready.

"use strict";
/*global THREE, Mustache, Stats, Detector, $, dat:false, QUALITY_PRESETS, applyQualityPresetValues */
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
var taaPass = null;
var shaderUniforms = null;
var baseDevicePixelRatio = Math.max(window.devicePixelRatio || 1.0, 1.0);
var isMobileClient = false;
var lastTaaCameraMat = new THREE.Matrix4().identity();

var applyRenderScaleFromSettings = function() {};
var resetTemporalAAHistory = function() {};
var QUALITY_BENCHMARK_STORAGE_KEY = 'black-hole-quality-benchmark-v3';
var QUALITY_BENCHMARK_SCHEMA_VERSION = 3;
var QUALITY_BENCHMARK_TARGET_FRAME_MS = 32.0;
var QUALITY_BENCHMARK_HIGH_PROBE_GATE_MS = 18.5;
var QUALITY_BENCHMARK_HIGH_PROBE_GATE_MS_ULTRA_CAPABLE = 24.0;
var QUALITY_BENCHMARK_HIGH_TARGET_FRAME_MS_ULTRA_CAPABLE = 40.0;
var qualityBenchmarkState = null;

function getGpuRendererName() {
    if (!renderer || typeof renderer.getContext !== 'function') return '';
    try {
        var gl = renderer.getContext();
        if (!gl) return '';

        var ext = gl.getExtension('WEBGL_debug_renderer_info');
        if (ext && ext.UNMASKED_RENDERER_WEBGL) {
            var unmasked = gl.getParameter(ext.UNMASKED_RENDERER_WEBGL);
            if (typeof unmasked === 'string' && unmasked.length > 0) return unmasked;
        }

        var masked = gl.getParameter(gl.RENDERER);
        return (typeof masked === 'string') ? masked : '';
    } catch (err) {
        return '';
    }
}

function isUltraCapableGpu(rendererName) {
    if (!rendererName) return false;
    var gpu = rendererName.toLowerCase();
    if (gpu.indexOf('swiftshader') !== -1) return false;
    if (gpu.indexOf('nvidia') === -1 || gpu.indexOf('rtx') === -1) return false;

    if (/rtx\s*4070\s*super/.test(gpu)) return true;
    var match = gpu.match(/rtx\s*(\d{4})/);
    if (!match) return false;

    var model = parseInt(match[1], 10);
    return isFinite(model) && model >= 4080;
}

function readStoredQualityPreset() {
    var raw = null;
    try {
        raw = window.localStorage.getItem(QUALITY_BENCHMARK_STORAGE_KEY);
    } catch (err) {
        return null;
    }
    if (!raw) return null;

    try {
        var parsed = JSON.parse(raw);
        if (!parsed || parsed.version !== QUALITY_BENCHMARK_SCHEMA_VERSION) return null;
        var storedQuality = parsed.quality;
        if (storedQuality === 'fast') storedQuality = 'mobile';
        if (!storedQuality || !QUALITY_PRESETS[storedQuality]) return null;
        return storedQuality;
    } catch (err) {
        return null;
    }
}

function storeQualityPreset(qualityName, avgFrameMs) {
    if (!qualityName || !QUALITY_PRESETS[qualityName]) return;

    var roundedMs = null;
    if (typeof avgFrameMs === 'number' && isFinite(avgFrameMs)) {
        roundedMs = Math.round(avgFrameMs * 100) / 100;
    }

    var payload = {
        version: QUALITY_BENCHMARK_SCHEMA_VERSION,
        quality: qualityName,
        avg_frame_ms: roundedMs,
        timestamp: new Date().toISOString()
    };

    try {
        window.localStorage.setItem(QUALITY_BENCHMARK_STORAGE_KEY, JSON.stringify(payload));
    } catch (err) {
        // localStorage may be blocked; in that case auto-benchmark simply reruns next visit.
    }
}

function applyQualityPresetRuntime(qualityName) {
    if (!shader || typeof applyQualityPresetValues !== 'function') return false;

    var preset = applyQualityPresetValues(shader.parameters, qualityName);
    if (!preset) return false;

    if (preset.hide_planet_controls) {
        $('.planet-controls').hide();
    } else {
        $('.planet-controls').show();
    }

    applyRenderScaleFromSettings();
    if (scene && typeof scene.updateShader === 'function') {
        scene.updateShader();
    }
    if (refreshAllControllersGlobal) {
        refreshAllControllersGlobal();
    }
    return true;
}

function resetQualityBenchmarkCounters(state) {
    state.frameCount = 0;
    state.sampleCount = 0;
    state.accumulatedDt = 0.0;
}

function finishQualityBenchmark(qualityName, avgFrameMs) {
    applyQualityPresetRuntime(qualityName);
    storeQualityPreset(qualityName, avgFrameMs);
    qualityBenchmarkState = null;
}

function beginQualityBenchmarkIfNeeded() {
    if (qualityBenchmarkState) return;
    if (readStoredQualityPreset()) return;

    qualityBenchmarkState = {
        phase: 'optimal',
        warmupFrames: 24,
        sampleFrames: 72,
        frameCount: 0,
        sampleCount: 0,
        accumulatedDt: 0.0,
        optimalAvgMs: null,
        highProbeGateMs: QUALITY_BENCHMARK_HIGH_PROBE_GATE_MS,
        highTargetFrameMs: QUALITY_BENCHMARK_TARGET_FRAME_MS
    };
    var gpuRendererName = getGpuRendererName();
    if (isUltraCapableGpu(gpuRendererName)) {
        qualityBenchmarkState.highProbeGateMs = QUALITY_BENCHMARK_HIGH_PROBE_GATE_MS_ULTRA_CAPABLE;
        qualityBenchmarkState.highTargetFrameMs = QUALITY_BENCHMARK_HIGH_TARGET_FRAME_MS_ULTRA_CAPABLE;
    }
    applyQualityPresetRuntime('optimal');
}

function advanceQualityBenchmark(frameDt) {
    if (!qualityBenchmarkState) return;

    qualityBenchmarkState.frameCount++;
    if (qualityBenchmarkState.frameCount <= qualityBenchmarkState.warmupFrames) return;

    qualityBenchmarkState.accumulatedDt += frameDt;
    qualityBenchmarkState.sampleCount++;

    if (qualityBenchmarkState.sampleCount < qualityBenchmarkState.sampleFrames) return;

    var avgFrameMs = (qualityBenchmarkState.accumulatedDt / qualityBenchmarkState.sampleCount) * 1000.0;

    if (qualityBenchmarkState.phase === 'optimal') {
        if (avgFrameMs <= qualityBenchmarkState.highProbeGateMs) {
            qualityBenchmarkState.phase = 'high';
            qualityBenchmarkState.optimalAvgMs = avgFrameMs;
            resetQualityBenchmarkCounters(qualityBenchmarkState);
            applyQualityPresetRuntime('high');
            return;
        }

        if (avgFrameMs <= QUALITY_BENCHMARK_TARGET_FRAME_MS) {
            finishQualityBenchmark('optimal', avgFrameMs);
            return;
        }

        qualityBenchmarkState.phase = 'mobile';
        resetQualityBenchmarkCounters(qualityBenchmarkState);
        applyQualityPresetRuntime('mobile');
        return;
    }

    if (qualityBenchmarkState.phase === 'mobile') {
        finishQualityBenchmark('mobile', avgFrameMs);
        return;
    }

    if (qualityBenchmarkState.phase === 'high') {
        if (avgFrameMs <= qualityBenchmarkState.highTargetFrameMs) {
            finishQualityBenchmark('high', avgFrameMs);
        } else {
            finishQualityBenchmark('optimal', qualityBenchmarkState.optimalAvgMs);
        }
    }
}

function isLikelyMobileDevice() {
    var ua = (window.navigator && window.navigator.userAgent) || '';
    var uaMobile = /(Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini|Mobile)/i.test(ua);
    var coarsePointer = !!(window.matchMedia && window.matchMedia('(pointer: coarse)').matches);
    var smallViewport = Math.min(window.innerWidth || 0, window.innerHeight || 0) <= 900;
    return uaMobile || (coarsePointer && smallViewport);
}

function clampResolutionScale(value) {
    return Math.max(0.35, Math.min(2.0, value || 1.0));
}

function halton(index, base) {
    var f = 1.0;
    var r = 0.0;
    var i = index;
    while (i > 0) {
        f /= base;
        r += f * (i % base);
        i = Math.floor(i / base);
    }
    return r;
}

function setupTemporalAA() {
    var ppVertexShader = [
        'varying vec2 vUv;',
        'void main() {',
        '    vUv = uv;',
        '    gl_Position = vec4(position, 1.0);',
        '}'
    ].join('\n');

    var blendFS = [
        'uniform sampler2D tCurrent;',
        'uniform sampler2D tHistory;',
        'uniform float historyWeight;',
        'uniform float historyValid;',
        'uniform float clipBox;',
        'varying vec2 vUv;',
        'void main() {',
        '    vec3 current = texture2D(tCurrent, vUv).rgb;',
        '    vec3 history = texture2D(tHistory, vUv).rgb;',
        '    history = clamp(history, current - vec3(clipBox), current + vec3(clipBox));',
        '    float lumaCurrent = dot(current, vec3(0.299, 0.587, 0.114));',
        '    float lumaHistory = dot(history, vec3(0.299, 0.587, 0.114));',
        '    float reactive = clamp(1.0 - abs(lumaCurrent - lumaHistory) * 5.0, 0.0, 1.0);',
        '    float w = historyWeight * historyValid * reactive;',
        '    gl_FragColor = vec4(mix(current, history, w), 1.0);',
        '}'
    ].join('\n');

    var copyFS = [
        'uniform sampler2D tDiffuse;',
        'varying vec2 vUv;',
        'void main() {',
        '    gl_FragColor = texture2D(tDiffuse, vUv);',
        '}'
    ].join('\n');

    var rtParams = {
        minFilter: THREE.LinearFilter,
        magFilter: THREE.LinearFilter,
        format: THREE.RGBAFormat
    };

    function createTarget(w, h) {
        return new THREE.WebGLRenderTarget(Math.max(1, w), Math.max(1, h), rtParams);
    }

    var ppScene = new THREE.Scene();
    var ppCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
    var ppMesh = new THREE.Mesh(new THREE.PlaneBufferGeometry(2, 2));
    ppScene.add(ppMesh);

    var blendMat = new THREE.ShaderMaterial({
        uniforms: {
            tCurrent: { type: 't', value: null },
            tHistory: { type: 't', value: null },
            historyWeight: { type: 'f', value: 0.0 },
            historyValid: { type: 'f', value: 0.0 },
            clipBox: { type: 'f', value: 0.08 }
        },
        vertexShader: ppVertexShader,
        fragmentShader: blendFS,
        depthWrite: false,
        depthTest: false
    });

    var copyMat = new THREE.ShaderMaterial({
        uniforms: {
            tDiffuse: { type: 't', value: null }
        },
        vertexShader: ppVertexShader,
        fragmentShader: copyFS,
        depthWrite: false,
        depthTest: false
    });

    var pass = {
        ppScene: ppScene,
        ppCamera: ppCamera,
        ppMesh: ppMesh,
        blendMat: blendMat,
        copyMat: copyMat,
        currentRT: createTarget(1, 1),
        historyRT: createTarget(1, 1),
        outputRT: createTarget(1, 1),
        historyValid: false,
        frameIndex: 0,
        jitter: new THREE.Vector2(0, 0),

        reset: function() {
            this.historyValid = false;
            this.frameIndex = 0;
        },

        resize: function(w, h) {
            this.currentRT.dispose();
            this.historyRT.dispose();
            this.outputRT.dispose();
            this.currentRT = createTarget(w, h);
            this.historyRT = createTarget(w, h);
            this.outputRT = createTarget(w, h);
            this.reset();
        },

        nextJitter: function() {
            var idx = (this.frameIndex % 8) + 1;
            this.jitter.set(halton(idx, 2) - 0.5, halton(idx, 3) - 0.5);
            this.frameIndex++;
            return this.jitter;
        },

        render: function(rdr, currentTarget, cameraDelta, taaSettings) {
            var settings = taaSettings || {};
            var baseHistoryWeight = Math.max(0.0, Math.min(0.98,
                settings.history_weight !== undefined ? settings.history_weight : 0.88));
            var baseClip = Math.max(0.01, Math.min(0.5,
                settings.clip_box !== undefined ? settings.clip_box : 0.06));
            var motionRejection = Math.max(0.0, Math.min(20.0,
                settings.motion_rejection !== undefined ? settings.motion_rejection : 8.0));
            var maxCameraDelta = Math.max(0.005, Math.min(0.5,
                settings.max_camera_delta !== undefined ? settings.max_camera_delta : 0.08));
            var motionClipScale = Math.max(0.0, Math.min(2.0,
                settings.motion_clip_scale !== undefined ? settings.motion_clip_scale : 0.6));

            var useHistory = this.historyValid && cameraDelta < maxCameraDelta;
            var motionAttenuation = Math.max(0.0, 1.0 - cameraDelta * motionRejection);
            var historyWeight = useHistory ? baseHistoryWeight * motionAttenuation : 0.0;
            var clip = baseClip + Math.min(cameraDelta * motionClipScale, 0.5);

            this.blendMat.uniforms.tCurrent.value = currentTarget;
            this.blendMat.uniforms.tHistory.value = this.historyRT;
            this.blendMat.uniforms.historyWeight.value = historyWeight;
            this.blendMat.uniforms.historyValid.value = useHistory ? 1.0 : 0.0;
            this.blendMat.uniforms.clipBox.value = clip;
            this.ppMesh.material = this.blendMat;
            rdr.render(this.ppScene, this.ppCamera, this.outputRT, true);

            this.copyMat.uniforms.tDiffuse.value = this.outputRT;
            this.ppMesh.material = this.copyMat;
            rdr.render(this.ppScene, this.ppCamera);

            rdr.render(this.ppScene, this.ppCamera, this.historyRT, true);
            this.historyValid = true;
        }
    };

    return pass;
}

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

// ─── Hover Approach State ─────────────────────────────────────────────────────
// Tracks the hovering animation where a static observer slowly descends toward
// the event horizon under powered flight.  Unlike the freefall dive, the
// observer has ZERO velocity at every radius (they fire thrusters to hover).
// This produces the pure gravitational blueshift described by GR: background
// light from infinity gains energy falling into the potential well.
// At radius r the frequency boost is 1/sqrt(1 - r_s/r), which diverges at
// the horizon — you cannot hover at r = r_s (infinite acceleration needed).
var hoverState = {
    active: false,
    paused: false,
    speed: 0.3,
    currentR: 11.0,
    direction: new THREE.Vector3(1, 0, 0),
    startPosition: new THREE.Vector3(10, 0, 0),
    startVelocity: new THREE.Vector3(0, 1, 0),
    prevMotionState: true,
    prevDistance: 11.0,
    minR: 1.005  // Cannot hover at the horizon (infinite proper acceleration)
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
    isMobileClient = isLikelyMobileDevice();
    if (isMobileClient) {
        document.body.classList.add('mobile-ui');
    }
    var storedQualityPreset = readStoredQualityPreset();
    var initialQualityPreset = storedQualityPreset || 'optimal';
    if (typeof applyQualityPresetValues === 'function') {
        applyQualityPresetValues(shader.parameters, initialQualityPreset);
    } else {
        shader.parameters.quality = initialQualityPreset;
    }

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
        taa_jitter: { type: "v2", value: new THREE.Vector2() },

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

        grav_blueshift_factor: { type: "f", value: 1.0 },

        star_texture: { type: "t", value: textures.stars },
        galaxy_texture: { type: "t", value: textures.galaxy },
        planet_texture: { type: "t", value: textures.moon },
        spectrum_texture: { type: "t", value: textures.spectra }
    };
    shaderUniforms = uniforms;

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

        // Gravitational blueshift factor for background sky:
        // sqrt(1 - r_s/r) = sqrt(1 - 1/r).  Light from infinity has its
        // frequency boosted by 1/grav_blueshift_factor when received by a
        // hovering observer at radius r.  Inside the horizon the concept of
        // a static observer doesn't apply, so we set the factor to 1.0 and
        // let the existing interior_boost shader code handle that regime.
        if (obsR > 1.0) {
            uniforms.grav_blueshift_factor.value =
                Math.sqrt(Math.max(1.0 - 1.0 / obsR, 0.001));
        } else {
            uniforms.grav_blueshift_factor.value = 1.0;
        }

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
        resetTemporalAAHistory();
    };

    scene.updateShader();

    var mesh = new THREE.Mesh( geometry, material );
    scene.add( mesh );

    renderer = new THREE.WebGLRenderer({
        antialias: true,
        powerPreference: 'high-performance'
    });
    renderer.domElement.style.touchAction = 'none';
    container.appendChild( renderer.domElement );

    // ============== BLOOM POST-PROCESSING ==============
    bloomPass = setupBloom();
    // ============== END BLOOM ==============
    taaPass = setupTemporalAA();

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
    cameraControls.zoomCallback = function(dollyDeltaY) {
        // Pinch out (positive delta) zooms in; pinch in zooms out.
        var zoomFactor = dollyDeltaY > 0 ? 0.93 : 1.08;
        var newDist = shader.parameters.observer.distance * zoomFactor;
        newDist = Math.max(1.5, Math.min(30, newDist));
        shader.parameters.observer.distance = newDist;
        updateCamera();
        shader.needsUpdate = true;
        if (distanceController) distanceController.updateDisplay();
        return true;
    };
    cameraControls.addEventListener( 'change', updateCamera );
    updateCamera();

    applyRenderScaleFromSettings();

    window.addEventListener( 'resize', onWindowResize, false );

    setupGUI();
    beginQualityBenchmarkIfNeeded();
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

    // Abort any active hover first — restores observer to pre-hover state so
    // diveState saves the correct original position/velocity below.
    if (hoverState.active) {
        resetHover();
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

// ─── Hover Approach Animation Functions ─────────────────────────────────────
// Physics: Powered hovering at each radius in Schwarzschild geometry.
// The observer is STATIONARY (v = 0) at every point — thrusters fire to
// counteract gravity.  Required proper acceleration:
//   a = M / (r² √(1 - r_s/r)) = 0.5 / (r² √(1 - 1/r))
// which diverges at r → r_s.  The pure gravitational blueshift of
// background light is f_obs/f_emit = 1/√(1 - r_s/r).
//
// Unlike the freefall dive, the observer has zero kinematic Doppler, so the
// full gravitational blueshift is visible: the sky progressively shifts toward
// blue/UV with an intensity boost of D³ (Liouville invariant) as the observer
// descends.  The approaching motion itself is an instantaneous quasi-static
// sequence of hovering positions (not a free-fall trajectory).
// ─────────────────────────────────────────────────────────────────────────────

function startHover() {
    // Toggle pause if already active
    if (hoverState.active && !hoverState.paused) {
        hoverState.paused = true;
        updateHoverUI();
        return;
    }
    if (hoverState.paused) {
        hoverState.paused = false;
        updateHoverUI();
        return;
    }

    // Abort any active dive first
    if (diveState.active || diveState.reachedSingularity) {
        resetDive();
    }

    // Save current observer state for reset
    hoverState.prevMotionState = shader.parameters.observer.motion;
    hoverState.prevDistance = shader.parameters.observer.distance;
    hoverState.startPosition = observer.position.clone();
    hoverState.startVelocity = observer.velocity.clone();

    // Disable orbital motion — hover controls the observer now
    shader.parameters.observer.motion = false;

    // Hover direction = radially inward from current position
    hoverState.direction = observer.position.clone().normalize();
    hoverState.currentR = observer.position.length();
    hoverState.active = true;
    hoverState.paused = false;

    // Set observer as stationary (hovering)
    observer.velocity.set(0, 0, 0);

    scene.updateShader();
    updateCamera();
    updateHoverUI();
    if (refreshAllControllersGlobal) refreshAllControllersGlobal();
}

function resetHover() {
    hoverState.active = false;
    hoverState.paused = false;

    // Restore pre-hover observer state
    shader.parameters.observer.motion = hoverState.prevMotionState;
    shader.parameters.observer.distance = hoverState.prevDistance;
    hoverState.currentR = hoverState.prevDistance;

    observer.position.copy(hoverState.startPosition);
    observer.velocity.copy(hoverState.startVelocity);

    scene.updateShader();
    updateCamera();
    shader.needsUpdate = true;
    updateHoverUI();
    if (refreshAllControllersGlobal) refreshAllControllersGlobal();
}

function seekHover(targetR) {
    if (!hoverState.active) return;
    targetR = Math.max(hoverState.minR, Math.min(hoverState.prevDistance, targetR));

    hoverState.currentR = targetR;
    hoverState.paused = true;

    // Update observer position (stationary, zero velocity)
    observer.position.copy(hoverState.direction.clone().multiplyScalar(targetR));
    observer.velocity.set(0, 0, 0);  // Hovering = stationary
    shader.parameters.observer.distance = targetR;

    shader.needsUpdate = true;
    updateCamera();
    updateHoverUI();
}

function updateHover(dt) {
    if (!hoverState.active || hoverState.paused) return;

    var r = hoverState.currentR;
    if (r <= hoverState.minR) {
        hoverState.paused = true;
        updateHoverUI();
        return;
    }

    // Controlled quasi-static descent: the approach rate scales as
    // (r - minR) so the observer naturally decelerates as they
    // approach the minimum hoverable radius.
    var approachRate = hoverState.speed *
        Math.max(r - hoverState.minR, 0.001) *
        shader.parameters.time_scale;
    var newR = Math.max(r - approachRate * dt, hoverState.minR);

    hoverState.currentR = newR;

    // Update observer position (stationary, zero velocity)
    observer.position.copy(hoverState.direction.clone().multiplyScalar(newR));
    observer.velocity.set(0, 0, 0);

    shader.parameters.observer.distance = newR;

    // Advance observer time with gravitational time dilation.
    // For a hovering observer: dτ/dt = √(1 - r_s/r) = √(1 - 1/r)
    var timeDilation = Math.sqrt(Math.max(1.0 - 1.0 / newR, 0.001));
    observer.time += dt * shader.parameters.time_scale / timeDilation;

    shader.needsUpdate = true;
    updateHoverUI();
}

function updateHoverUI() {
    var radiusEl = document.getElementById('hover-radius');
    var blueshiftEl = document.getElementById('hover-blueshift');
    var accelEl = document.getElementById('hover-accel');
    var statusEl = document.getElementById('hover-status');
    var btnEl = document.getElementById('hover-start-btn');
    var resetBtn = document.getElementById('hover-reset-btn');
    var horizonBar = document.getElementById('hover-horizon-bar');

    if (!radiusEl) return;

    var r = hoverState.currentR;

    // Gravitational blueshift factor: f_obs/f_emit = 1/√(1 - 1/r)
    var gravFactor = Math.sqrt(Math.max(1.0 - 1.0 / r, 0.001));
    var blueshift = 1.0 / gravFactor;

    // Required proper acceleration to hover: a = M/(r²√(1 - r_s/r))
    // In units where M = 0.5, r_s = 1:
    var properAccel = 0.5 / (r * r * gravFactor);

    radiusEl.innerHTML = 'r = ' + r.toFixed(3) + ' r<sub>s</sub>';
    blueshiftEl.innerHTML = 'z<sub>grav</sub> = ' + blueshift.toFixed(2) + '\u00d7';
    accelEl.innerHTML = 'a = ' + (properAccel < 100 ? properAccel.toFixed(2) : properAccel.toFixed(0)) +
        ' c\u00b2/r<sub>s</sub>';

    // Update horizon proximity bar
    if (horizonBar) {
        var progress = Math.max(0, Math.min(100,
            (1.0 - r / Math.max(hoverState.prevDistance, 1)) * 100));
        horizonBar.style.width = progress + '%';
        if (r < 1.5) {
            horizonBar.className = 'hover-horizon-fill near';
        } else {
            horizonBar.className = 'hover-horizon-fill normal';
        }
    }

    // Show effective speed
    var speedEl = document.getElementById('hover-speed-val');
    if (speedEl && hoverState.active) {
        var effSpd = hoverState.speed;
        speedEl.textContent = (effSpd < 0.1 ? effSpd.toFixed(2) : effSpd.toFixed(1)) + '\u00d7';
    }

    if (!hoverState.active) {
        statusEl.textContent = 'Ready';
        statusEl.className = 'hover-status ready';
        btnEl.textContent = '\u25b6 START HOVER';
        btnEl.disabled = false;
    } else if (hoverState.paused && r <= hoverState.minR + 0.01) {
        statusEl.textContent = '\u26a0 Minimum hover radius';
        statusEl.className = 'hover-status min-radius';
        btnEl.textContent = '\u25b6 START HOVER';
        btnEl.disabled = true;
    } else if (hoverState.paused) {
        statusEl.textContent = '\u23f8 Hovering at r = ' + r.toFixed(2);
        statusEl.className = 'hover-status paused';
        btnEl.textContent = '\u25b6 RESUME';
    } else if (r > 3.0) {
        statusEl.textContent = '\u2193 Descending — mild blueshift';
        statusEl.className = 'hover-status descending';
        btnEl.textContent = '\u23f8 PAUSE';
    } else if (r > 1.5) {
        statusEl.textContent = '\u26a1 Strong blueshift zone';
        statusEl.className = 'hover-status strong';
        btnEl.textContent = '\u23f8 PAUSE';
    } else {
        statusEl.textContent = '\ud83d\udca0 Extreme blueshift!';
        statusEl.className = 'hover-status extreme';
        btnEl.textContent = '\u23f8 PAUSE';
    }

    if (resetBtn) {
        resetBtn.disabled = !hoverState.active;
    }
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

function resizeRendererAndPasses() {
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

applyRenderScaleFromSettings = function() {
    resizeRendererAndPasses();
    resetTemporalAAHistory();
    if (updateUniforms) updateUniforms();
    if (shader) shader.needsUpdate = true;
};

resetTemporalAAHistory = function() {
    if (taaPass) taaPass.reset();
    if (shaderUniforms && shaderUniforms.taa_jitter) {
        shaderUniforms.taa_jitter.value.set(0, 0);
    }
    lastTaaCameraMat.identity();
};

function onWindowResize( event ) {
    resizeRendererAndPasses();
    resetTemporalAAHistory();
    updateUniforms();
}

var lastCameraMat = new THREE.Matrix4().identity();

// ─── Frame timing ─────────────────────────────────────────────────────────────
// Always called once per RAF tick (inside animate()), never inside render().
// Capping at MAX_FRAME_DT means a tab-switch or slow frame can never hand a
// multi-second delta to the physics integration.
var getFrameDuration = (function() {
    var MAX_FRAME_DT = 0.1; // seconds — max delta per frame
    var _now = (typeof performance !== 'undefined' && performance.now)
        ? function() { return performance.now(); }
        : function() { return new Date().getTime(); };
    var lastTimestamp = _now();

    function resetClock() { lastTimestamp = _now(); }

    // Reset on tab-show so the first resumed frame gets ~0 dt, not accumulated time.
    document.addEventListener('visibilitychange', function() {
        if (!document.hidden) resetClock();
    });
    // Also handle bfcache page restore (navigating back/forward).
    window.addEventListener('pageshow', function(e) {
        if (e.persisted) resetClock();
    });

    return function() {
        var now = _now();
        var diff = (now - lastTimestamp) / 1000.0;
        lastTimestamp = now;
        return Math.min(diff, MAX_FRAME_DT);
    };
})();
// ─────────────────────────────────────────────────────────────────────────────

function animate() {
    requestAnimationFrame( animate );

    // ── Advance simulation time unconditionally every RAF frame ───────────────
    // This MUST happen outside the lazy-render gate so observer.time (and the
    // shader's time uniform) always tracks real-world elapsed time, even when
    // nothing in the scene changes visually (static camera, no orbital motion).
    // Before this was inside render(), so hidden-tab pauses or post-dive states
    // with hasMovingParts()==false silently froze the disk animation.
    var dt = getFrameDuration();
    if (diveState.active && !diveState.paused && !diveState.reachedSingularity) {
        updateDive(dt);
        updateCamera();
    } else if (hoverState.active && !hoverState.paused) {
        updateHover(dt);
        updateCamera();
    } else {
        observer.move(dt);
        if (shader.parameters.observer.motion) updateCamera();
    }
    advanceQualityBenchmark(dt);
    // ─────────────────────────────────────────────────────────────────────────

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

function renderSceneToTarget(target) {
    if (shader.parameters.bloom.enabled && bloomPass) {
        bloomPass.render(renderer, scene, camera, shader.parameters.bloom, target);
    } else if (target) {
        renderer.render(scene, camera, target, true);
    } else {
        renderer.render(scene, camera);
    }
}

function render() {
    var taaEnabled = !!shader.parameters.taa_enabled && !!taaPass;
    if (shaderUniforms && shaderUniforms.taa_jitter) {
        if (taaEnabled) {
            var jitter = taaPass.nextJitter();
            shaderUniforms.taa_jitter.value.set(jitter.x, jitter.y);
        } else {
            shaderUniforms.taa_jitter.value.set(0, 0);
        }
    }

    // Time advancement has already been done in animate(); render() only draws.
    updateUniforms();

    if (taaEnabled) {
        renderSceneToTarget(taaPass.currentRT);
        taaPass.render(
            renderer,
            taaPass.currentRT,
            frobeniusDistance(camera.matrixWorldInverse, lastTaaCameraMat),
            shader.parameters.taa
        );
        lastTaaCameraMat.copy(camera.matrixWorldInverse);
    } else {
        renderSceneToTarget(null);
        if (taaPass && taaPass.historyValid) taaPass.reset();
        lastTaaCameraMat.copy(camera.matrixWorldInverse);
    }
    updateAxesGizmo();
}
