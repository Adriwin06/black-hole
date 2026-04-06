"use strict";

import { THREE } from '../vendor.js';

export function halton(index, base) {
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

export function setupTemporalAA() {
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

            var tmp = this.historyRT;
            this.historyRT = this.outputRT;
            this.outputRT = tmp;
            this.historyValid = true;
        }
    };

    return pass;
}
