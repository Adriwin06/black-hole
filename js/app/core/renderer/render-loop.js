"use strict";

import { THREE } from '../../vendor.js';
import {
    renderer,
    camera,
    scene,
    stats,
    shader,
    observer,
    bloomPass,
    taaPass,
    shaderUniforms,
    lastTaaCameraMat,
    rendererContextLost,
    updateUniforms,
    resetTemporalAAHistory
} from '../runtime/runtime-state.js';
import { diveState, hoverState } from '../scenarios/scenario-state.js';
import {
    updatePresentation,
    getPresentationState,
    updatePresentationOverlay
} from '../../presentation/runtime/presentation-controller.js';
import {
    updateDive,
    advanceTimelineDrivenDiveObserverTime
} from '../scenarios/dive.js';
import {
    updateHover,
    advanceTimelineDrivenHoverObserverTime
} from '../scenarios/hover.js';
import { updateCamera, frobeniusDistance } from '../../scene/camera.js';
import { advanceQualityBenchmark } from './renderer-quality.js';
import { updateAnimationTimelineCaptureFrame } from '../scenarios/animation-capture.js';
import { resizeRendererAndPasses, updateAxesGizmo } from './renderer.js';
import { registerBlackHoleRuntimeApi } from '../runtime/runtime-registry.js';

var lastCameraMat = new THREE.Matrix4().identity();
var resetRendererFrameClock = function() {};
var rendererOfflineSteppingActive = false;

var getFrameDuration = (function() {
    var MAX_FRAME_DT = 0.1;
    var nowFn = (typeof performance !== 'undefined' && performance.now)
        ? function() { return performance.now(); }
        : function() { return new Date().getTime(); };
    var lastTimestamp = nowFn();

    function resetClock() {
        lastTimestamp = nowFn();
    }
    resetRendererFrameClock = resetClock;

    document.addEventListener('visibilitychange', function() {
        if (!document.hidden) resetClock();
    });
    window.addEventListener('pageshow', function(e) {
        if (e.persisted) resetClock();
    });

    return function() {
        var now = nowFn();
        var diff = (now - lastTimestamp) / 1000.0;
        lastTimestamp = now;
        return Math.min(diff, MAX_FRAME_DT);
    };
})();

function stepRendererSimulation(dt, skipBenchmark) {
    var presentationRuntimeState = getPresentationState();
    if (presentationRuntimeState &&
        presentationRuntimeState.playing) {
        updatePresentation(dt);
    }
    var presentationDrivesObserverTime = !!(
        presentationRuntimeState &&
        presentationRuntimeState.playing &&
        presentationRuntimeState.drives_observer_time
    );
    if (diveState.active && !diveState.reachedSingularity) {
        if (!diveState.paused && !diveState.timelineDriven) {
            updateDive(dt);
        } else if (diveState.timelineDriven &&
            presentationRuntimeState &&
            presentationRuntimeState.playing &&
            !presentationDrivesObserverTime) {
            advanceTimelineDrivenDiveObserverTime(dt);
        }
        if (!diveState.paused || diveState.timelineDriven) {
            updateCamera();
        }
    } else if (hoverState.active) {
        if (!hoverState.paused && !hoverState.timelineDriven) {
            updateHover(dt);
        } else if (hoverState.timelineDriven &&
            presentationRuntimeState &&
            presentationRuntimeState.playing &&
            !presentationDrivesObserverTime) {
            advanceTimelineDrivenHoverObserverTime(dt);
        }
        if (!hoverState.paused || hoverState.timelineDriven) {
            updateCamera();
        }
    } else {
        observer.move(dt);
        if (shader.parameters.observer.motion) updateCamera();
    }
    if (!skipBenchmark) {
        advanceQualityBenchmark(dt);
    }
}

function drawRendererFrame(forceRender) {
    camera.updateMatrixWorld();
    camera.matrixWorldInverse.getInverse(camera.matrixWorld);
    updatePresentationOverlay();

    if (forceRender || shader.needsUpdate || shader.hasMovingParts() ||
        frobeniusDistance(camera.matrixWorldInverse, lastCameraMat) > 1e-10) {
        shader.needsUpdate = false;
        render();
        lastCameraMat = camera.matrixWorldInverse.clone();
    }
    stats.update();
}

export function setRendererOfflineSteppingActive(active) {
    rendererOfflineSteppingActive = !!active;
    if (!rendererOfflineSteppingActive) {
        resetRendererFrameClock();
    }
}

export function stepRendererForOfflineRecording(dt) {
    var frameDt = parseFloat(dt);
    if (!isFinite(frameDt) || frameDt <= 0) frameDt = 1.0 / 60.0;
    stepRendererSimulation(frameDt, true);
    drawRendererFrame(true);
}

export function animate() {
    requestAnimationFrame(animate);

    if (rendererOfflineSteppingActive) {
        stats.update();
        return;
    }

    var dt = getFrameDuration();
    stepRendererSimulation(dt, false);
    updateAnimationTimelineCaptureFrame();
    drawRendererFrame(false);
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

export var blackHoleRendererRuntime = {
    setOfflineSteppingActive: setRendererOfflineSteppingActive,
    stepOfflineFrame: stepRendererForOfflineRecording,
    resetFrameClock: function() {
        resetRendererFrameClock();
    },
    isContextLost: function() {
        return rendererContextLost;
    },
    resizeForOfflineRecording: function(w, h) {
        if (!renderer) return false;
        renderer.setPixelRatio(1);
        renderer.setSize(w, h, false);
        renderer.domElement.style.objectFit = 'contain';
        var rw = renderer.domElement.width;
        var rh = renderer.domElement.height;
        if (bloomPass) bloomPass.resize(rw, rh);
        if (taaPass) taaPass.resize(rw, rh);
        if (camera) {
            camera.aspect = w / h;
            camera.updateProjectionMatrix();
        }
        var taaWarmupFrames = shader.parameters.taa_enabled ? 64 : 0;
        for (var i = 0; i < taaWarmupFrames; i++) {
            drawRendererFrame(true);
        }
        return true;
    },
    restoreWindowSizeAfterRecording: function() {
        renderer.domElement.style.objectFit = '';
        resizeRendererAndPasses();
        resetTemporalAAHistory();
        if (camera) {
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
        }
    }
};

window.blackHoleRendererRuntime = blackHoleRendererRuntime;
registerBlackHoleRuntimeApi('renderer', blackHoleRendererRuntime);


