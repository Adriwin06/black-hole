"use strict";

import { observer, shader, scene } from '../runtime/runtime-state.js';
import { diveState, hoverState, animationTimelineCaptureState } from './scenario-state.js';
import { updateCamera } from '../../scene/camera.js';
import { refreshRendererUiBindings } from '../runtime/ui-bindings.js';
import {
    readAnimationCaptureAnchorPosition,
    readAnimationCaptureAnchorVelocity,
    finalizeAnimationTimelineCapture
} from './animation-capture.js';
import { resetHover } from './hover.js';

export function advanceTimelineDrivenDiveObserverTime(dt) {
    if (!isFinite(dt) || dt <= 0 || !diveState.active || diveState.reachedSingularity) {
        return;
    }
    var r = Math.max(diveState.currentR, 0.08);
    var effectiveSpeed = diveState.cinematic
        ? diveState.speed * cinematicFactor(r)
        : diveState.speed;
    observer.time += dt * effectiveSpeed * shader.parameters.time_scale;
}

export function startDive(options) {
    options = options || {};
    var restartRequested = !!options.restart;

    if (!restartRequested && diveState.active && !diveState.paused) {
        diveState.paused = true;
        diveState.timelineDriven = false;
        updateDiveUI();
        return;
    }
    if (!restartRequested && diveState.paused) {
        diveState.paused = false;
        diveState.timelineDriven = false;
        updateDiveUI();
        return;
    }
    if (restartRequested && (diveState.active || diveState.reachedSingularity)) {
        resetDive();
    }

    if (hoverState.active) {
        resetHover();
    }

    var anchorPosition = readAnimationCaptureAnchorPosition(options.anchorPosition) ||
        observer.position.clone();
    if (anchorPosition.lengthSq() < 1e-10) {
        anchorPosition.set(Math.max(shader.parameters.observer.distance, 1.0), 0, 0);
    }
    var anchorVelocity = readAnimationCaptureAnchorVelocity(options.anchorVelocity) ||
        observer.velocity.clone();
    var anchorRadius = anchorPosition.length();

    observer.position.copy(anchorPosition);
    observer.velocity.copy(anchorVelocity);
    shader.parameters.observer.distance = anchorRadius;
    if (typeof options.observerTime === 'number' && isFinite(options.observerTime)) {
        observer.time = options.observerTime;
    }

    diveState.prevMotionState = (options.prevMotionState !== undefined)
        ? !!options.prevMotionState
        : shader.parameters.observer.motion;
    diveState.prevDistance = (typeof options.prevDistance === 'number' &&
        isFinite(options.prevDistance))
        ? options.prevDistance
        : anchorRadius;
    diveState.startPosition = anchorPosition.clone();
    diveState.startVelocity = anchorVelocity.clone();
    diveState.startRenderSettings = {
        n_steps: shader.parameters.n_steps,
        max_revolutions: shader.parameters.max_revolutions,
        rk4_integration: shader.parameters.rk4_integration
    };

    shader.parameters.observer.motion = false;

    diveState.direction = observer.position.clone().normalize();
    diveState.currentR = observer.position.length();
    diveState.active = true;
    diveState.paused = false;
    diveState.timelineDriven = false;
    diveState.reachedSingularity = false;

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
    refreshRendererUiBindings();
}

export function resetDive() {
    if (animationTimelineCaptureState.active &&
        animationTimelineCaptureState.mode === 'dive') {
        finalizeAnimationTimelineCapture('dive');
    }
    diveState.active = false;
    diveState.paused = false;
    diveState.timelineDriven = false;
    diveState.reachedSingularity = false;

    shader.parameters.observer.motion = diveState.prevMotionState;
    shader.parameters.observer.distance = diveState.prevDistance;
    diveState.currentR = diveState.prevDistance;

    observer.position.copy(diveState.startPosition);
    observer.velocity.copy(diveState.startVelocity);

    if (diveState.startRenderSettings) {
        shader.parameters.n_steps = diveState.startRenderSettings.n_steps;
        shader.parameters.max_revolutions = diveState.startRenderSettings.max_revolutions;
        shader.parameters.rk4_integration = diveState.startRenderSettings.rk4_integration;
        diveState.startRenderSettings = null;
    }

    scene.updateShader();
    updateCamera();
    shader.needsUpdate = true;
    updateDiveUI();
    updateDiveFade();
    refreshRendererUiBindings();
}

function cinematicFactor(r) {
    var farBoost = 2.0 * Math.max(r - 3.0, 0.0) / 7.0;
    var photonSlow = 3.0 * Math.exp(-Math.pow((r - 1.5) / 0.30, 2));
    var horizonSlow = 5.0 * Math.exp(-Math.pow((r - 1.0) / 0.22, 2));
    return (1.0 + farBoost) / (1.0 + photonSlow + horizonSlow);
}

export function seekDive(targetR, options) {
    options = options || {};
    if (!diveState.active && !diveState.reachedSingularity) return;
    targetR = Math.max(0.08, Math.min(diveState.prevDistance, targetR));

    if (diveState.reachedSingularity && targetR > 0.12) {
        diveState.reachedSingularity = false;
        diveState.active = true;
    }

    diveState.currentR = targetR;
    if (targetR <= 0.09) {
        diveState.reachedSingularity = true;
    }
    diveState.paused = true;
    diveState.timelineDriven = !!options.timelineDriven;

    observer.position.copy(diveState.direction.clone().multiplyScalar(targetR));
    var velocity = Math.min(Math.sqrt(1.0 / targetR), 0.998);
    observer.velocity.copy(diveState.direction.clone().multiplyScalar(-velocity));
    shader.parameters.observer.distance = targetR;

    shader.needsUpdate = true;
    updateCamera();
    updateDiveUI();
}

export function updateDive(dt) {
    if (!diveState.active || diveState.paused || diveState.reachedSingularity) return;

    var r = diveState.currentR;
    if (r < 0.08) {
        diveState.reachedSingularity = true;
        updateDiveUI();
        return;
    }

    var effectiveSpeed = diveState.cinematic
        ? diveState.speed * cinematicFactor(r)
        : diveState.speed;
    var fallDt = dt * effectiveSpeed * shader.parameters.time_scale;
    var k1 = -Math.sqrt(1.0 / r) * fallDt;
    var rMid = Math.max(r + k1 * 0.5, 0.01);
    var k2 = -Math.sqrt(1.0 / rMid) * fallDt;
    var newR = Math.max(r + k2, 0.08);

    diveState.currentR = newR;
    observer.position.copy(diveState.direction.clone().multiplyScalar(newR));

    var velocity = Math.min(Math.sqrt(1.0 / newR), 0.998);
    observer.velocity.copy(diveState.direction.clone().multiplyScalar(-velocity));

    shader.parameters.observer.distance = newR;
    observer.time += dt * effectiveSpeed * shader.parameters.time_scale;

    if (newR < 1.0 && r >= 1.0) {
        scene.updateShader();
    }

    shader.needsUpdate = true;
    updateDiveUI();
}

export function updateDiveUI() {
    var radiusEl = document.getElementById('dive-radius');
    var velocityEl = document.getElementById('dive-velocity');
    var statusEl = document.getElementById('dive-status');
    var btnEl = document.getElementById('dive-start-btn');
    var resetBtn = document.getElementById('dive-reset-btn');
    var horizonBar = document.getElementById('dive-horizon-bar');

    if (!radiusEl) return;

    var r = diveState.currentR;
    var velocity = r > 0.01 ? Math.min(Math.sqrt(1.0 / r), 0.999) : 0.999;

    radiusEl.innerHTML = 'r = ' + r.toFixed(3) + ' r<sub>s</sub>';
    velocityEl.textContent = 'v = ' + velocity.toFixed(3) + ' c';

    var speedEl = document.getElementById('dive-speed-val');
    if (speedEl && diveState.active) {
        var effectiveSpeed = diveState.cinematic
            ? diveState.speed * cinematicFactor(r) : diveState.speed;
        speedEl.textContent = (effectiveSpeed < 0.1
            ? effectiveSpeed.toFixed(2) : effectiveSpeed.toFixed(1)) + 'x';
    }

    if (horizonBar) {
        var progress = Math.max(0, Math.min(100,
            (1.0 - r / Math.max(diveState.prevDistance, 1)) * 100));
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

export function updateDiveFade() {
    // The interior darkening is driven by the actual ray tracing result.
}


