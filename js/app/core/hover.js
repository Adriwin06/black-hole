"use strict";

import { observer, shader, scene } from './runtime-state.js';
import { diveState, hoverState, animationTimelineCaptureState } from './scenario-state.js';
import { updateCamera } from '../scene/camera.js';
import { refreshRendererUiBindings } from './ui-bindings.js';
import {
    readAnimationCaptureAnchorPosition,
    readAnimationCaptureAnchorVelocity,
    finalizeAnimationTimelineCapture
} from './animation-capture.js';
import { resetDive } from './dive.js';

export function advanceTimelineDrivenHoverObserverTime(dt) {
    if (!isFinite(dt) || dt <= 0 || !hoverState.active) return;
    var r = Math.max(hoverState.currentR, hoverState.minR, 1.0001);
    var timeDilation = Math.sqrt(Math.max(1.0 - 1.0 / r, 0.001));
    observer.time += dt * shader.parameters.time_scale / timeDilation;
}

export function startHover(options) {
    options = options || {};
    var restartRequested = !!options.restart;

    if (!restartRequested && hoverState.active && !hoverState.paused) {
        hoverState.paused = true;
        hoverState.timelineDriven = false;
        updateHoverUI();
        return;
    }
    if (!restartRequested && hoverState.paused) {
        hoverState.paused = false;
        hoverState.timelineDriven = false;
        updateHoverUI();
        return;
    }
    if (restartRequested && hoverState.active) {
        resetHover();
    }

    if (diveState.active || diveState.reachedSingularity) {
        resetDive();
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

    hoverState.prevMotionState = (options.prevMotionState !== undefined)
        ? !!options.prevMotionState
        : shader.parameters.observer.motion;
    hoverState.prevDistance = (typeof options.prevDistance === 'number' &&
        isFinite(options.prevDistance))
        ? options.prevDistance
        : anchorRadius;
    hoverState.startPosition = anchorPosition.clone();
    hoverState.startVelocity = anchorVelocity.clone();

    shader.parameters.observer.motion = false;

    hoverState.direction = observer.position.clone().normalize();
    hoverState.currentR = observer.position.length();
    hoverState.active = true;
    hoverState.paused = false;
    hoverState.timelineDriven = false;

    observer.velocity.set(0, 0, 0);

    scene.updateShader();
    updateCamera();
    updateHoverUI();
    refreshRendererUiBindings();
}

export function resetHover() {
    if (animationTimelineCaptureState.active &&
        animationTimelineCaptureState.mode === 'hover') {
        finalizeAnimationTimelineCapture('hover');
    }
    hoverState.active = false;
    hoverState.paused = false;
    hoverState.timelineDriven = false;

    shader.parameters.observer.motion = hoverState.prevMotionState;
    shader.parameters.observer.distance = hoverState.prevDistance;
    hoverState.currentR = hoverState.prevDistance;

    observer.position.copy(hoverState.startPosition);
    observer.velocity.copy(hoverState.startVelocity);

    scene.updateShader();
    updateCamera();
    shader.needsUpdate = true;
    updateHoverUI();
    refreshRendererUiBindings();
}

export function seekHover(targetR, options) {
    options = options || {};
    if (!hoverState.active) return;
    targetR = Math.max(hoverState.minR, Math.min(hoverState.prevDistance, targetR));

    hoverState.currentR = targetR;
    hoverState.paused = true;
    hoverState.timelineDriven = !!options.timelineDriven;

    observer.position.copy(hoverState.direction.clone().multiplyScalar(targetR));
    observer.velocity.set(0, 0, 0);
    shader.parameters.observer.distance = targetR;

    shader.needsUpdate = true;
    updateCamera();
    updateHoverUI();
}

export function updateHover(dt) {
    if (!hoverState.active || hoverState.paused) return;

    var r = hoverState.currentR;
    if (r <= hoverState.minR) {
        hoverState.paused = true;
        updateHoverUI();
        return;
    }

    var approachRate = hoverState.speed *
        Math.max(r - hoverState.minR, 0.001) *
        shader.parameters.time_scale;
    var newR = Math.max(r - approachRate * dt, hoverState.minR);

    hoverState.currentR = newR;

    observer.position.copy(hoverState.direction.clone().multiplyScalar(newR));
    observer.velocity.set(0, 0, 0);

    shader.parameters.observer.distance = newR;

    var timeDilation = Math.sqrt(Math.max(1.0 - 1.0 / newR, 0.001));
    observer.time += dt * shader.parameters.time_scale / timeDilation;

    shader.needsUpdate = true;
    updateHoverUI();
}

export function updateHoverUI() {
    var radiusEl = document.getElementById('hover-radius');
    var blueshiftEl = document.getElementById('hover-blueshift');
    var accelEl = document.getElementById('hover-accel');
    var statusEl = document.getElementById('hover-status');
    var btnEl = document.getElementById('hover-start-btn');
    var resetBtn = document.getElementById('hover-reset-btn');
    var horizonBar = document.getElementById('hover-horizon-bar');

    if (!radiusEl) return;

    var r = hoverState.currentR;
    var gravFactor = Math.sqrt(Math.max(1.0 - 1.0 / r, 0.001));
    var blueshift = 1.0 / gravFactor;
    var properAccel = 0.5 / (r * r * gravFactor);

    radiusEl.innerHTML = 'r = ' + r.toFixed(3) + ' r<sub>s</sub>';
    blueshiftEl.innerHTML = 'D<sub>grav</sub> = ' + blueshift.toFixed(2) + '\u00d7';
    accelEl.innerHTML = 'a = ' +
        (properAccel < 100 ? properAccel.toFixed(2) : properAccel.toFixed(0)) +
        ' c\u00b2/r<sub>s</sub>';

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

    var speedEl = document.getElementById('hover-speed-val');
    if (speedEl && hoverState.active) {
        var effectiveSpeed = hoverState.speed;
        speedEl.textContent = (effectiveSpeed < 0.1
            ? effectiveSpeed.toFixed(2) : effectiveSpeed.toFixed(1)) + '\u00d7';
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
        statusEl.textContent = '\u2193 Descending - mild blueshift';
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
