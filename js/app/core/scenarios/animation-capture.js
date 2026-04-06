"use strict";

import { THREE } from '../../vendor.js';
import { camera, cameraPan, observer, shader } from '../runtime/runtime-state.js';
import { diveState, hoverState, animationTimelineCaptureState } from './scenario-state.js';
import { getTimelinePanelBinding } from '../runtime/ui-bindings.js';
import { startDive, updateDiveUI } from './dive.js';
import { startHover, updateHoverUI } from './hover.js';
import {
    getPresentationState,
    pausePresentation
} from '../../presentation/runtime/presentation-controller.js';

function cloneVector3Plain(vec) {
    return {
        x: vec ? vec.x : 0,
        y: vec ? vec.y : 0,
        z: vec ? vec.z : 0
    };
}

function cloneQuaternionPlain(quat) {
    return {
        x: quat ? quat.x : 0,
        y: quat ? quat.y : 0,
        z: quat ? quat.z : 0,
        w: quat ? quat.w : 1
    };
}

function cloneAnimationTimelineCaptureSample(sample) {
    return {
        t: sample.t,
        radius: sample.radius,
        observerTime: sample.observerTime,
        cameraPanX: sample.cameraPanX,
        cameraPanY: sample.cameraPanY,
        cameraPosition: cloneVector3Plain(sample.cameraPosition),
        cameraQuaternion: cloneQuaternionPlain(sample.cameraQuaternion)
    };
}

function animationCaptureQuaternionDot(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

function normalizeAnimationCaptureQuaternion(quat) {
    var len = Math.sqrt(
        quat.x * quat.x +
        quat.y * quat.y +
        quat.z * quat.z +
        quat.w * quat.w
    );
    if (len < 1e-8) {
        return { x: 0, y: 0, z: 0, w: 1 };
    }
    return {
        x: quat.x / len,
        y: quat.y / len,
        z: quat.z / len,
        w: quat.w / len
    };
}

function alignAnimationCaptureQuaternion(refQuat, sampleQuat) {
    var q = cloneQuaternionPlain(sampleQuat);
    if (animationCaptureQuaternionDot(refQuat, q) < 0) {
        q.x = -q.x;
        q.y = -q.y;
        q.z = -q.z;
        q.w = -q.w;
    }
    return q;
}

function smoothAnimationTimelineCaptureSamples(samples, radius) {
    var src = Array.isArray(samples) ? samples : [];
    radius = Math.max(0, Math.floor(radius || 0));
    if (!src.length || radius <= 0) {
        return src.slice();
    }

    var out = [];
    for (var i = 0; i < src.length; i++) {
        var base = cloneAnimationTimelineCaptureSample(src[i]);
        if (i === 0 || i === src.length - 1) {
            out.push(base);
            continue;
        }

        var panX = 0, panY = 0;
        var posX = 0, posY = 0, posZ = 0;
        var quatX = 0, quatY = 0, quatZ = 0, quatW = 0;
        var totalWeight = 0;
        var refQuat = normalizeAnimationCaptureQuaternion(base.cameraQuaternion);

        for (var j = Math.max(0, i - radius); j <= Math.min(src.length - 1, i + radius); j++) {
            var neighbor = src[j];
            if (!neighbor) continue;
            var weight = radius + 1 - Math.abs(j - i);
            var neighborQuat = alignAnimationCaptureQuaternion(refQuat, neighbor.cameraQuaternion);

            panX += (isFinite(neighbor.cameraPanX) ? neighbor.cameraPanX : 0) * weight;
            panY += (isFinite(neighbor.cameraPanY) ? neighbor.cameraPanY : 0) * weight;
            posX += (neighbor.cameraPosition && isFinite(neighbor.cameraPosition.x)
                ? neighbor.cameraPosition.x : 0) * weight;
            posY += (neighbor.cameraPosition && isFinite(neighbor.cameraPosition.y)
                ? neighbor.cameraPosition.y : 0) * weight;
            posZ += (neighbor.cameraPosition && isFinite(neighbor.cameraPosition.z)
                ? neighbor.cameraPosition.z : 0) * weight;
            quatX += neighborQuat.x * weight;
            quatY += neighborQuat.y * weight;
            quatZ += neighborQuat.z * weight;
            quatW += neighborQuat.w * weight;
            totalWeight += weight;
        }

        if (totalWeight > 0) {
            base.cameraPanX = panX / totalWeight;
            base.cameraPanY = panY / totalWeight;
            base.cameraPosition.x = posX / totalWeight;
            base.cameraPosition.y = posY / totalWeight;
            base.cameraPosition.z = posZ / totalWeight;
            base.cameraQuaternion = normalizeAnimationCaptureQuaternion({
                x: quatX / totalWeight,
                y: quatY / totalWeight,
                z: quatZ / totalWeight,
                w: quatW / totalWeight
            });
        }

        out.push(base);
    }

    return out;
}

export function setAnimationTimelineCaptureCameraSmoothingEnabled(enabled) {
    animationTimelineCaptureState.cameraSmoothingEnabled = !!enabled;
    try {
        localStorage.setItem(
            'black-hole.anim-capture.camera-smoothing',
            animationTimelineCaptureState.cameraSmoothingEnabled ? '1' : '0'
        );
    } catch (e) {}
    updateAnimationTimelineCaptureUi();
}

function setAnimationTimelineCaptureFeedback(mode, text, tone) {
    if (!animationTimelineCaptureState.feedback[mode]) return;
    animationTimelineCaptureState.feedback[mode].text = text || 'Idle';
    animationTimelineCaptureState.feedback[mode].tone = tone || '';
}

export function updateAnimationTimelineCaptureUi() {
    function syncMode(mode, btnId, statusId) {
        var btn = document.getElementById(btnId);
        var status = document.getElementById(statusId);
        var smoothToggle = document.getElementById(mode + '-capture-smooth');
        var isActive = animationTimelineCaptureState.active &&
            animationTimelineCaptureState.mode === mode;
        var otherActive = animationTimelineCaptureState.active &&
            animationTimelineCaptureState.mode !== mode;
        var feedback = animationTimelineCaptureState.feedback[mode] ||
            { text: 'Idle', tone: '' };

        if (btn) {
            btn.classList.toggle('is-recording', isActive);
            btn.disabled = otherActive;
            btn.innerHTML = isActive
                ? '&#9632; STOP &amp; SAVE'
                : '&#9679; RECORD TO TIMELINE';
        }
        if (status) {
            var text = feedback.text || 'Idle';
            var tone = feedback.tone || '';
            if (isActive) {
                text = 'Recording ' +
                    animationTimelineCaptureState.lastElapsed.toFixed(2) + 's';
                tone = 'is-recording';
            } else if (otherActive) {
                text = 'Other capture active';
                tone = 'is-warning';
            }
            status.textContent = text;
            status.className = 'anim-capture-status' + (tone ? ' ' + tone : '');
        }
        if (smoothToggle) {
            smoothToggle.checked = !!animationTimelineCaptureState.cameraSmoothingEnabled;
        }
    }

    syncMode('dive', 'dive-capture-btn', 'dive-capture-status');
    syncMode('hover', 'hover-capture-btn', 'hover-capture-status');
}

export function readAnimationCaptureAnchorPosition(rawPosition) {
    if (!rawPosition || typeof rawPosition !== 'object') return null;
    var x = parseFloat(rawPosition.x);
    var y = parseFloat(rawPosition.y);
    var z = parseFloat(rawPosition.z);
    if (!isFinite(x) || !isFinite(y) || !isFinite(z)) return null;
    var out = new THREE.Vector3(x, y, z);
    return out.lengthSq() > 1e-10 ? out : null;
}

export function readAnimationCaptureAnchorVelocity(rawVelocity) {
    if (!rawVelocity || typeof rawVelocity !== 'object') return null;
    var x = parseFloat(rawVelocity.x);
    var y = parseFloat(rawVelocity.y);
    var z = parseFloat(rawVelocity.z);
    if (!isFinite(x) || !isFinite(y) || !isFinite(z)) return null;
    return new THREE.Vector3(x, y, z);
}

function buildAnimationTimelineCaptureSample(mode, elapsedSeconds) {
    var radius = mode === 'dive' ? diveState.currentR : hoverState.currentR;
    return {
        t: Math.max(0, elapsedSeconds),
        radius: radius,
        observerTime: (observer && typeof observer.time === 'number') ? observer.time : 0,
        cameraPanX: cameraPan ? cameraPan.x : 0,
        cameraPanY: cameraPan ? cameraPan.y : 0,
        cameraPosition: cloneVector3Plain(camera && camera.position ? camera.position : null),
        cameraQuaternion: cloneQuaternionPlain(camera && camera.quaternion ? camera.quaternion : null)
    };
}

function animationTimelineCaptureSamplesEqual(a, b) {
    if (!a || !b) return false;
    return Math.abs(a.radius - b.radius) < 1e-5 &&
        Math.abs(a.observerTime - b.observerTime) < 1e-5 &&
        Math.abs(a.cameraPanX - b.cameraPanX) < 1e-5 &&
        Math.abs(a.cameraPanY - b.cameraPanY) < 1e-5 &&
        Math.abs(a.cameraPosition.x - b.cameraPosition.x) < 1e-5 &&
        Math.abs(a.cameraPosition.y - b.cameraPosition.y) < 1e-5 &&
        Math.abs(a.cameraPosition.z - b.cameraPosition.z) < 1e-5 &&
        Math.abs(a.cameraQuaternion.x - b.cameraQuaternion.x) < 1e-5 &&
        Math.abs(a.cameraQuaternion.y - b.cameraQuaternion.y) < 1e-5 &&
        Math.abs(a.cameraQuaternion.z - b.cameraQuaternion.z) < 1e-5 &&
        Math.abs(a.cameraQuaternion.w - b.cameraQuaternion.w) < 1e-5;
}

function pushAnimationTimelineCaptureSample(mode, elapsedSeconds, force) {
    if (!animationTimelineCaptureState.active ||
        animationTimelineCaptureState.mode !== mode) {
        return false;
    }

    var sample = buildAnimationTimelineCaptureSample(mode, elapsedSeconds);
    var samples = animationTimelineCaptureState.samples;
    var last = samples.length ? samples[samples.length - 1] : null;

    if (last && Math.abs(last.t - sample.t) < 1e-4) {
        samples[samples.length - 1] = sample;
        animationTimelineCaptureState.lastElapsed = sample.t;
        animationTimelineCaptureState.nextSampleTime =
            sample.t + animationTimelineCaptureState.sampleInterval;
        return true;
    }
    if (!force && last && animationTimelineCaptureSamplesEqual(last, sample)) {
        return false;
    }

    samples.push(sample);
    animationTimelineCaptureState.lastElapsed = sample.t;
    animationTimelineCaptureState.nextSampleTime =
        sample.t + animationTimelineCaptureState.sampleInterval;
    return true;
}

export function finalizeAnimationTimelineCapture(mode) {
    if (!animationTimelineCaptureState.active ||
        animationTimelineCaptureState.mode !== mode) {
        return false;
    }

    var elapsed = Math.max(
        animationTimelineCaptureState.lastElapsed,
        (performance.now() - animationTimelineCaptureState.startedAtMs) / 1000.0
    );
    pushAnimationTimelineCaptureSample(mode, elapsed, true);
    var captureSamples = animationTimelineCaptureState.samples.slice();
    if (animationTimelineCaptureState.cameraSmoothingEnabled) {
        captureSamples = smoothAnimationTimelineCaptureSamples(captureSamples, 2);
    }

    var payload = {
        mode: mode,
        duration: animationTimelineCaptureState.lastElapsed,
        startPosition: cloneVector3Plain(animationTimelineCaptureState.startPosition),
        startVelocity: cloneVector3Plain(animationTimelineCaptureState.startVelocity),
        prevMotionState: !!animationTimelineCaptureState.prevMotionState,
        prevDistance: animationTimelineCaptureState.prevDistance,
        startObserverTime: animationTimelineCaptureState.startObserverTime,
        samples: captureSamples
    };

    animationTimelineCaptureState.active = false;
    animationTimelineCaptureState.mode = '';
    animationTimelineCaptureState.startedAtMs = 0;
    animationTimelineCaptureState.lastElapsed = 0;
    animationTimelineCaptureState.nextSampleTime = 0;
    animationTimelineCaptureState.samples = [];
    animationTimelineCaptureState.startPosition = null;
    animationTimelineCaptureState.startVelocity = null;
    animationTimelineCaptureState.startObserverTime = 0.0;

    var timelineBinding = getTimelinePanelBinding();
    var result = null;
    if (timelineBinding &&
        typeof timelineBinding.insertAnimationCapture === 'function') {
        result = timelineBinding.insertAnimationCapture(payload);
    }

    if (result && result.ok) {
        setAnimationTimelineCaptureFeedback(
            mode,
            'Saved ' + result.sampleCount + ' samples @ t=' +
                result.startTime.toFixed(2) + 's' +
                (animationTimelineCaptureState.cameraSmoothingEnabled ? ' (smoothed)' : ''),
            'is-ready'
        );
    } else {
        setAnimationTimelineCaptureFeedback(
            mode,
            (result && result.error) ? result.error : 'Timeline unavailable',
            'is-warning'
        );
    }

    updateAnimationTimelineCaptureUi();
    return !!(result && result.ok);
}

export function startAnimationTimelineCapture(mode) {
    if (mode !== 'dive' && mode !== 'hover') return false;
    if (animationTimelineCaptureState.active) {
        if (animationTimelineCaptureState.mode === mode) {
            return finalizeAnimationTimelineCapture(mode);
        }
        setAnimationTimelineCaptureFeedback(
            mode,
            'Stop the current capture first',
            'is-warning'
        );
        updateAnimationTimelineCaptureUi();
        return false;
    }

    var presentationRuntimeState = (typeof getPresentationState === 'function')
        ? getPresentationState()
        : null;
    if (presentationRuntimeState && presentationRuntimeState.recording) {
        setAnimationTimelineCaptureFeedback(
            mode,
            'Stop timeline recording first',
            'is-warning'
        );
        updateAnimationTimelineCaptureUi();
        return false;
    }
    if (presentationRuntimeState && presentationRuntimeState.playing &&
        typeof pausePresentation === 'function') {
        pausePresentation();
    }

    var modeWasAlreadyActive = (mode === 'dive')
        ? (diveState.active && !diveState.reachedSingularity)
        : hoverState.active;

    if (mode === 'dive') {
        if (!diveState.active || diveState.reachedSingularity) {
            startDive({ restart: true });
        } else if (diveState.paused) {
            diveState.timelineDriven = false;
            diveState.paused = false;
            updateDiveUI();
        }
        if (!diveState.active) {
            setAnimationTimelineCaptureFeedback(
                mode,
                'Unable to start live dive',
                'is-warning'
            );
            updateAnimationTimelineCaptureUi();
            return false;
        }
    } else {
        if (!hoverState.active) {
            startHover({ restart: true });
        } else if (hoverState.paused) {
            hoverState.timelineDriven = false;
            hoverState.paused = false;
            updateHoverUI();
        }
        if (!hoverState.active) {
            setAnimationTimelineCaptureFeedback(
                mode,
                'Unable to start live hover',
                'is-warning'
            );
            updateAnimationTimelineCaptureUi();
            return false;
        }
    }

    animationTimelineCaptureState.active = true;
    animationTimelineCaptureState.mode = mode;
    animationTimelineCaptureState.startedAtMs = performance.now();
    animationTimelineCaptureState.lastElapsed = 0.0;
    animationTimelineCaptureState.nextSampleTime = 0.0;
    animationTimelineCaptureState.sampleInterval = 1.0 / 30.0;
    animationTimelineCaptureState.samples = [];
    if (mode === 'dive' && !modeWasAlreadyActive) {
        animationTimelineCaptureState.startPosition = diveState.startPosition.clone();
        animationTimelineCaptureState.startVelocity = diveState.startVelocity.clone();
        animationTimelineCaptureState.prevMotionState = diveState.prevMotionState;
        animationTimelineCaptureState.prevDistance = diveState.prevDistance;
    } else if (mode === 'hover' && !modeWasAlreadyActive) {
        animationTimelineCaptureState.startPosition = hoverState.startPosition.clone();
        animationTimelineCaptureState.startVelocity = hoverState.startVelocity.clone();
        animationTimelineCaptureState.prevMotionState = hoverState.prevMotionState;
        animationTimelineCaptureState.prevDistance = hoverState.prevDistance;
    } else {
        animationTimelineCaptureState.startPosition = observer.position.clone();
        animationTimelineCaptureState.startVelocity = observer.velocity.clone();
        animationTimelineCaptureState.prevMotionState = shader.parameters.observer.motion;
        animationTimelineCaptureState.prevDistance = shader.parameters.observer.distance;
    }
    animationTimelineCaptureState.startObserverTime = observer.time;

    setAnimationTimelineCaptureFeedback(mode, 'Recording 0.00s', 'is-recording');
    pushAnimationTimelineCaptureSample(mode, 0.0, true);
    updateAnimationTimelineCaptureUi();
    return true;
}

export function toggleAnimationTimelineCapture(mode) {
    if (animationTimelineCaptureState.active &&
        animationTimelineCaptureState.mode === mode) {
        return finalizeAnimationTimelineCapture(mode);
    }
    return startAnimationTimelineCapture(mode);
}

export function updateAnimationTimelineCaptureFrame() {
    if (!animationTimelineCaptureState.active) return;

    var mode = animationTimelineCaptureState.mode;
    if (mode === 'dive' && !diveState.active && !diveState.reachedSingularity) {
        finalizeAnimationTimelineCapture(mode);
        return;
    }
    if (mode === 'hover' && !hoverState.active) {
        finalizeAnimationTimelineCapture(mode);
        return;
    }

    var elapsed = Math.max(
        0.0,
        (performance.now() - animationTimelineCaptureState.startedAtMs) / 1000.0
    );
    if (elapsed + 1e-6 >= animationTimelineCaptureState.nextSampleTime) {
        pushAnimationTimelineCaptureSample(mode, elapsed, false);
        updateAnimationTimelineCaptureUi();
    }
}


