// Role: Presentation timeline and recording controller.
//       Provides presets, keyframe evaluation, scripted events, annotations,
//       and realtime/offline capture hooks for slideshow-ready sequences.

import { $ } from '../../vendor.js';
import {
    camera,
    shader,
    scene,
    renderer,
    cameraControls
} from '../../core/runtime/runtime-state.js';
import { diveState, hoverState } from '../../core/scenarios/scenario-state.js';
import { initializeCamera } from '../../scene/camera.js';
import { startDive, resetDive } from '../../core/scenarios/dive.js';
import { startHover, resetHover } from '../../core/scenarios/hover.js';
import {
    PRESENTATION_PRESETS,
    PRESENTATION_PRESET_ORDER,
    presentationPresetLoadState,
    registerPresentationPreset,
    ensurePresentationPresetsLoaded
} from '../presets/presentation-presets.js';
import {
    presentationAnnotationState,
    presentationParamHudState,
    configurePresentationOverlay,
    syncPresentationTimelineUiConfig,
    normalizePresentationAnnotationsConfig,
    normalizePresentationParamHudConfig,
    applyPresentationTimelineUiConfig,
    setPresentationAnnotation,
    clearPresentationAnnotation,
    setPresentationAnnotationsEnabled,
    setPresentationAnnotationsIncludedInRecording,
    getPresentationAnnotationsState,
    setPresentationParamHudEnabled,
    setPresentationParamHudIncludedInRecording,
    getPresentationParamHudState,
    setParamHudLayout,
    isParamInHud,
    addParamToHud,
    removeParamFromHud,
    toggleParamInHud,
    clearParamHud,
    ensurePresentationAnnotationCanvas,
    updatePresentationOverlay
} from './presentation-overlay.js';
import {
    presentationEasing,
    samplePresentationTrack
} from './track-sampling.js';
import {
    resolvePresentationPath,
    presentationPathNeedsCompile,
    refreshPresentationUiBindings,
    setPresentationInteractionLock,
    setPresentationPathValue,
    getPresentationPathValue
} from './path-bindings.js';
import {
    choosePresentationMimeType,
    normalizePresentationRecordingQualityPreset,
    normalizePresentationRecordingMode,
    normalizePresentationRecordingResolutionPreset,
    resolvePresentationRecordingResolution,
    getPresentationRendererRuntimeApi,
    isRendererContextLost,
    getPresentationWebMMuxerApi,
    getOfflinePresentationRecordingSupportState,
    isOfflinePresentationRecordingSupported,
    stopPresentationCaptureStreamTracks,
    downloadPresentationCaptureDataUrl,
    downloadPresentationRecordingBlob
} from './presentation-recording-support.js';
import {
    capturePresentationQualitySnapshot,
    restorePresentationQualitySnapshot,
    applyPresentationRecordingQualityPreset
} from './presentation-recording-quality.js';
export {
    PRESENTATION_PRESETS,
    PRESENTATION_PRESET_ORDER,
    registerPresentationPreset,
    ensurePresentationPresetsLoaded,
    setPresentationAnnotationsEnabled,
    setPresentationAnnotationsIncludedInRecording,
    getPresentationAnnotationsState,
    setPresentationParamHudEnabled,
    setPresentationParamHudIncludedInRecording,
    getPresentationParamHudState,
    setParamHudLayout,
    isParamInHud,
    addParamToHud,
    removeParamFromHud,
    toggleParamInHud,
    clearParamHud,
    updatePresentationOverlay,
    getPresentationPathValue
};

// â”€â”€â”€ Presentation Timeline + Capture â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
// Timeline system for scripted camera/parameter animations suitable for slides
// and demos. Supports keyframed parameter tracks, timed events (dive/hover),
// and optional canvas recording via MediaRecorder.
var presentationState = {
    active: false,
    paused: true,
    loop: false,
    time: 0.0,
    duration: 0.0,
    timeline: null,
    eventCursor: 0,
    compileRequested: false
};

var presentationCaptureState = {
    active: false,
    recorder: null,
    stream: null,
    chunks: [],
    mode: 'realtime',
    preferredMode: 'offline',
    fps: 60,
    bitrateMbps: 20.0,
    filenamePrefix: 'black-hole-presentation',
    autoStopOnPresentationEnd: true,
    mimeType: '',
    qualityPreset: 'current',
    resolutionPreset: 'current',
    outputWidth: 0,
    outputHeight: 0,
    backgroundThrottleDetected: false,
    restoreQualitySnapshot: null,
    includeAnnotationsInRecording: false,
    captureCanvas: null,
    captureCtx: null,
    compositeCanvas: null,
    compositeCtx: null,
    compositeRaf: 0,
    offlineJob: null,
    offlineUnavailableReason: '',
    rendererResizedForRecording: false
};

var presentationUiRefreshAccumulator = 0.0;

function clonePresentationData(value) {
    return JSON.parse(JSON.stringify(value));
}

function presentationClamp(v, lo, hi) {
    return Math.max(lo, Math.min(hi, v));
}

function syncPresentationTimelineUiConfigToState() {
    if (!presentationState.timeline) return;

    var annotations = getPresentationAnnotationsState();
    var paramHud = getPresentationParamHudState();

    presentationState.timeline.loop = !!presentationState.loop;
    presentationState.timeline.annotations = {
        enabled: !!annotations.enabled,
        includeInRecording: !!annotations.includeInRecording
    };
    presentationState.timeline.paramHud = {
        enabled: !!paramHud.enabled,
        includeInRecording: !!paramHud.includeInRecording,
        anchorX: paramHud.anchorX,
        anchorY: paramHud.anchorY,
        fontSize: paramHud.fontSize,
        items: clonePresentationData(paramHud.items)
    };
}

configurePresentationOverlay({
    onUiConfigChange: syncPresentationTimelineUiConfigToState,
    getPathValue: getPresentationPathValue
});


function normalizePresentationTimeline(timeline) {
    if (!timeline) return null;

    var raw = clonePresentationData(timeline);
    var out = {
        name: raw.name || 'Custom',
        loop: !!raw.loop,
        duration: 0.0,
        tracks: [],
        events: [],
        annotationTracks: [],
        annotations: normalizePresentationAnnotationsConfig(raw.annotations),
        paramHud: normalizePresentationParamHudConfig(raw.paramHud)
    };

    if (Array.isArray(raw.annotationTracks)) {
        for (var ai = 0; ai < raw.annotationTracks.length; ai++) {
            var at = raw.annotationTracks[ai];
            if (at && typeof at === 'object') out.annotationTracks.push({ label: at.label || ('Annotation ' + (ai + 1)) });
        }
    }

    var maxTime = 0.0;
    // Paths that are UI/performance meta-settings and must never be driven by
    // a timeline (they would silently override the user's chosen setting on play).
    // These are all the parameters owned by the quality preset system.
    var TIMELINE_EXCLUDED_PATHS = {
        'quality': true,
        'n_steps': true,
        'sample_count': true,
        'max_revolutions': true,
        'rk4_integration': true,
        'cinematic_tonemap': true,
        'resolution_scale': true,
        'taa_enabled': true,
        'taa.history_weight': true,
        'taa.clip_box': true,
        'taa.motion_rejection': true,
        'taa.max_camera_delta': true,
        'taa.motion_clip_scale': true
    };
    function normalizeExcludedTimelinePath(path) {
        if (typeof path !== 'string') return '';
        var clean = path.trim();
        if (clean.indexOf('params.') === 0) clean = clean.substring('params.'.length);
        if (clean.indexOf('shader.parameters.') === 0) {
            clean = clean.substring('shader.parameters.'.length);
        }
        return clean;
    }

    var tracks = Array.isArray(raw.tracks) ? raw.tracks : [];
    for (var i = 0; i < tracks.length; i++) {
        var track = tracks[i];
        if (!track || typeof track.path !== 'string') continue;
        if (TIMELINE_EXCLUDED_PATHS[normalizeExcludedTimelinePath(track.path)]) continue;

        var keys = Array.isArray(track.keys) ? track.keys : [];
        var normalizedKeys = [];
        for (var k = 0; k < keys.length; k++) {
            var key = keys[k];
            if (!key) continue;
            var t = parseFloat(key.t);
            if (!isFinite(t)) continue;
            t = Math.max(0.0, t);
            if (t > maxTime) maxTime = t;
            normalizedKeys.push({
                t: t,
                v: key.v,
                ease: key.ease || 'linear'
            });
        }
        if (!normalizedKeys.length) continue;
        normalizedKeys.sort(function(a, b) { return a.t - b.t; });
        out.tracks.push({
            path: track.path,
            compile: !!track.compile,
            keys: normalizedKeys
        });
    }

    var events = Array.isArray(raw.events) ? raw.events : [];
    for (var j = 0; j < events.length; j++) {
        var ev = events[j];
        if (!ev || typeof ev.action !== 'string') continue;
        if (ev.action === 'set' &&
            typeof ev.path === 'string' &&
            TIMELINE_EXCLUDED_PATHS[normalizeExcludedTimelinePath(ev.path)]) {
            continue;
        }
        var et = parseFloat(ev.t);
        if (!isFinite(et)) continue;
        et = Math.max(0.0, et);
        if (et > maxTime) maxTime = et;
        var normEv = {
            t: et,
            action: ev.action,
            path: ev.path,
            value: ev.value,
            compile: !!ev.compile,
            note: ev.note
        };
        if (typeof ev.channel === 'number') normEv.channel = ev.channel;
        if (ev._pairOf !== undefined) normEv._pairOf = ev._pairOf;
        var extraKeys = Object.keys(ev);
        for (var ek = 0; ek < extraKeys.length; ek++) {
            var extraKey = extraKeys[ek];
            if (extraKey === 't' || extraKey === 'action' ||
                extraKey === 'path' || extraKey === 'value' ||
                extraKey === 'compile' || extraKey === 'note' ||
                extraKey === 'channel' || extraKey === '_pairOf') {
                continue;
            }
            if (ev[extraKey] !== undefined) {
                normEv[extraKey] = clonePresentationData(ev[extraKey]);
            }
        }
        out.events.push(normEv);
    }
    out.events.sort(function(a, b) { return a.t - b.t; });

    var duration = parseFloat(raw.duration);
    if (!isFinite(duration) || duration <= 0.0) {
        duration = maxTime;
    }
    out.duration = Math.max(duration, maxTime, 0.001);
    return out;
}

export function presentationTimelineHasTrack(path) {
    if (!presentationState.timeline || typeof path !== 'string') return false;
    var clean = path.trim();
    if (!clean) return false;

    var tracks = presentationState.timeline.tracks;
    for (var i = 0; i < tracks.length; i++) {
        if (!tracks[i] || typeof tracks[i].path !== 'string') continue;
        if (tracks[i].path.trim() === clean) return true;
    }
    return false;
}

function flushPresentationShaderCompile() {
    if (!presentationState.compileRequested) return;
    if (scene && typeof scene.updateShader === 'function') {
        scene.updateShader();
    }
    presentationState.compileRequested = false;
}

function applyPresentationTracks(timeSeconds) {
    if (!presentationState.timeline) return;

    var compileNeeded = false;
    var anyChanged = false;
    var tracks = presentationState.timeline.tracks;
    for (var i = 0; i < tracks.length; i++) {
        var track = tracks[i];
        var sampledValue = samplePresentationTrack(track, timeSeconds);
        var changed = setPresentationPathValue(track.path, sampledValue);
        if (changed) {
            anyChanged = true;
            if (track.compile || presentationPathNeedsCompile(track.path)) {
                compileNeeded = true;
            }
        }
    }

    if (compileNeeded && scene && typeof scene.updateShader === 'function') {
        scene.updateShader();
    }

    if (anyChanged) shader.needsUpdate = true;
}

function executePresentationEvent(event) {
    if (!event || !event.action) return;

    switch (event.action) {
        case 'set':
            if (typeof event.path === 'string') {
                var changed = setPresentationPathValue(event.path, event.value);
                if (changed &&
                    (event.compile || presentationPathNeedsCompile(event.path))) {
                    presentationState.compileRequested = true;
                }
                if (changed) refreshPresentationUiBindings();
            }
            break;
        case 'startDive':
            startDive({
                restart: true,
                anchorPosition: event.position,
                anchorVelocity: event.velocity,
                prevMotionState: event.prevMotionState,
                prevDistance: event.prevDistance,
                observerTime: event.observerTime
            });
            break;
        case 'resetDive':
            resetDive();
            break;
        case 'pauseDive':
            if (diveState.active) {
                diveState.paused = true;
                diveState.timelineDriven = false;
                updateDiveUI();
            }
            break;
        case 'startHover':
            startHover({
                restart: true,
                anchorPosition: event.position,
                anchorVelocity: event.velocity,
                prevMotionState: event.prevMotionState,
                prevDistance: event.prevDistance,
                observerTime: event.observerTime
            });
            break;
        case 'resetHover':
            resetHover();
            break;
        case 'pauseHover':
            if (hoverState.active) {
                hoverState.paused = true;
                hoverState.timelineDriven = false;
                updateHoverUI();
            }
            break;
        case 'updateShader':
            presentationState.compileRequested = true;
            break;
        case 'annotation':
            setPresentationAnnotation(event.note, event.channel || 0);
            break;
        case 'clearAnnotation':
            clearPresentationAnnotation(typeof event.channel === 'number' ? event.channel : undefined);
            break;
    }
}

function processPresentationEvents(fromTime, toTime) {
    if (!presentationState.timeline) return;
    var events = presentationState.timeline.events;
    while (presentationState.eventCursor < events.length &&
        events[presentationState.eventCursor].t <= toTime + 1e-6) {
        var eventTime = events[presentationState.eventCursor].t;
        if (eventTime > fromTime + 1e-6) {
            executePresentationEvent(events[presentationState.eventCursor]);
        }
        presentationState.eventCursor++;
    }
    flushPresentationShaderCompile();
}

export function setPresentationTimeline(timeline) {
    var normalized = normalizePresentationTimeline(timeline);
    if (!normalized) return false;

    presentationState.timeline = normalized;
    presentationState.loop = normalized.loop;
    presentationState.duration = normalized.duration;
    presentationState.time = 0.0;
    presentationState.eventCursor = 0;
    presentationState.compileRequested = false;
    presentationState.active = false;
    presentationState.paused = true;
    presentationUiRefreshAccumulator = 0.0;
    setPresentationInteractionLock(false);

    if (diveState.active || diveState.reachedSingularity) resetDive();
    if (hoverState.active) resetHover();
    clearPresentationAnnotation();
    applyPresentationTimelineUiConfig(normalized);

    applyPresentationTracks(0.0);
    shader.needsUpdate = true;

    // Notify bottom timeline panel (if open) to resync
    try {
        window.dispatchEvent(new CustomEvent('presentation:timeline-panel-sync'));
    } catch (e) {}

    return true;
}

export function listPresentationPresets() {
    return PRESENTATION_PRESET_ORDER.slice();
}

export function getPresentationTimeline() {
    if (!presentationState.timeline) return null;
    syncPresentationTimelineUiConfig();
    return clonePresentationData(presentationState.timeline);
}

export function loadPresentationPreset(name) {
    if (!PRESENTATION_PRESETS[name]) {
        ensurePresentationPresetsLoaded();
        return false;
    }
    var preset = clonePresentationData(PRESENTATION_PRESETS[name]);
    if (!preset.name) preset.name = name;
    return setPresentationTimeline(preset);
}

export function seekPresentation(timeSeconds) {
    if (!presentationState.timeline) return false;

    var t = parseFloat(timeSeconds);
    if (!isFinite(t)) return false;
    t = Math.max(0.0, Math.min(presentationState.duration, t));

    presentationState.time = t;
    presentationState.eventCursor = 0;

    clearPresentationAnnotation();
    if (diveState.active || diveState.reachedSingularity) resetDive();
    if (hoverState.active) resetHover();

    applyPresentationTracks(0.0);
    processPresentationEvents(-1.0, t);
    applyPresentationTracks(t);

    shader.needsUpdate = true;
    refreshPresentationUiBindings();
    return true;
}

export function playPresentation(fromStart) {
    if (!presentationState.timeline) return false;

    var shouldRestart = !!fromStart ||
        presentationState.time >= presentationState.duration - 1e-6;
    if (shouldRestart) {
        if (typeof initializeCamera === 'function' && camera) {
            initializeCamera(camera);
            if (typeof cameraControls !== 'undefined' && cameraControls && cameraControls.target) {
                cameraControls.target.set(0, 0, 0);
            }
        }
        if (diveState.active || diveState.reachedSingularity) resetDive();
        if (hoverState.active) resetHover();
        seekPresentation(0.0);
        presentationState.eventCursor = 0;
        processPresentationEvents(-1.0, 0.0); // fire t=0 events once
        applyPresentationTracks(0.0);
    } else if (presentationState.time <= 1e-6) {
        presentationState.eventCursor = 0;
        processPresentationEvents(-1.0, 0.0);
        applyPresentationTracks(0.0);
    }

    presentationState.active = true;
    presentationState.paused = false;
    setPresentationInteractionLock(true);
    shader.needsUpdate = true;
    return true;
}

export function pausePresentation() {
    if (!presentationState.timeline) return false;
    presentationState.active = false;
    presentationState.paused = true;
    setPresentationInteractionLock(false);
    return true;
}

export function stopPresentation() {
    if (!presentationState.timeline) return false;
    presentationState.active = false;
    presentationState.paused = true;
    setPresentationInteractionLock(false);
    if (diveState.active || diveState.reachedSingularity) resetDive();
    if (hoverState.active) resetHover();
    seekPresentation(0.0);
    clearPresentationAnnotation();
    refreshPresentationUiBindings();
    return true;
}

export function setPresentationLoop(enabled) {
    presentationState.loop = !!enabled;
    syncPresentationTimelineUiConfig();
    return presentationState.loop;
}

export function getPresentationBackgroundThrottleState() {
    var visible = true;
    var focused = true;

    if (typeof document !== 'undefined') {
        if (typeof document.visibilityState === 'string') {
            visible = (document.visibilityState === 'visible');
        } else if (typeof document.hidden === 'boolean') {
            visible = !document.hidden;
        }

        if (typeof document.hasFocus === 'function') {
            try {
                focused = !!document.hasFocus();
            } catch (err) {
                focused = true;
            }
        }
    }

    var throttleRisk = !visible;
    var reason = '';
    if (throttleRisk) {
        reason = 'Tab/window is hidden or minimized; browser may throttle or pause rendering.';
    }

    return {
        visible: visible,
        focused: focused,
        throttleRisk: throttleRisk,
        reason: reason
    };
}

export function getPresentationState() {
    var offlineJob = presentationCaptureState.offlineJob;
    var backgroundState = getPresentationBackgroundThrottleState();
    if (presentationCaptureState.active && backgroundState.throttleRisk) {
        presentationCaptureState.backgroundThrottleDetected = true;
    }
    var backgroundThrottleDetected =
        !!backgroundState.throttleRisk || !!presentationCaptureState.backgroundThrottleDetected;
    var backgroundThrottleReason = backgroundState.reason || '';
    if (!backgroundThrottleReason && presentationCaptureState.backgroundThrottleDetected) {
        backgroundThrottleReason = 'Background throttling was detected during this recording.';
    }
    var offlineFramesDone = offlineJob ? (offlineJob.frameCount || 0) : 0;
    var offlineFramesTotal = offlineJob ? (offlineJob.totalFrames || 0) : 0;
    var offlineElapsedSeconds = 0;
    if (offlineJob && offlineJob.wallStartMs) {
        offlineElapsedSeconds = Math.max(0.0, (Date.now() - offlineJob.wallStartMs) / 1000.0);
    }
    var offlineSinceLastFrameSeconds = 0;
    if (offlineJob && offlineJob.lastFrameWallMs) {
        offlineSinceLastFrameSeconds = Math.max(0.0, (Date.now() - offlineJob.lastFrameWallMs) / 1000.0);
    }
    var offlineRenderFps = (offlineElapsedSeconds > 0.0)
        ? (offlineFramesDone / offlineElapsedSeconds)
        : 0.0;
    var offlineProgress = (offlineFramesTotal > 0)
        ? presentationClamp(offlineFramesDone / offlineFramesTotal, 0.0, 1.0)
        : 0.0;
    var offlineEtaSeconds = -1;
    if (offlineFramesTotal > 0 && offlineRenderFps > 1e-6) {
        offlineEtaSeconds = Math.max(0.0, (offlineFramesTotal - offlineFramesDone) / offlineRenderFps);
    }
    var offlinePhase = (offlineJob && offlineJob.phase) ? offlineJob.phase : 'idle';
    var offlineFinalizingProgress = -1;
    if (offlineJob && offlinePhase === 'finalizing-encode') {
        var q = (offlineJob.encoder && typeof offlineJob.encoder.encodeQueueSize === 'number')
            ? offlineJob.encoder.encodeQueueSize
            : 0;
        var q0 = Math.max(1, offlineJob.finalizingStartQueue || q || 1);
        offlineFinalizingProgress = 0.1 + 0.7 * (1.0 - presentationClamp(q / q0, 0.0, 1.0));
    } else if (offlineJob && offlinePhase === 'finalizing-mux') {
        offlineFinalizingProgress = 0.9;
    } else if (offlineJob && offlinePhase === 'finalizing-download') {
        offlineFinalizingProgress = 0.97;
    } else if (offlineJob && offlinePhase === 'done') {
        offlineFinalizingProgress = 1.0;
    }

    return {
        loaded: !!presentationState.timeline,
        name: presentationState.timeline ? presentationState.timeline.name : '',
        presets_loaded: !!presentationPresetLoadState.loaded,
        presets_loading: !!presentationPresetLoadState.loading,
        playing: presentationState.active && !presentationState.paused,
        drives_observer_time: presentationTimelineHasTrack('observerState.time'),
        loop: !!presentationState.loop,
        time: presentationState.time,
        duration: presentationState.duration,
        recording: !!presentationCaptureState.active,
        recording_mode: presentationCaptureState.mode || 'realtime',
        recording_mode_preferred: presentationCaptureState.preferredMode || 'offline',
        recording_offline_supported: isOfflinePresentationRecordingSupported(),
        recording_offline_unavailable_reason: presentationCaptureState.offlineUnavailableReason || '',
        recording_background_visible: !!backgroundState.visible,
        recording_background_focused: !!backgroundState.focused,
        recording_background_throttle_risk: !!backgroundState.throttleRisk,
        recording_background_throttle_detected: backgroundThrottleDetected,
        recording_background_throttle_reason: backgroundThrottleReason,
        recording_offline_phase: offlinePhase,
        recording_offline_frames_done: offlineFramesDone,
        recording_offline_frames_total: offlineFramesTotal,
        recording_offline_progress: offlineProgress,
        recording_offline_finalizing_progress: offlineFinalizingProgress,
        recording_offline_elapsed_s: offlineElapsedSeconds,
        recording_offline_since_last_frame_s: offlineSinceLastFrameSeconds,
        recording_offline_render_fps: offlineRenderFps,
        recording_offline_eta_s: offlineEtaSeconds,
        recording_offline_encode_queue: (offlineJob && offlineJob.encoder &&
            typeof offlineJob.encoder.encodeQueueSize === 'number')
            ? offlineJob.encoder.encodeQueueSize
            : 0,
        recording_offline_timeline_done_s: offlineFramesDone / Math.max(presentationCaptureState.fps || 60, 1),
        recording_offline_timeline_total_s: (offlineFramesTotal > 0)
            ? (offlineFramesTotal / Math.max(presentationCaptureState.fps || 60, 1))
            : 0,
        recording_quality_preset: presentationCaptureState.qualityPreset || 'current',
        recording_resolution_preset: presentationCaptureState.resolutionPreset || 'current',
        recording_output_width: presentationCaptureState.outputWidth || 0,
        recording_output_height: presentationCaptureState.outputHeight || 0,
        annotations_enabled: !!presentationAnnotationState.enabled,
        annotations_in_recording: !!presentationAnnotationState.includeInRecording,
        param_hud_enabled: !!presentationParamHudState.enabled,
        param_hud_in_recording: !!presentationParamHudState.includeInRecording,
        param_hud_count: presentationParamHudState.items.length
    };
}

export function updatePresentation(dt) {
    if (!presentationState.timeline || !presentationState.active || presentationState.paused) return;

    var previousTime = presentationState.time;
    var nextTime = previousTime + dt;
    var duration = Math.max(presentationState.duration, 0.001);

    if (nextTime < duration) {
        processPresentationEvents(previousTime, nextTime);
        presentationState.time = nextTime;
        applyPresentationTracks(nextTime);
        if (presentationParamHudState.enabled && presentationParamHudState.items.length > 0) {
            updatePresentationOverlay();
        }
        presentationUiRefreshAccumulator += dt;
        if (presentationUiRefreshAccumulator >= 0.2) {
            presentationUiRefreshAccumulator = 0.0;
            refreshPresentationUiBindings();
        }
        return;
    }

    // Final frame of the segment
    processPresentationEvents(previousTime, duration);
    applyPresentationTracks(duration);
    if (presentationParamHudState.enabled && presentationParamHudState.items.length > 0) {
        updatePresentationOverlay();
    }

    if (presentationState.loop) {
        nextTime = nextTime % duration;
        presentationState.time = nextTime;
        presentationState.eventCursor = 0;
        processPresentationEvents(-1.0, nextTime);
        applyPresentationTracks(nextTime);
        refreshPresentationUiBindings();
        return;
    }

    presentationState.time = duration;
    presentationState.active = false;
    presentationState.paused = true;
    setPresentationInteractionLock(false);
    refreshPresentationUiBindings();

    if (presentationCaptureState.active && presentationCaptureState.autoStopOnPresentationEnd) {
        stopPresentationRecording();
    }
}

function setPresentationRendererOfflineStepping(enabled) {
    var runtimeApi = getPresentationRendererRuntimeApi();
    if (!runtimeApi) return false;
    runtimeApi.setOfflineSteppingActive(!!enabled);
    return true;
}

function stepPresentationRendererOfflineFrame(dt) {
    var runtimeApi = getPresentationRendererRuntimeApi();
    if (!runtimeApi) return false;
    runtimeApi.stepOfflineFrame(dt);
    return true;
}

function syncPresentationCompositeSize(canvas) {
    if (!canvas || !renderer || !renderer.domElement) return;
    canvas.width = Math.max(1, renderer.domElement.width || 1);
    canvas.height = Math.max(1, renderer.domElement.height || 1);
}

function ensurePresentationCaptureCanvas(width, height) {
    if (!presentationCaptureState.captureCanvas || !presentationCaptureState.captureCtx) {
        presentationCaptureState.captureCanvas = document.createElement('canvas');
        presentationCaptureState.captureCtx = presentationCaptureState.captureCanvas.getContext('2d');
    }

    if (!presentationCaptureState.captureCtx) return false;
    presentationCaptureState.captureCanvas.width = Math.max(1, Math.floor(width || 1));
    presentationCaptureState.captureCanvas.height = Math.max(1, Math.floor(height || 1));
    presentationCaptureState.outputWidth = presentationCaptureState.captureCanvas.width;
    presentationCaptureState.outputHeight = presentationCaptureState.captureCanvas.height;
    return true;
}

function syncPresentationCaptureCanvasForCurrentResolution() {
    if (!presentationCaptureState.captureCanvas || !renderer || !renderer.domElement) return;
    if (presentationCaptureState.resolutionPreset !== 'current') return;
    if (presentationCaptureState.active && presentationCaptureState.mode === 'offline') return;

    var width = Math.max(1, renderer.domElement.width || 1);
    var height = Math.max(1, renderer.domElement.height || 1);
    if (presentationCaptureState.captureCanvas.width !== width ||
        presentationCaptureState.captureCanvas.height !== height) {
        presentationCaptureState.captureCanvas.width = width;
        presentationCaptureState.captureCanvas.height = height;
    }
    presentationCaptureState.outputWidth = presentationCaptureState.captureCanvas.width;
    presentationCaptureState.outputHeight = presentationCaptureState.captureCanvas.height;
}

function drawPresentationCaptureFrame() {
    if (!presentationCaptureState.captureCanvas || !presentationCaptureState.captureCtx) return false;
    if (!renderer || !renderer.domElement) return false;

    syncPresentationCaptureCanvasForCurrentResolution();

    var source = renderer.domElement;
    var includeOverlayInRecording = presentationCaptureState.includeAnnotationsInRecording ||
        (presentationParamHudState.includeInRecording && presentationParamHudState.items.length > 0);
    if (includeOverlayInRecording) {
        if (!drawPresentationCompositeFrame(
            presentationCaptureState.compositeCanvas,
            presentationCaptureState.compositeCtx
        )) {
            return false;
        }
        source = presentationCaptureState.compositeCanvas;
    }

    var targetCanvas = presentationCaptureState.captureCanvas;
    var ctx = presentationCaptureState.captureCtx;
    var w = targetCanvas.width;
    var h = targetCanvas.height;
    ctx.clearRect(0, 0, w, h);

    var sourceWidth = Math.max(1,
        source.width || source.videoWidth || source.naturalWidth || source.clientWidth || w);
    var sourceHeight = Math.max(1,
        source.height || source.videoHeight || source.naturalHeight || source.clientHeight || h);

    var sourceAspect = sourceWidth / sourceHeight;
    var targetAspect = w / h;
    var epsilon = 1e-6;

    if (Math.abs(sourceAspect - targetAspect) <= epsilon) {
        ctx.drawImage(source, 0, 0, sourceWidth, sourceHeight, 0, 0, w, h);
        return true;
    }

    // Keep source aspect ratio when recording to a different output aspect.
    // This prevents visible flattening/stretching when, for example, recording
    // 16:10 viewport content to a 16:9 file.
    var drawWidth = w;
    var drawHeight = h;
    var drawX = 0;
    var drawY = 0;

    if (sourceAspect > targetAspect) {
        drawWidth = w;
        drawHeight = Math.max(1, Math.round(w / sourceAspect));
        drawY = Math.round((h - drawHeight) * 0.5);
    } else {
        drawHeight = h;
        drawWidth = Math.max(1, Math.round(h * sourceAspect));
        drawX = Math.round((w - drawWidth) * 0.5);
    }

    ctx.fillStyle = '#000';
    ctx.fillRect(0, 0, w, h);
    ctx.drawImage(source, 0, 0, sourceWidth, sourceHeight, drawX, drawY, drawWidth, drawHeight);
    return true;
}

function drawPresentationCompositeFrame(canvas, ctx) {
    if (!canvas || !ctx || !renderer || !renderer.domElement) return false;

    if (canvas.width !== renderer.domElement.width ||
        canvas.height !== renderer.domElement.height) {
        syncPresentationCompositeSize(canvas);
    }

    if (typeof updatePresentationOverlay === 'function') {
        updatePresentationOverlay();
    }

    var w = canvas.width;
    var h = canvas.height;
    ctx.clearRect(0, 0, w, h);
    ctx.drawImage(renderer.domElement, 0, 0, w, h);
    var showOverlayCanvas = (presentationAnnotationState.enabled ||
        (presentationParamHudState.enabled && presentationParamHudState.items.length > 0));
    if (showOverlayCanvas && presentationAnnotationState.canvas) {
        var annCanvas = presentationAnnotationState.canvas;
        var annCtx = presentationAnnotationState.ctx;
        // When the annotation canvas is a different resolution from the composite
        // (e.g. recording at 2560Ã—1440 while the window is a different aspect ratio)
        // temporarily resize it to the composite dimensions so that text is laid out
        // at the correct positions and scale, then restore it to the window size.
        if (annCtx && (annCanvas.width !== w || annCanvas.height !== h)) {
            var savedStyleW = annCanvas.style.width;
            var savedStyleH = annCanvas.style.height;
            var savedW = annCanvas.width;
            var savedH = annCanvas.height;
            annCanvas.style.width = w + 'px';
            annCanvas.style.height = h + 'px';
            annCanvas.width = w;
            annCanvas.height = h;
            annCtx.setTransform(1, 0, 0, 1, 0, 0);
            if (typeof updatePresentationOverlay === 'function') {
                updatePresentationOverlay();
            }
            ctx.drawImage(annCanvas, 0, 0, w, h);
            // Restore annotation canvas to window dimensions for the DOM overlay.
            annCanvas.style.width = savedStyleW;
            annCanvas.style.height = savedStyleH;
            annCanvas.width = savedW;
            annCanvas.height = savedH;
            var dpr = Math.max((typeof window !== 'undefined' && window.devicePixelRatio) || 1, 1);
            annCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
            if (typeof updatePresentationOverlay === 'function') {
                updatePresentationOverlay();
            }
        } else {
            ctx.drawImage(annCanvas, 0, 0, w, h);
        }
    }
    return true;
}

function stopPresentationCompositeCapture() {
    if (presentationCaptureState.compositeRaf) {
        cancelAnimationFrame(presentationCaptureState.compositeRaf);
    }
    presentationCaptureState.compositeRaf = 0;
    presentationCaptureState.compositeCanvas = null;
    presentationCaptureState.compositeCtx = null;
}

function clearPresentationCaptureBuffers() {
    stopPresentationCompositeCapture();
    presentationCaptureState.includeAnnotationsInRecording = false;
    presentationCaptureState.captureCanvas = null;
    presentationCaptureState.captureCtx = null;
    presentationCaptureState.outputWidth = 0;
    presentationCaptureState.outputHeight = 0;
}

function cleanupPresentationRecordingState() {
    stopPresentationCaptureStreamTracks(presentationCaptureState.stream);

    setPresentationRendererOfflineStepping(false);

    if (presentationCaptureState.rendererResizedForRecording) {
        var restoreApi = getPresentationRendererRuntimeApi();
        if (restoreApi && typeof restoreApi.restoreWindowSizeAfterRecording === 'function') {
            restoreApi.restoreWindowSizeAfterRecording();
        }
        presentationCaptureState.rendererResizedForRecording = false;
    }

    presentationCaptureState.active = false;
    presentationCaptureState.recorder = null;
    presentationCaptureState.stream = null;
    presentationCaptureState.chunks = [];
    clearPresentationCaptureBuffers();
    presentationCaptureState.backgroundThrottleDetected = false;
    presentationCaptureState.offlineJob = null;

    restorePresentationQualitySnapshot(presentationCaptureState.restoreQualitySnapshot);
    presentationCaptureState.restoreQualitySnapshot = null;
}

function finalizeOfflinePresentationRecording() {
    var offlineJob = presentationCaptureState.offlineJob;
    if (!offlineJob || offlineJob.finalizing) return;
    offlineJob.finalizing = true;
    if (offlineJob.failed) {
        offlineJob.phase = 'failed';
    } else {
        offlineJob.phase = 'finalizing-encode';
        offlineJob.finalizingStartQueue =
            (offlineJob.encoder && typeof offlineJob.encoder.encodeQueueSize === 'number')
                ? offlineJob.encoder.encodeQueueSize
                : 0;
    }
    refreshPresentationUiBindings();

    function finishCleanup() {
        if (offlineJob && !offlineJob.failed) {
            offlineJob.phase = 'done';
            refreshPresentationUiBindings();
        }
        cleanupPresentationRecordingState();
        refreshPresentationUiBindings();
    }

    if (offlineJob.failed) {
        if (offlineJob.encoder && typeof offlineJob.encoder.close === 'function') {
            try { offlineJob.encoder.close(); } catch (closeErr) {}
        }
        finishCleanup();
        return;
    }

    function finalizeMuxedOutput() {
        if (!offlineJob.encoder || typeof offlineJob.encoder.flush !== 'function') {
            finishCleanup();
            return;
        }

        // Guard against two failure modes caused by GPU context loss / hardware encoder crash:
        //   1. flush() throws synchronously (encoder state already 'closed')
        //   2. flush() returns a Promise that never settles (HW encoder died silently,
        //      no error callback fired). Without a watchdog this hangs the UI at 10% forever.
        var flushSettled = false;
        var flushWatchdog = setTimeout(function() {
            if (flushSettled) return;
            flushSettled = true;
            offlineJob.flushWatchdog = null;
            console.warn('Offline encoder flush timed out â€” GPU context may have been lost.');
            offlineJob.failed = true;
            offlineJob.phase = 'failed';
            refreshPresentationUiBindings();
            if (offlineJob.encoder && typeof offlineJob.encoder.close === 'function') {
                try { offlineJob.encoder.close(); } catch (closeErr) {}
            }
            finishCleanup();
        }, 30000);

        var flushPromise;
        try {
            offlineJob.flushWatchdog = flushWatchdog;
            flushPromise = offlineJob.encoder.flush();
        } catch (syncFlushErr) {
            clearTimeout(flushWatchdog);
            offlineJob.flushWatchdog = null;
            flushSettled = true;
            console.warn('Offline encoder flush threw synchronously:', syncFlushErr);
            offlineJob.failed = true;
            offlineJob.phase = 'failed';
            refreshPresentationUiBindings();
            if (offlineJob.encoder && typeof offlineJob.encoder.close === 'function') {
                try { offlineJob.encoder.close(); } catch (closeErr) {}
            }
            finishCleanup();
            return;
        }

        Promise.resolve(flushPromise).then(function() {
            clearTimeout(flushWatchdog);
            offlineJob.flushWatchdog = null;
            if (flushSettled) return;
            flushSettled = true;
            offlineJob.phase = 'finalizing-mux';
            refreshPresentationUiBindings();

            // Yield one tick so the UI can paint the finalization phase
            setTimeout(function() {
                try {
                    if (offlineJob.muxer && typeof offlineJob.muxer.finalize === 'function') {
                        offlineJob.muxer.finalize();
                    }
                    offlineJob.phase = 'finalizing-download';
                    refreshPresentationUiBindings();

                    setTimeout(function() {
                        if (offlineJob.encoder && typeof offlineJob.encoder.close === 'function') {
                            try { offlineJob.encoder.close(); } catch (closeErr) {}
                        }

                        if (offlineJob.writableFileStream) {
                            // FileSystemWritableFileStreamTarget path: data already streamed to
                            // disk during encoding. Just close the writable stream (which also
                            // flushes any pending buffered writes) and we're done â€” no Blob ever
                            // accumulates in RAM.
                            offlineJob.writableFileStream.close().then(function() {
                                finishCleanup();
                            }).catch(function(closeErr) {
                                console.warn('Failed to close recording file stream:', closeErr);
                                offlineJob.failed = true;
                                finishCleanup();
                            });
                        } else {
                            // ArrayBufferTarget path: build a Blob from the in-memory buffer
                            // and trigger the browser download dialog.
                            try {
                                var buffer = offlineJob.target && offlineJob.target.buffer;
                                if (buffer) {
                                    var mime = 'video/webm';
                                    var blob = new Blob([buffer], { type: mime });
                                    downloadPresentationRecordingBlob(
                                        blob,
                                        mime,
                                        presentationCaptureState.filenamePrefix
                                    );
                                }
                            } catch (downloadErr) {
                                console.warn('Offline recording download preparation failed:', downloadErr);
                                offlineJob.failed = true;
                            }
                            finishCleanup();
                        }
                    }, 0);
                } catch (muxErr) {
                    console.warn('Offline recording finalize failed:', muxErr);
                    offlineJob.failed = true;
                    if (offlineJob.encoder && typeof offlineJob.encoder.close === 'function') {
                        try { offlineJob.encoder.close(); } catch (closeErr2) {}
                    }
                    finishCleanup();
                }
            }, 0);
        }).catch(function(flushErr) {
            clearTimeout(flushWatchdog);
            offlineJob.flushWatchdog = null;
            if (flushSettled) return;
            flushSettled = true;
            console.warn('Offline encoder flush failed:', flushErr);
            offlineJob.failed = true;
            offlineJob.phase = 'failed';
            refreshPresentationUiBindings();
            if (offlineJob.encoder && typeof offlineJob.encoder.close === 'function') {
                try { offlineJob.encoder.close(); } catch (closeErr3) {}
            }
            finishCleanup();
        });
    }

    if (offlineJob.encoder && typeof offlineJob.encoder.flush === 'function') {
        finalizeMuxedOutput();
        return;
    }

    finishCleanup();
}

function runOfflinePresentationRecordingLoop() {
    var offlineJob = presentationCaptureState.offlineJob;
    if (!presentationCaptureState.active || !offlineJob) return;

    var frameDt = 1.0 / Math.max(presentationCaptureState.fps, 1);

    function encodeNextFrame() {
        if (!presentationCaptureState.active || !presentationCaptureState.offlineJob ||
            presentationCaptureState.offlineJob !== offlineJob) {
            return;
        }

        if (offlineJob.stopRequested || offlineJob.failed) {
            finalizeOfflinePresentationRecording();
            return;
        }

        if (isRendererContextLost()) {
            console.warn('Offline recording aborted: WebGL context lost (GPU driver reset / TDR).');
            offlineJob.failed = true;
            finalizeOfflinePresentationRecording();
            return;
        }

        if (!stepPresentationRendererOfflineFrame(frameDt)) {
            console.warn('Offline recording aborted: renderer stepping API unavailable.');
            offlineJob.failed = true;
            finalizeOfflinePresentationRecording();
            return;
        }

        if (!drawPresentationCaptureFrame()) {
            console.warn('Offline recording aborted: failed to draw capture frame.');
            offlineJob.failed = true;
            finalizeOfflinePresentationRecording();
            return;
        }

        if (!offlineJob.encoder || offlineJob.encoder.state === 'closed') {
            offlineJob.failed = true;
            finalizeOfflinePresentationRecording();
            return;
        }

        var frame = null;
        try {
            frame = new VideoFrame(presentationCaptureState.captureCanvas, {
                timestamp: offlineJob.nextTimestampUs,
                duration: offlineJob.frameDurationUs
            });
        } catch (err) {
            console.warn('Offline recording aborted: failed to create video frame.', err);
            offlineJob.failed = true;
            finalizeOfflinePresentationRecording();
            return;
        }

        try {
            var keyEvery = Math.max(1, Math.round(presentationCaptureState.fps));
            offlineJob.encoder.encode(frame, {
                keyFrame: (offlineJob.frameCount % keyEvery) === 0
            });
        } catch (encodeErr) {
            console.warn('Offline recording aborted: failed to encode frame.', encodeErr);
            frame.close();
            offlineJob.failed = true;
            finalizeOfflinePresentationRecording();
            return;
        }
        frame.close();

        offlineJob.frameCount += 1;
        offlineJob.nextTimestampUs += offlineJob.frameDurationUs;
        offlineJob.lastFrameWallMs = Date.now();

        if (offlineJob.totalFrames > 0 && offlineJob.frameCount >= offlineJob.totalFrames) {
            offlineJob.stopRequested = true;
            finalizeOfflinePresentationRecording();
            return;
        }

        var autoStopReachedEnd = false;
        if (offlineJob.totalFrames <= 0 &&
            presentationCaptureState.autoStopOnPresentationEnd &&
            typeof getPresentationState === 'function') {
            var state = getPresentationState();
            autoStopReachedEnd = !!state.loaded && !state.playing;
        }

        if (offlineJob.stopRequested || autoStopReachedEnd) {
            finalizeOfflinePresentationRecording();
            return;
        }

        // Always yield via setTimeout â€” never call encodeNextFrame() directly.
        // This prevents the JS thread from hammering the GPU in tight synchronous bursts,
        // which is the primary cause of VRAM spikes and TDR (GPU driver reset) on long renders.
        //
        // When the encode queue is backing up (hardware encoder falling behind), add a real
        // delay so the encoder has time to drain VRAM before the next frame is submitted.
        // Each 2K VideoFrame is ~14 MB of raw data; without back-pressure they accumulate fast.
        var queueSize = offlineJob.encoder ? offlineJob.encoder.encodeQueueSize : 0;
        var frameDelay = (queueSize >= 3) ? 50 : (queueSize >= 1) ? 8 : 0;
        setTimeout(encodeNextFrame, frameDelay);
    }

    encodeNextFrame();
}

export function capturePresentationScreenshot(options) {
    if (!renderer || !renderer.domElement) return false;
    if (presentationCaptureState.active) {
        presentationCaptureState.offlineUnavailableReason =
            'Offline screenshot capture is unavailable while recording is active.';
        refreshPresentationUiBindings();
        return false;
    }

    options = options || {};
    var qualityPreset = normalizePresentationRecordingQualityPreset(
        (options.qualityPreset === undefined) ? 'cinematic' : options.qualityPreset
    );
    if (qualityPreset === 'current') qualityPreset = 'cinematic';

    var resolutionPreset = normalizePresentationRecordingResolutionPreset(
        (options.recordingResolution === undefined)
            ? presentationCaptureState.resolutionPreset
            : options.recordingResolution
    );
    var resolvedResolution = resolvePresentationRecordingResolution(resolutionPreset);
    var includeAnnotations = (options.includeAnnotationsInScreenshot === undefined)
        ? !!presentationAnnotationState.includeInRecording
        : !!options.includeAnnotationsInScreenshot;
    var filenamePrefix = options.filenamePrefix || 'black-hole-screenshot';
    var previousCaptureQualityPreset = presentationCaptureState.qualityPreset;
    var previousCaptureResolutionPreset = presentationCaptureState.resolutionPreset;

    var previousQualitySnapshot = capturePresentationQualitySnapshot();
    var qualityOverridden = false;

    function fail(reason) {
        presentationCaptureState.offlineUnavailableReason =
            reason || 'Offline screenshot capture failed.';
        refreshPresentationUiBindings();
        return false;
    }

    if (!previousQualitySnapshot) {
        return fail('Offline screenshot capture failed: renderer quality state is unavailable.');
    }

    if (!applyPresentationRecordingQualityPreset(qualityPreset)) {
        return fail('Offline screenshot capture failed: unable to apply Offline quality preset.');
    }
    qualityOverridden = true;

    presentationCaptureState.qualityPreset = qualityPreset;
    presentationCaptureState.resolutionPreset = resolvedResolution.preset;
    presentationCaptureState.includeAnnotationsInRecording = includeAnnotations;
    presentationCaptureState.offlineUnavailableReason = '';

    try {
        if (includeAnnotations) {
            ensurePresentationAnnotationCanvas();
            updatePresentationOverlay();
            var overlayCompositeCanvas = document.createElement('canvas');
            var overlayCompositeCtx = overlayCompositeCanvas.getContext('2d');
            if (!overlayCompositeCtx) {
                return fail('Offline screenshot capture failed: annotation compositing is unavailable.');
            }
            syncPresentationCompositeSize(overlayCompositeCanvas);
            presentationCaptureState.compositeCanvas = overlayCompositeCanvas;
            presentationCaptureState.compositeCtx = overlayCompositeCtx;
        } else {
            stopPresentationCompositeCapture();
        }

        if (!ensurePresentationCaptureCanvas(resolvedResolution.width, resolvedResolution.height)) {
            return fail('Offline screenshot capture failed: capture canvas is unavailable.');
        }

        // Force a fresh frame render at cinematic quality before copying to PNG.
        if (typeof render === 'function') {
            render();
        } else if (!stepPresentationRendererOfflineFrame(1.0 / 60.0)) {
            return fail('Offline screenshot capture failed: renderer frame API is unavailable.');
        }

        if (!drawPresentationCaptureFrame()) {
            return fail('Offline screenshot capture failed: could not read rendered pixels.');
        }

        var captureCanvas = presentationCaptureState.captureCanvas;
        if (!captureCanvas) {
            return fail('Offline screenshot capture failed: capture canvas is unavailable.');
        }

        var mimeType = 'image/png';
        var downloadStarted = false;

        if (typeof captureCanvas.toBlob === 'function') {
            captureCanvas.toBlob(function(blob) {
                if (blob && blob.size > 0) {
                    downloadPresentationRecordingBlob(blob, mimeType, filenamePrefix);
                    return;
                }
                try {
                    var fallbackUrl = captureCanvas.toDataURL(mimeType);
                    downloadPresentationCaptureDataUrl(fallbackUrl, mimeType, filenamePrefix);
                } catch (fallbackErr) {
                    console.warn('Offline screenshot fallback download failed:', fallbackErr);
                }
            }, mimeType);
            downloadStarted = true;
        } else if (typeof captureCanvas.toDataURL === 'function') {
            var dataUrl = captureCanvas.toDataURL(mimeType);
            downloadStarted =
                downloadPresentationCaptureDataUrl(dataUrl, mimeType, filenamePrefix);
        }

        if (!downloadStarted) {
            return fail('Offline screenshot capture failed: PNG export is unsupported.');
        }

        presentationCaptureState.offlineUnavailableReason = '';
        refreshPresentationUiBindings();
        return true;
    } catch (err) {
        console.warn('Offline screenshot capture failed:', err);
        return fail('Offline screenshot capture failed.');
    } finally {
        clearPresentationCaptureBuffers();
        presentationCaptureState.qualityPreset = previousCaptureQualityPreset;
        presentationCaptureState.resolutionPreset = previousCaptureResolutionPreset;
        if (qualityOverridden && previousQualitySnapshot) {
            restorePresentationQualitySnapshot(previousQualitySnapshot);
        } else {
            refreshPresentationUiBindings();
        }
    }
}

export function stopPresentationRecording() {
    if (!presentationCaptureState.active) return false;

    if (presentationCaptureState.mode === 'offline' && presentationCaptureState.offlineJob) {
        presentationCaptureState.offlineJob.stopRequested = true;
        return true;
    }

    if (!presentationCaptureState.recorder) {
        cleanupPresentationRecordingState();
        return true;
    }

    if (presentationCaptureState.recorder.state !== 'inactive') {
        presentationCaptureState.recorder.stop();
    } else {
        cleanupPresentationRecordingState();
    }
    return true;
}

export function startPresentationRecording(options) {
    if (!renderer || !renderer.domElement) return false;
    if (presentationCaptureState.active) return false;

    options = options || {};
    var fps = parseFloat(options.fps);
    if (!isFinite(fps) || fps <= 0) fps = presentationCaptureState.fps;
    fps = Math.max(10, Math.min(120, fps));

    var bitrate = parseFloat(options.bitrateMbps);
    if (!isFinite(bitrate) || bitrate <= 0) bitrate = presentationCaptureState.bitrateMbps;
    bitrate = Math.max(2.0, Math.min(80.0, bitrate));

    var includeAnnotationsInRecording = (options.includeAnnotationsInRecording === undefined)
        ? presentationAnnotationState.includeInRecording
        : !!options.includeAnnotationsInRecording;
    var includeOverlayInRecording = includeAnnotationsInRecording ||
        (presentationParamHudState.includeInRecording && presentationParamHudState.items.length > 0);

    var requestedMode = normalizePresentationRecordingMode(
        (options.recordingMode === undefined)
            ? presentationCaptureState.preferredMode
            : options.recordingMode
    );
    presentationCaptureState.preferredMode = requestedMode;

    var recordingMode = requestedMode;
    presentationCaptureState.offlineUnavailableReason = '';

    if (recordingMode === 'offline') {
        var offlineSupportState = getOfflinePresentationRecordingSupportState();
        if (!offlineSupportState.supported) {
            presentationCaptureState.offlineUnavailableReason = offlineSupportState.reason;
            console.warn('Offline recording unavailable:', offlineSupportState.reason);
            refreshPresentationUiBindings();
            return false;
        }
        if (!presentationState.timeline) {
            presentationCaptureState.offlineUnavailableReason =
                'Offline recording requires a loaded presentation timeline.';
            console.warn('Offline recording unavailable:',
                presentationCaptureState.offlineUnavailableReason);
            refreshPresentationUiBindings();
            return false;
        }
    } else if (typeof MediaRecorder === 'undefined') {
        presentationCaptureState.offlineUnavailableReason =
            'Realtime capture requires the MediaRecorder API.';
        refreshPresentationUiBindings();
        return false;
    }

    var qualityPreset = normalizePresentationRecordingQualityPreset(
        (options.qualityPreset === undefined)
            ? presentationCaptureState.qualityPreset
            : options.qualityPreset
    );
    var resolutionPreset = normalizePresentationRecordingResolutionPreset(
        (options.recordingResolution === undefined)
            ? presentationCaptureState.resolutionPreset
            : options.recordingResolution
    );
    var previousQualitySnapshot = null;
    var qualityOverridden = false;
    if (qualityPreset !== 'current') {
        previousQualitySnapshot = capturePresentationQualitySnapshot();
        qualityOverridden = !!previousQualitySnapshot &&
            applyPresentationRecordingQualityPreset(qualityPreset);
        if (!qualityOverridden) {
            previousQualitySnapshot = null;
            qualityPreset = 'current';
        }
    }

    function rollbackRecordingQualityOverride() {
        if (qualityOverridden && previousQualitySnapshot) {
            restorePresentationQualitySnapshot(previousQualitySnapshot);
        }
    }

    var filenamePrefix = options.filenamePrefix || presentationCaptureState.filenamePrefix;
    var autoStop = (options.autoStopOnPresentationEnd === undefined)
        ? presentationCaptureState.autoStopOnPresentationEnd
        : !!options.autoStopOnPresentationEnd;

    var stream = null;
    var compositeTick = null;
    var recorder = null;
    var mimeType = 'video/webm';
    var offlineJob = null;
    if (includeOverlayInRecording) {
        ensurePresentationAnnotationCanvas();
        updatePresentationOverlay();

        var overlayCompositeCanvas = document.createElement('canvas');
        var overlayCompositeCtx = overlayCompositeCanvas.getContext('2d');
        if (!overlayCompositeCtx) {
            rollbackRecordingQualityOverride();
            return false;
        }
        syncPresentationCompositeSize(overlayCompositeCanvas);
        presentationCaptureState.compositeCanvas = overlayCompositeCanvas;
        presentationCaptureState.compositeCtx = overlayCompositeCtx;
    } else {
        stopPresentationCompositeCapture();
    }

    var resolvedResolution = resolvePresentationRecordingResolution(resolutionPreset);
    if (!ensurePresentationCaptureCanvas(resolvedResolution.width, resolvedResolution.height)) {
        rollbackRecordingQualityOverride();
        return false;
    }
    presentationCaptureState.resolutionPreset = resolvedResolution.preset;
    if (!drawPresentationCaptureFrame()) {
        rollbackRecordingQualityOverride();
        return false;
    }

    if (recordingMode === 'offline') {
        var width = Math.max(1, presentationCaptureState.captureCanvas.width || 1);
        var height = Math.max(1, presentationCaptureState.captureCanvas.height || 1);
        if (resolvedResolution.preset !== 'current') {
            var offlineRuntimeApi = getPresentationRendererRuntimeApi();
            if (offlineRuntimeApi && typeof offlineRuntimeApi.resizeForOfflineRecording === 'function') {
                offlineRuntimeApi.resizeForOfflineRecording(width, height);
                presentationCaptureState.rendererResizedForRecording = true;
            }
        }
        var muxApi = getPresentationWebMMuxerApi();
        var codecVariants = [
            { encoder: 'vp09.00.10.08', muxer: 'V_VP9' },
            { encoder: 'vp8', muxer: 'V_VP8' },
            { encoder: 'av01.0.08M.08', muxer: 'V_AV1' }
        ];

        // Phase 1: find a supported codec using lightweight test encoders.
        // We do this before creating the real muxer+target so that we can use
        // FileSystemWritableFileStreamTarget (which streams to disk) without
        // having to discard any partially-written file on codec fallback.
        var selectedVariant = null;
        for (var cv = 0; cv < codecVariants.length; cv++) {
            var testEncoder = null;
            try {
                testEncoder = new VideoEncoder({ output: function() {}, error: function() {} });
                testEncoder.configure({
                    codec: codecVariants[cv].encoder,
                    width: width,
                    height: height,
                    bitrate: Math.round(bitrate * 1000000.0),
                    framerate: fps
                });
                selectedVariant = codecVariants[cv];
                testEncoder.close();
                break;
            } catch (codecErr) {
                if (testEncoder && typeof testEncoder.close === 'function') {
                    try { testEncoder.close(); } catch (closeErr) {}
                }
            }
        }

        if (!selectedVariant) {
            presentationCaptureState.offlineUnavailableReason =
                'No supported WebCodecs video encoder was found for offline rendering.';
            rollbackRecordingQualityOverride();
            refreshPresentationUiBindings();
            return false;
        }

        // Phase 2: create the real target, muxer, and encoder with the winning codec.
        // Prefer FileSystemWritableFileStreamTarget (streams directly to disk, no giant
        // ArrayBuffer accumulating in RAM) when the caller passed a writable file stream.
        var writableFileStream = options.writableFileStream || null;
        var target = null;
        var muxer = null;
        var encoder = null;
        try {
            if (writableFileStream && typeof muxApi.FileSystemWritableFileStreamTarget === 'function') {
                target = new muxApi.FileSystemWritableFileStreamTarget(writableFileStream);
            } else {
                writableFileStream = null;
                target = new muxApi.ArrayBufferTarget();
            }
            muxer = new muxApi.Muxer({
                target: target,
                video: {
                    codec: selectedVariant.muxer,
                    width: width,
                    height: height,
                    frameRate: fps
                }
            });
            encoder = new VideoEncoder({
                output: function(chunk, meta) {
                    try {
                        if (offlineJob && offlineJob.muxer) {
                            offlineJob.muxer.addVideoChunk(chunk, meta);
                        }
                    } catch (muxErr) {
                        console.warn('Offline recording mux error:', muxErr);
                        if (offlineJob) offlineJob.failed = true;
                    }
                },
                error: function(err) {
                    console.warn('Offline recording encoder error:', err);
                    if (offlineJob) {
                        offlineJob.failed = true;
                        // If the error fires while flush() is already pending (e.g. after a
                        // GPU context loss), kick the watchdog so cleanup runs immediately
                        // rather than waiting the full timeout.
                        if (offlineJob.flushWatchdog) {
                            clearTimeout(offlineJob.flushWatchdog);
                            offlineJob.flushWatchdog = null;
                            offlineJob.phase = 'failed';
                            refreshPresentationUiBindings();
                            cleanupPresentationRecordingState();
                            refreshPresentationUiBindings();
                        }
                    }
                }
            });
            encoder.configure({
                codec: selectedVariant.encoder,
                width: width,
                height: height,
                bitrate: Math.round(bitrate * 1000000.0),
                framerate: fps
            });
        } catch (setupErr) {
            console.warn('Offline recording setup failed:', setupErr);
            if (encoder && typeof encoder.close === 'function') { try { encoder.close(); } catch(e) {} }
            presentationCaptureState.offlineUnavailableReason = 'Failed to initialize offline encoder.';
            rollbackRecordingQualityOverride();
            refreshPresentationUiBindings();
            return false;
        }

        offlineJob = {
            stopRequested: false,
            finalizing: false,
            failed: false,
            phase: 'preparing',
            encoder: encoder,
            muxer: muxer,
            target: target,
            writableFileStream: writableFileStream,
            codec: selectedVariant.encoder,
            frameDurationUs: Math.max(1, Math.round(1000000.0 / fps)),
            nextTimestampUs: 0,
            frameCount: 0,
            totalFrames: 0,
            totalTimelineSeconds: 0,
            wallStartMs: 0,
            lastFrameWallMs: 0,
            finalizingStartQueue: 0,
            compositeCanvas: presentationCaptureState.compositeCanvas,
            compositeCtx: presentationCaptureState.compositeCtx,
            captureCanvas: presentationCaptureState.captureCanvas
        };
    } else {
        if (!presentationCaptureState.captureCanvas ||
            typeof presentationCaptureState.captureCanvas.captureStream !== 'function') {
            rollbackRecordingQualityOverride();
            return false;
        }
        stream = presentationCaptureState.captureCanvas.captureStream(fps);
        compositeTick = function() {
            if (!presentationCaptureState.active) return;
            if (!drawPresentationCaptureFrame()) return;
            presentationCaptureState.compositeRaf = requestAnimationFrame(compositeTick);
        };

        mimeType = choosePresentationMimeType();
        var recorderConfig = {};
        if (mimeType) recorderConfig.mimeType = mimeType;
        recorderConfig.videoBitsPerSecond = Math.round(bitrate * 1000000.0);

        try {
            recorder = new MediaRecorder(stream, recorderConfig);
        } catch (err) {
            try {
                recorder = new MediaRecorder(stream);
                mimeType = recorder.mimeType || '';
            } catch (err2) {
                stopPresentationCaptureStreamTracks(stream);
                rollbackRecordingQualityOverride();
                return false;
            }
        }
        mimeType = mimeType || recorder.mimeType || 'video/webm';
    }

    presentationCaptureState.active = true;
    presentationCaptureState.recorder = recorder;
    presentationCaptureState.stream = stream;
    presentationCaptureState.chunks = [];
    presentationCaptureState.mode = recordingMode;
    presentationCaptureState.fps = fps;
    presentationCaptureState.bitrateMbps = bitrate;
    presentationCaptureState.filenamePrefix = filenamePrefix;
    presentationCaptureState.autoStopOnPresentationEnd = autoStop;
    presentationCaptureState.mimeType = mimeType;
    presentationCaptureState.qualityPreset = qualityPreset;
    presentationCaptureState.resolutionPreset = resolvedResolution.preset;
    presentationCaptureState.outputWidth = presentationCaptureState.captureCanvas
        ? presentationCaptureState.captureCanvas.width
        : resolvedResolution.width;
    presentationCaptureState.outputHeight = presentationCaptureState.captureCanvas
        ? presentationCaptureState.captureCanvas.height
        : resolvedResolution.height;
    presentationCaptureState.restoreQualitySnapshot =
        qualityOverridden ? previousQualitySnapshot : null;
    presentationCaptureState.includeAnnotationsInRecording = includeAnnotationsInRecording;
    presentationCaptureState.backgroundThrottleDetected = false;
    presentationCaptureState.offlineJob = offlineJob;

    if (recordingMode === 'offline') {
        setPresentationRendererOfflineStepping(true);
        if (typeof playPresentation === 'function') {
            var playState = getPresentationState();
            if (playState.loaded && !playState.playing) {
                playPresentation(false);
            }
        }
        if (offlineJob) {
            var startTimelineTime = presentationState.time;
            var totalDuration = Math.max(presentationState.duration, 0.0);
            var totalTimelineSeconds = 0.0;
            if (presentationCaptureState.autoStopOnPresentationEnd) {
                if (presentationState.loop) {
                    totalTimelineSeconds = totalDuration;
                } else {
                    totalTimelineSeconds = Math.max(0.0, totalDuration - startTimelineTime);
                }
                if (totalTimelineSeconds <= 1e-6) totalTimelineSeconds = 1.0 / Math.max(fps, 1);
                offlineJob.totalTimelineSeconds = totalTimelineSeconds;
                offlineJob.totalFrames = Math.max(1, Math.ceil(totalTimelineSeconds * fps));
            } else {
                offlineJob.totalTimelineSeconds = 0.0;
                offlineJob.totalFrames = 0;
            }
            offlineJob.wallStartMs = Date.now();
            offlineJob.lastFrameWallMs = offlineJob.wallStartMs;
            offlineJob.phase = 'rendering';
        }
        runOfflinePresentationRecordingLoop();
        shader.needsUpdate = true;
        refreshPresentationUiBindings();
        return true;
    }

    if (compositeTick) {
        presentationCaptureState.compositeRaf = requestAnimationFrame(compositeTick);
    }

    recorder.ondataavailable = function(event) {
        if (event && event.data && event.data.size > 0) {
            presentationCaptureState.chunks.push(event.data);
        }
    };

    recorder.onstop = function() {
        var finalMime = presentationCaptureState.mimeType || 'video/webm';
        if (presentationCaptureState.chunks.length > 0) {
            var blob = new Blob(presentationCaptureState.chunks, { type: finalMime });
            downloadPresentationRecordingBlob(
                blob,
                finalMime,
                presentationCaptureState.filenamePrefix
            );
        }
        cleanupPresentationRecordingState();
        refreshPresentationUiBindings();
    };

    recorder.onerror = function(err) {
        console.warn('Presentation recorder error:', err);
    };

    try {
        recorder.start(250);
    } catch (startErr) {
        console.warn('Presentation recorder failed to start:', startErr);
        cleanupPresentationRecordingState();
        return false;
    }

    shader.needsUpdate = true;
    refreshPresentationUiBindings();
    return true;
}

ensurePresentationPresetsLoaded();

if (typeof window !== 'undefined') {
    window.blackHolePresentation = {
        listPresets: listPresentationPresets,
        ensurePresetsLoaded: ensurePresentationPresetsLoaded,
        loadPreset: loadPresentationPreset,
        getTimeline: getPresentationTimeline,
        setTimeline: setPresentationTimeline,
        play: playPresentation,
        pause: pausePresentation,
        stop: stopPresentation,
        seek: seekPresentation,
        setLoop: setPresentationLoop,
        state: getPresentationState,
        hasTrack: presentationTimelineHasTrack,
        setAnnotationsEnabled: setPresentationAnnotationsEnabled,
        setAnnotationsIncludedInRecording: setPresentationAnnotationsIncludedInRecording,
        annotationState: getPresentationAnnotationsState,
        showAnnotation: setPresentationAnnotation,
        clearAnnotation: clearPresentationAnnotation,
        getPathValue: getPresentationPathValue,
        setParamHudEnabled: setPresentationParamHudEnabled,
        setParamHudIncludedInRecording: setPresentationParamHudIncludedInRecording,
        paramHudState: getPresentationParamHudState,
        setParamHudLayout: setParamHudLayout,
        addParamToHud: addParamToHud,
        removeParamFromHud: removeParamFromHud,
        toggleParamInHud: toggleParamInHud,
        isParamInHud: isParamInHud,
        clearParamHud: clearParamHud,
        startRecording: startPresentationRecording,
        stopRecording: stopPresentationRecording,
        captureScreenshot: capturePresentationScreenshot
    };
}




