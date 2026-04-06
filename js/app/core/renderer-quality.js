"use strict";

import { $ } from '../vendor.js';
import {
    shader,
    scene,
    applyRenderScaleFromSettings
} from './runtime-state.js';
import { refreshRendererUiBindings } from './ui-bindings.js';
import { QUALITY_PRESETS, applyQualityPresetValues } from '../ui/quality-presets.js';

var QUALITY_BENCHMARK_STORAGE_KEY = 'black-hole-quality-benchmark-v4';
var QUALITY_BENCHMARK_SCHEMA_VERSION = 4;
var QUALITY_BENCHMARK_TARGET_FRAME_MS = 32.0;
var qualityBenchmarkState = null;

export function readStoredQualityPreset() {
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

export function storeQualityPreset(qualityName, avgFrameMs) {
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
    } catch (err) {}
}

export function applyQualityPresetRuntime(qualityName) {
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
    refreshRendererUiBindings();
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

export function beginQualityBenchmarkIfNeeded() {
    if (qualityBenchmarkState) return;
    if (readStoredQualityPreset()) return;

    qualityBenchmarkState = {
        phase: 'optimal',
        warmupFrames: 24,
        sampleFrames: 72,
        frameCount: 0,
        sampleCount: 0,
        accumulatedDt: 0.0
    };
    applyQualityPresetRuntime('optimal');
}

export function advanceQualityBenchmark(frameDt) {
    if (!qualityBenchmarkState) return;

    qualityBenchmarkState.frameCount++;
    if (qualityBenchmarkState.frameCount <= qualityBenchmarkState.warmupFrames) return;

    qualityBenchmarkState.accumulatedDt += frameDt;
    qualityBenchmarkState.sampleCount++;

    if (qualityBenchmarkState.sampleCount < qualityBenchmarkState.sampleFrames) return;

    var avgFrameMs =
        (qualityBenchmarkState.accumulatedDt / qualityBenchmarkState.sampleCount) * 1000.0;

    if (qualityBenchmarkState.phase === 'optimal') {
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
    }
}

export function isLikelyMobileDevice() {
    var ua = (window.navigator && window.navigator.userAgent) || '';
    var uaMobile = /(Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini|Mobile)/i
        .test(ua);
    var coarsePointer = !!(
        window.matchMedia &&
        window.matchMedia('(pointer: coarse)').matches
    );
    var smallViewport = Math.min(window.innerWidth || 0, window.innerHeight || 0) <= 900;
    return uaMobile || (coarsePointer && smallViewport);
}

export function clampResolutionScale(value) {
    return Math.max(0.35, Math.min(2.0, value || 1.0));
}
