import {
    scene,
    shader,
    applyRenderScaleFromSettings
} from '../../core/runtime/runtime-state.js';
import { applyQualityPresetValues } from '../../ui/quality-presets.js';
import { refreshPresentationUiBindings } from './path-bindings.js';

export function capturePresentationQualitySnapshot() {
    if (!shader || !shader.parameters) return null;
    var p = shader.parameters;
    return {
        quality: p.quality,
        n_steps: p.n_steps,
        sample_count: p.sample_count,
        max_revolutions: p.max_revolutions,
        rk4_integration: p.rk4_integration,
        cinematic_tonemap: p.cinematic_tonemap,
        resolution_scale: p.resolution_scale,
        taa_enabled: p.taa_enabled,
        taa: {
            history_weight: p.taa.history_weight,
            clip_box: p.taa.clip_box,
            motion_rejection: p.taa.motion_rejection,
            max_camera_delta: p.taa.max_camera_delta,
            motion_clip_scale: p.taa.motion_clip_scale
        }
    };
}

export function restorePresentationQualitySnapshot(snapshot) {
    if (!snapshot || !shader || !shader.parameters) return false;
    var p = shader.parameters;
    p.quality = snapshot.quality;
    p.n_steps = snapshot.n_steps;
    p.sample_count = snapshot.sample_count;
    p.max_revolutions = snapshot.max_revolutions;
    p.rk4_integration = snapshot.rk4_integration;
    p.cinematic_tonemap = snapshot.cinematic_tonemap;
    p.resolution_scale = snapshot.resolution_scale;
    p.taa_enabled = snapshot.taa_enabled;
    p.taa.history_weight = snapshot.taa.history_weight;
    p.taa.clip_box = snapshot.taa.clip_box;
    p.taa.motion_rejection = snapshot.taa.motion_rejection;
    p.taa.max_camera_delta = snapshot.taa.max_camera_delta;
    p.taa.motion_clip_scale = snapshot.taa.motion_clip_scale;

    if (typeof applyRenderScaleFromSettings === 'function') {
        applyRenderScaleFromSettings();
    }
    if (scene && typeof scene.updateShader === 'function') {
        scene.updateShader();
    }
    refreshPresentationUiBindings();
    return true;
}

export function applyPresentationRecordingQualityPreset(presetName) {
    if (!presetName || presetName === 'current') return false;
    if (!shader || !shader.parameters) return false;
    if (typeof applyQualityPresetValues !== 'function') return false;

    var preset = applyQualityPresetValues(shader.parameters, presetName);
    if (!preset) return false;

    if (typeof applyRenderScaleFromSettings === 'function') {
        applyRenderScaleFromSettings();
    }
    if (scene && typeof scene.updateShader === 'function') {
        scene.updateShader();
    }
    refreshPresentationUiBindings();
    return true;
}
