import { renderer } from '../../core/runtime/runtime-state.js';
import { getBlackHoleRuntimeApi } from '../../core/runtime/runtime-registry.js';
import { QUALITY_PRESETS } from '../../ui/quality-presets.js';

export function choosePresentationMimeType() {
    if (typeof MediaRecorder === 'undefined') return '';
    if (typeof MediaRecorder.isTypeSupported !== 'function') {
        return 'video/webm';
    }

    var candidates = [
        'video/webm;codecs=vp9',
        'video/webm;codecs=vp8',
        'video/webm'
    ];
    for (var i = 0; i < candidates.length; i++) {
        if (MediaRecorder.isTypeSupported(candidates[i])) return candidates[i];
    }
    return '';
}

function presentationCaptureFilename(prefix, mimeType) {
    var now = new Date();

    function pad2(v) {
        return (v < 10 ? '0' : '') + v;
    }

    var stamp = now.getFullYear().toString() +
        pad2(now.getMonth() + 1) +
        pad2(now.getDate()) + '-' +
        pad2(now.getHours()) +
        pad2(now.getMinutes()) +
        pad2(now.getSeconds());
    var ext = 'webm';
    if (typeof mimeType === 'string' && mimeType) {
        var cleanMime = mimeType.toLowerCase();
        if (cleanMime.indexOf('png') !== -1) {
            ext = 'png';
        } else if (cleanMime.indexOf('jpeg') !== -1 || cleanMime.indexOf('jpg') !== -1) {
            ext = 'jpg';
        } else if (cleanMime.indexOf('mp4') !== -1) {
            ext = 'mp4';
        } else if (cleanMime.indexOf('webm') !== -1) {
            ext = 'webm';
        }
    }
    return (prefix || 'black-hole-presentation') + '-' + stamp + '.' + ext;
}

export function normalizePresentationRecordingQualityPreset(value) {
    if (typeof value !== 'string') return 'current';
    var clean = value.trim().toLowerCase();
    if (!clean || clean === 'current') return 'current';
    if (QUALITY_PRESETS && QUALITY_PRESETS[clean]) {
        return clean;
    }
    return 'current';
}

export function normalizePresentationRecordingMode(value) {
    if (typeof value !== 'string') return 'offline';
    var clean = value.trim().toLowerCase();
    if (clean === 'realtime' || clean === 'screen' || clean === 'live') {
        return 'realtime';
    }
    return 'offline';
}

export function normalizePresentationRecordingResolutionPreset(value) {
    if (typeof value !== 'string') return 'current';
    var clean = value.trim().toLowerCase();
    if (!clean || clean === 'current') return 'current';
    var match = /^(\d{3,5})x(\d{3,5})$/.exec(clean);
    if (!match) return 'current';

    var w = parseInt(match[1], 10);
    var h = parseInt(match[2], 10);
    if (!isFinite(w) || !isFinite(h)) return 'current';
    if (w < 160 || h < 90 || w > 8192 || h > 8192) return 'current';
    return w + 'x' + h;
}

export function resolvePresentationRecordingResolution(preset) {
    var normalized = normalizePresentationRecordingResolutionPreset(preset);
    if (normalized === 'current') {
        var currentWidth = Math.max(1, renderer && renderer.domElement ? (renderer.domElement.width || 1) : 1);
        var currentHeight = Math.max(1, renderer && renderer.domElement ? (renderer.domElement.height || 1) : 1);
        if ((currentWidth % 2) !== 0 && currentWidth > 2) currentWidth -= 1;
        if ((currentHeight % 2) !== 0 && currentHeight > 2) currentHeight -= 1;
        return {
            preset: 'current',
            width: currentWidth,
            height: currentHeight
        };
    }

    var match = /^(\d{3,5})x(\d{3,5})$/.exec(normalized);
    var width = match ? parseInt(match[1], 10) : 1920;
    var height = match ? parseInt(match[2], 10) : 1080;
    width = Math.max(160, Math.min(8192, width));
    height = Math.max(90, Math.min(8192, height));

    if ((width % 2) !== 0) width -= 1;
    if ((height % 2) !== 0) height -= 1;
    width = Math.max(160, width);
    height = Math.max(90, height);

    return {
        preset: normalized,
        width: width,
        height: height
    };
}

export function getPresentationRendererRuntimeApi() {
    var runtimeApi = null;
    if (typeof getBlackHoleRuntimeApi === 'function') {
        runtimeApi = getBlackHoleRuntimeApi('renderer');
    }
    if (!runtimeApi && typeof window !== 'undefined') {
        runtimeApi = window.blackHoleRendererRuntime;
    }
    if (!runtimeApi) return null;
    if (typeof runtimeApi.setOfflineSteppingActive !== 'function') return null;
    if (typeof runtimeApi.stepOfflineFrame !== 'function') return null;
    return runtimeApi;
}

export function isRendererContextLost() {
    var runtimeApi = getPresentationRendererRuntimeApi();
    return runtimeApi && typeof runtimeApi.isContextLost === 'function' && runtimeApi.isContextLost();
}

export function getPresentationWebMMuxerApi() {
    if (typeof window === 'undefined' || !window.WebMMuxer) return null;
    var muxApi = window.WebMMuxer;
    if (typeof muxApi.Muxer !== 'function') return null;
    if (typeof muxApi.ArrayBufferTarget !== 'function') return null;
    return muxApi;
}

export function getOfflinePresentationRecordingSupportState() {
    if (!renderer || !renderer.domElement) {
        return { supported: false, reason: 'Renderer not initialized yet.' };
    }
    if (typeof VideoFrame === 'undefined') {
        return { supported: false, reason: 'VideoFrame API is unavailable.' };
    }
    if (typeof VideoEncoder === 'undefined') {
        return { supported: false, reason: 'VideoEncoder API is unavailable.' };
    }
    if (!getPresentationWebMMuxerApi()) {
        return { supported: false, reason: 'WebM muxer library is unavailable.' };
    }
    if (!getPresentationRendererRuntimeApi()) {
        return { supported: false, reason: 'Renderer offline stepping API is unavailable.' };
    }
    return { supported: true, reason: '' };
}

export function isOfflinePresentationRecordingSupported() {
    return getOfflinePresentationRecordingSupportState().supported;
}

export function stopPresentationCaptureStreamTracks(stream) {
    if (!stream || typeof stream.getTracks !== 'function') return;
    var tracks = stream.getTracks();
    for (var i = 0; i < tracks.length; i++) {
        if (tracks[i] && typeof tracks[i].stop === 'function') {
            tracks[i].stop();
        }
    }
}

export function downloadPresentationCaptureDataUrl(dataUrl, mime, filenamePrefix) {
    if (!dataUrl || typeof dataUrl !== 'string' || typeof document === 'undefined') return false;
    var a = document.createElement('a');
    a.href = dataUrl;
    a.download = presentationCaptureFilename(filenamePrefix, mime || 'image/png');
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    return true;
}

export function downloadPresentationRecordingBlob(blob, mime, filenamePrefix) {
    if (!blob || blob.size <= 0 || typeof document === 'undefined') return;
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url;
    a.download = presentationCaptureFilename(filenamePrefix, mime);
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function() { URL.revokeObjectURL(url); }, 1500);
}
