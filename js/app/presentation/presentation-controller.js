// Role: Presentation timeline and recording controller.
//       Provides presets, keyframe evaluation, scripted events, annotations,
//       and realtime/offline capture hooks for slideshow-ready sequences.

import { $, THREE } from '../vendor.js';
import {
    camera,
    observer,
    shader,
    scene,
    renderer,
    cameraControls,
    cameraPan,
    distanceController,
    refreshAllControllersGlobal
} from '../core/runtime-state.js';
import { diveState, hoverState } from '../core/scenario-state.js';
import { initializeCamera, updateCamera } from '../scene/camera.js';
import { startDive, resetDive, seekDive } from '../core/dive.js';
import { startHover, resetHover, seekHover } from '../core/hover.js';
import { applyQualityPresetValues, QUALITY_PRESETS } from '../ui/quality-presets.js';
import { getBlackHoleRuntimeApi, getBlackHoleUiBinding } from '../core/runtime-registry.js';

// ─── Presentation Timeline + Capture ─────────────────────────────────────────
// Timeline system for scripted camera/parameter animations suitable for slides
// and demos. Supports keyframed parameter tracks, timed events (dive/hover),
// and optional canvas recording via MediaRecorder.
var PRESENTATION_PRESET_MANIFEST_PATH = 'js/app/presentation/presets/manifest.json';
export var PRESENTATION_PRESETS = {};
export var PRESENTATION_PRESET_ORDER = [];

var presentationPresetLoadState = {
    loading: false,
    loaded: false,
    error: null,
    promise: null
};

function requestPresentationJson(path) {
    return new Promise(function(resolve, reject) {
        if (typeof $ !== 'undefined' && $ && typeof $.getJSON === 'function') {
            $.getJSON(path)
                .done(function(data) { resolve(data); })
                .fail(function(jqXHR, textStatus, err) {
                    reject(new Error('Failed to load ' + path + ': ' + (err || textStatus || 'unknown error')));
                });
            return;
        }

        if (typeof XMLHttpRequest === "undefined") {
            reject(new Error('No JSON loader available for ' + path));
            return;
        }

        var req = new XMLHttpRequest();
        req.open('GET', path, true);
        req.onreadystatechange = function() {
            if (req.readyState !== 4) return;
            if (req.status >= 200 && req.status < 300) {
                try {
                    resolve(JSON.parse(req.responseText));
                } catch (parseErr) {
                    reject(parseErr);
                }
                return;
            }
            reject(new Error('Failed to load ' + path + ': HTTP ' + req.status));
        };
        req.send();
    });
}

export function registerPresentationPreset(preset, fallbackName) {
    if (!preset || typeof preset !== 'object') return false;

    var name = (typeof preset.name === 'string' && preset.name.trim())
        ? preset.name.trim()
        : (fallbackName || '');
    if (!name) return false;

    var copy = clonePresentationData(preset);
    copy.name = name;
    PRESENTATION_PRESETS[name] = copy;
    if (PRESENTATION_PRESET_ORDER.indexOf(name) === -1) {
        PRESENTATION_PRESET_ORDER.push(name);
    }
    return true;
}

export function ensurePresentationPresetsLoaded() {
    if (presentationPresetLoadState.loaded) {
        return Promise.resolve(PRESENTATION_PRESETS);
    }
    if (presentationPresetLoadState.promise) {
        return presentationPresetLoadState.promise;
    }

    presentationPresetLoadState.loading = true;
    presentationPresetLoadState.error = null;

    presentationPresetLoadState.promise = requestPresentationJson(PRESENTATION_PRESET_MANIFEST_PATH)
        .then(function(manifest) {
            if (!Array.isArray(manifest)) {
                throw new Error('Invalid presentation preset manifest format.');
            }

            PRESENTATION_PRESETS = {};
            PRESENTATION_PRESET_ORDER = [];

            var jobs = manifest.map(function(entry) {
                if (!entry || typeof entry.file !== 'string') {
                    return Promise.resolve(false);
                }
                return requestPresentationJson(entry.file)
                    .then(function(preset) {
                        var fallbackName = (typeof entry.name === 'string') ? entry.name : '';
                        if (!registerPresentationPreset(preset, fallbackName)) {
                            throw new Error('Invalid preset data in ' + entry.file);
                        }
                        return true;
                    });
            });
            return Promise.all(jobs);
        })
        .then(function() {
            presentationPresetLoadState.loading = false;
            presentationPresetLoadState.loaded = true;
            return PRESENTATION_PRESETS;
        })
        .catch(function(err) {
            presentationPresetLoadState.loading = false;
            presentationPresetLoadState.loaded = false;
            presentationPresetLoadState.error = err;
            console.warn('Presentation presets load failed:', err);
            return PRESENTATION_PRESETS;
        });

    return presentationPresetLoadState.promise;
}


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

var presentationAnnotationState = {
    enabled: true,
    includeInRecording: false,
    notes: {},     // { channel: note } keyed by channel number
    fadeMeta: {},  // { channel: { startTime, duration } } — active fades
    fadeRafId: 0,  // requestAnimationFrame id (0 = not running)
    canvas: null,
    ctx: null,
    resizeBound: false
};

// ── Parameter HUD — live numeric / boolean readouts drawn on the overlay canvas ──
var presentationParamHudState = {
    enabled: true,
    includeInRecording: false,
    items: [],      // array of { path, label }
    anchorX: 0.0,   // 0–1 fractional position (left edge of box)
    anchorY: 1.0,   // 0–1 fractional position (bottom edge of box, so 1=bottom)
    fontSize: 11    // px
};

var presentationUiRefreshAccumulator = 0.0;

function clonePresentationData(value) {
    return JSON.parse(JSON.stringify(value));
}

function presentationClamp(v, lo, hi) {
    return Math.max(lo, Math.min(hi, v));
}

function buildCurrentPresentationAnnotationsConfig() {
    return {
        enabled: !!presentationAnnotationState.enabled,
        includeInRecording: !!presentationAnnotationState.includeInRecording
    };
}

function normalizePresentationAnnotationsConfig(raw) {
    var out = buildCurrentPresentationAnnotationsConfig();
    if (!raw || typeof raw !== 'object') return out;
    if (raw.enabled !== undefined) out.enabled = !!raw.enabled;
    if (raw.includeInRecording !== undefined) {
        out.includeInRecording = !!raw.includeInRecording;
    }
    return out;
}

function normalizePresentationParamHudItems(items, fallback) {
    var src = Array.isArray(items) ? items : (Array.isArray(fallback) ? fallback : []);
    var out = [];
    var seen = {};
    for (var i = 0; i < src.length; i++) {
        var item = src[i];
        if (!item || typeof item !== 'object') continue;
        var path = (typeof item.path === 'string') ? item.path.trim() : '';
        if (!path || seen[path]) continue;
        seen[path] = true;
        out.push({
            path: path,
            label: (typeof item.label === 'string' && item.label.trim()) ? item.label.trim() : path
        });
    }
    return out;
}

function buildCurrentPresentationParamHudConfig() {
    return {
        enabled: !!presentationParamHudState.enabled,
        includeInRecording: !!presentationParamHudState.includeInRecording,
        anchorX: presentationParamHudState.anchorX,
        anchorY: presentationParamHudState.anchorY,
        fontSize: presentationParamHudState.fontSize,
        items: normalizePresentationParamHudItems(presentationParamHudState.items)
    };
}

function normalizePresentationParamHudConfig(raw) {
    var defaults = buildCurrentPresentationParamHudConfig();
    var out = {
        enabled: defaults.enabled,
        includeInRecording: defaults.includeInRecording,
        anchorX: defaults.anchorX,
        anchorY: defaults.anchorY,
        fontSize: defaults.fontSize,
        items: normalizePresentationParamHudItems(defaults.items)
    };
    if (!raw || typeof raw !== 'object') return out;
    if (raw.enabled !== undefined) out.enabled = !!raw.enabled;
    if (raw.includeInRecording !== undefined) {
        out.includeInRecording = !!raw.includeInRecording;
    }
    if (typeof raw.anchorX === 'number' && isFinite(raw.anchorX)) {
        out.anchorX = presentationClamp(raw.anchorX, 0, 1);
    }
    if (typeof raw.anchorY === 'number' && isFinite(raw.anchorY)) {
        out.anchorY = presentationClamp(raw.anchorY, 0, 1);
    }
    if (typeof raw.fontSize === 'number' && isFinite(raw.fontSize)) {
        out.fontSize = Math.max(8, Math.min(48, Math.round(raw.fontSize)));
    }
    out.items = normalizePresentationParamHudItems(raw.items, defaults.items);
    return out;
}

function syncPresentationTimelineUiConfig() {
    if (!presentationState.timeline) return;
    presentationState.timeline.loop = !!presentationState.loop;
    presentationState.timeline.annotations = buildCurrentPresentationAnnotationsConfig();
    presentationState.timeline.paramHud = buildCurrentPresentationParamHudConfig();
}

function applyPresentationTimelineUiConfig(timeline) {
    var annotations = normalizePresentationAnnotationsConfig(timeline && timeline.annotations);
    var paramHud = normalizePresentationParamHudConfig(timeline && timeline.paramHud);

    presentationAnnotationState.enabled = !!annotations.enabled;
    presentationAnnotationState.includeInRecording = !!annotations.includeInRecording;

    presentationParamHudState.enabled = !!paramHud.enabled;
    presentationParamHudState.includeInRecording = !!paramHud.includeInRecording;
    presentationParamHudState.anchorX = paramHud.anchorX;
    presentationParamHudState.anchorY = paramHud.anchorY;
    presentationParamHudState.fontSize = paramHud.fontSize;
    presentationParamHudState.items = normalizePresentationParamHudItems(paramHud.items);

    updatePresentationOverlay();
}

function parseColorHex(hex) {
    if (typeof hex !== 'string') return null;
    var m = /^#([0-9a-f]{6})$/i.exec(hex);
    if (!m) return null;
    var n = parseInt(m[1], 16);
    return {
        r: (n >> 16) & 255,
        g: (n >> 8) & 255,
        b: n & 255
    };
}

function colorWithAlpha(hex, alpha, fallback) {
    var rgb = parseColorHex(hex);
    if (!rgb) return fallback || 'rgba(120, 190, 255, ' + alpha + ')';
    return 'rgba(' + rgb.r + ', ' + rgb.g + ', ' + rgb.b + ', ' + alpha + ')';
}

function ensurePresentationAnnotationCanvas() {
    if (typeof document === 'undefined') return null;
    if (!presentationAnnotationState.canvas) {
        var canvas = document.createElement('canvas');
        canvas.id = 'presentation-annotation-layer';
        canvas.setAttribute('aria-hidden', 'true');
        document.body.appendChild(canvas);
        presentationAnnotationState.canvas = canvas;
        presentationAnnotationState.ctx = canvas.getContext('2d');
    }
    if (!presentationAnnotationState.resizeBound && typeof window !== 'undefined') {
        window.addEventListener('resize', resizePresentationAnnotationCanvas);
        presentationAnnotationState.resizeBound = true;
    }
    resizePresentationAnnotationCanvas();
    return presentationAnnotationState.canvas;
}

function resizePresentationAnnotationCanvas() {
    var canvas = presentationAnnotationState.canvas;
    var ctx = presentationAnnotationState.ctx;
    if (!canvas || !ctx || typeof window === 'undefined') return;

    var width = Math.max(window.innerWidth || 1, 1);
    var height = Math.max(window.innerHeight || 1, 1);
    var dpr = Math.max(window.devicePixelRatio || 1, 1);

    canvas.style.width = width + 'px';
    canvas.style.height = height + 'px';
    canvas.width = Math.max(1, Math.floor(width * dpr));
    canvas.height = Math.max(1, Math.floor(height * dpr));
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
}

function setPresentationAnnotation(note, channel) {
    var ch = (typeof channel === 'number' && channel >= 0) ? channel : 0;
    if (!note || typeof note !== 'object') {
        delete presentationAnnotationState.notes[ch];
        delete presentationAnnotationState.fadeMeta[ch];
        updatePresentationOverlay();
        return false;
    }
    presentationAnnotationState.notes[ch] = clonePresentationData(note);
    var fd = parseFloat(note.fadeIn);
    if (isFinite(fd) && fd > 0) {
        presentationAnnotationState.fadeMeta[ch] = {
            startTime: performance.now(),
            duration: fd * 1000
        };
        startAnnotationFade();
    } else {
        delete presentationAnnotationState.fadeMeta[ch];
    }
    updatePresentationOverlay();
    return true;
}

function clearPresentationAnnotation(channel) {
    if (typeof channel === 'number' && channel >= 0) {
        delete presentationAnnotationState.notes[channel];
        delete presentationAnnotationState.fadeMeta[channel];
    } else {
        presentationAnnotationState.notes = {};
        presentationAnnotationState.fadeMeta = {};
    }
    updatePresentationOverlay();
}

function getChannelFadeAlpha(ch) {
    var meta = presentationAnnotationState.fadeMeta[ch];
    if (!meta) return 1.0;
    var elapsed = performance.now() - meta.startTime;
    if (elapsed >= meta.duration) {
        delete presentationAnnotationState.fadeMeta[ch];
        return 1.0;
    }
    return elapsed / meta.duration;
}

function hasPendingAnnotationFades() {
    return Object.keys(presentationAnnotationState.fadeMeta).length > 0;
}

function annotationFadeTick() {
    presentationAnnotationState.fadeRafId = 0;
    updatePresentationOverlay();
    if (hasPendingAnnotationFades()) {
        presentationAnnotationState.fadeRafId = requestAnimationFrame(annotationFadeTick);
    }
}

function startAnnotationFade() {
    if (presentationAnnotationState.fadeRafId) return; // already ticking
    presentationAnnotationState.fadeRafId = requestAnimationFrame(annotationFadeTick);
}

export function setPresentationAnnotationsEnabled(enabled) {
    presentationAnnotationState.enabled = !!enabled;
    syncPresentationTimelineUiConfig();
    updatePresentationOverlay();
    return presentationAnnotationState.enabled;
}

export function setPresentationAnnotationsIncludedInRecording(enabled) {
    presentationAnnotationState.includeInRecording = !!enabled;
    syncPresentationTimelineUiConfig();
    return presentationAnnotationState.includeInRecording;
}

export function getPresentationAnnotationsState() {
    return {
        enabled: !!presentationAnnotationState.enabled,
        includeInRecording: !!presentationAnnotationState.includeInRecording,
        active: Object.keys(presentationAnnotationState.notes).length > 0
    };
}

// ── Parameter HUD API ────────────────────────────────────────────────────────

export function setPresentationParamHudEnabled(enabled) {
    presentationParamHudState.enabled = !!enabled;
    syncPresentationTimelineUiConfig();
    updatePresentationOverlay();
    return presentationParamHudState.enabled;
}

export function setPresentationParamHudIncludedInRecording(enabled) {
    presentationParamHudState.includeInRecording = !!enabled;
    syncPresentationTimelineUiConfig();
    return presentationParamHudState.includeInRecording;
}

export function getPresentationParamHudState() {
    return {
        enabled: !!presentationParamHudState.enabled,
        includeInRecording: !!presentationParamHudState.includeInRecording,
        anchorX: presentationParamHudState.anchorX,
        anchorY: presentationParamHudState.anchorY,
        fontSize: presentationParamHudState.fontSize,
        items: clonePresentationData(presentationParamHudState.items)
    };
}

export function setParamHudLayout(opts) {
    if (!opts || typeof opts !== 'object') return;
    if (typeof opts.anchorX === 'number' && isFinite(opts.anchorX)) {
        presentationParamHudState.anchorX = Math.max(0, Math.min(1, opts.anchorX));
    }
    if (typeof opts.anchorY === 'number' && isFinite(opts.anchorY)) {
        presentationParamHudState.anchorY = Math.max(0, Math.min(1, opts.anchorY));
    }
    if (typeof opts.fontSize === 'number' && isFinite(opts.fontSize)) {
        presentationParamHudState.fontSize = Math.max(8, Math.min(48, Math.round(opts.fontSize)));
    }
    syncPresentationTimelineUiConfig();
    updatePresentationOverlay();
}

export function isParamInHud(path) {
    for (var i = 0; i < presentationParamHudState.items.length; i++) {
        if (presentationParamHudState.items[i].path === path) return true;
    }
    return false;
}

export function addParamToHud(path, label) {
    if (!path || typeof path !== 'string') return false;
    if (isParamInHud(path)) return false;
    presentationParamHudState.items.push({ path: path, label: label || path });
    syncPresentationTimelineUiConfig();
    updatePresentationOverlay();
    return true;
}

export function removeParamFromHud(path) {
    var before = presentationParamHudState.items.length;
    presentationParamHudState.items = presentationParamHudState.items.filter(function(item) {
        return item.path !== path;
    });
    if (presentationParamHudState.items.length !== before) {
        syncPresentationTimelineUiConfig();
        updatePresentationOverlay();
        return true;
    }
    return false;
}

export function toggleParamInHud(path, label) {
    if (isParamInHud(path)) {
        removeParamFromHud(path);
        return false;
    }
    addParamToHud(path, label);
    return true;
}

export function clearParamHud() {
    presentationParamHudState.items = [];
    syncPresentationTimelineUiConfig();
    updatePresentationOverlay();
}

function formatParamHudValue(val) {
    if (val === undefined || val === null) return '\u2014';
    if (typeof val === 'boolean') return val ? 'true' : 'false';
    if (typeof val === 'number') {
        if (!isFinite(val)) return String(val);
        // Snap floating-point noise near zero to zero
        if (Math.abs(val) < 1e-9) return '0';
        var abs = Math.abs(val);
        // Choose decimal places so we get ~4 significant figures, no sci notation
        var decimals;
        if (abs >= 1000)       decimals = 0;
        else if (abs >= 100)   decimals = 1;
        else if (abs >= 10)    decimals = 2;
        else if (abs >= 1)     decimals = 3;
        else if (abs >= 0.1)   decimals = 4;
        else if (abs >= 0.01)  decimals = 5;
        else                   decimals = 6;
        var s = val.toFixed(decimals);
        // Strip trailing decimal zeros
        if (s.indexOf('.') !== -1) s = s.replace(/\.?0+$/, '');
        return s;
    }
    return String(val);
}

function drawParamHudOnCanvas(ctx, viewWidth, viewHeight) {
    var items = presentationParamHudState.items;
    if (!items.length) return;

    // Collect rows with current live values
    var rows = [];
    for (var i = 0; i < items.length; i++) {
        var item = items[i];
        var val = getPresentationPathValue(item.path);
        rows.push({ label: item.label || item.path, value: formatParamHudValue(val) });
    }
    if (!rows.length) return;

    // Layout constants (scale with user-chosen font size)
    var fs = Math.max(8, Math.min(48, presentationParamHudState.fontSize || 11));
    var fontStr = fs + 'px Consolas, "Courier New", monospace';
    var paddingX = Math.round(fs * 0.9);
    var paddingY = Math.round(fs * 0.7);
    var rowH = Math.round(fs * 1.55);
    var gap = Math.round(fs * 0.7);

    ctx.save();
    ctx.font = fontStr;
    var maxLabelW = 0, maxValueW = 0;
    for (var r = 0; r < rows.length; r++) {
        maxLabelW = Math.max(maxLabelW, ctx.measureText(rows[r].label + ':').width);
        maxValueW = Math.max(maxValueW, ctx.measureText(rows[r].value).width);
    }

    var boxW = paddingX * 2 + maxLabelW + gap + maxValueW;
    var boxH = paddingY * 2 + rows.length * rowH;

    // Anchor: anchorX is box left as fraction of view width,
    // anchorY is box top as fraction of view height (0=top, 1=bottom-aligned).
    // When anchorY===1 the box stays just above the bottom (72px margin).
    var ax = presentationParamHudState.anchorX;
    var ay = presentationParamHudState.anchorY;
    var minMarginX = 8;
    var minMarginY = 8;
    var x, y;
    if (ay >= 1.0) {
        // Legacy bottom-docked behaviour
        x = ax * viewWidth;
        y = viewHeight - boxH - 72;
    } else {
        x = ax * viewWidth;
        y = ay * viewHeight;
    }
    // Clamp so box stays within viewport
    x = Math.max(minMarginX, Math.min(viewWidth  - boxW - minMarginX, x));
    y = Math.max(minMarginY, Math.min(viewHeight - boxH - minMarginY, y));

    // Background
    ctx.shadowBlur = 0;
    drawRoundedRectPath(ctx, x, y, boxW, boxH, Math.round(fs * 0.5));
    ctx.fillStyle = 'rgba(6, 14, 28, 0.82)';
    ctx.fill();
    ctx.lineWidth = 1;
    ctx.strokeStyle = 'rgba(80, 140, 200, 0.45)';
    ctx.stroke();

    // Rows
    ctx.font = fontStr;
    for (var r2 = 0; r2 < rows.length; r2++) {
        var ry = y + paddingY + r2 * rowH + rowH - Math.round(fs * 0.25);
        ctx.fillStyle = '#7bbce8';
        ctx.fillText(rows[r2].label + ':', x + paddingX, ry);
        ctx.fillStyle = '#f0f5ff';
        ctx.fillText(rows[r2].value, x + paddingX + maxLabelW + gap, ry);
    }

    ctx.restore();
}

function wrapCanvasTextLines(ctx, text, maxWidth) {
    var clean = (text || '').toString().replace(/\s+/g, ' ').trim();
    if (!clean) return [];
    var words = clean.split(' ');
    var lines = [];
    var line = words[0];
    for (var i = 1; i < words.length; i++) {
        var test = line + ' ' + words[i];
        if (ctx.measureText(test).width <= maxWidth) {
            line = test;
        } else {
            lines.push(line);
            line = words[i];
        }
    }
    lines.push(line);
    return lines;
}

function drawRoundedRectPath(ctx, x, y, w, h, r) {
    var radius = Math.max(0, Math.min(r, Math.min(w, h) * 0.5));
    ctx.beginPath();
    ctx.moveTo(x + radius, y);
    ctx.lineTo(x + w - radius, y);
    ctx.quadraticCurveTo(x + w, y, x + w, y + radius);
    ctx.lineTo(x + w, y + h - radius);
    ctx.quadraticCurveTo(x + w, y + h, x + w - radius, y + h);
    ctx.lineTo(x + radius, y + h);
    ctx.quadraticCurveTo(x, y + h, x, y + h - radius);
    ctx.lineTo(x, y + radius);
    ctx.quadraticCurveTo(x, y, x + radius, y);
    ctx.closePath();
}

function clampAnnotationAnchorToViewport(x, y, viewWidth, viewHeight) {
    var marginX = Math.min(40, Math.max(12, viewWidth * 0.03));
    var marginY = Math.min(40, Math.max(12, viewHeight * 0.04));
    return {
        x: presentationClamp(x, marginX, Math.max(marginX, viewWidth - marginX)),
        y: presentationClamp(y, marginY, Math.max(marginY, viewHeight - marginY))
    };
}

function rectIntersectsCircle(rx, ry, rw, rh, cx, cy, radius) {
    var nearestX = presentationClamp(cx, rx, rx + rw);
    var nearestY = presentationClamp(cy, ry, ry + rh);
    var dx = nearestX - cx;
    var dy = nearestY - cy;
    return (dx * dx + dy * dy) <= (radius * radius);
}

function getPresentationOverlaySafeMargins(viewWidth) {
    var safe = { left: 14, right: 14 };
    if (typeof document === 'undefined') return safe;

    function includeRect(rect) {
        if (!rect || rect.width < 8 || rect.height < 8) return;
        var gap = 12;
        if (rect.left >= viewWidth * 0.5) {
            safe.right = Math.max(safe.right, Math.max(0, viewWidth - rect.left) + gap);
        } else if (rect.right <= viewWidth * 0.5) {
            safe.left = Math.max(safe.left, Math.max(0, rect.right) + gap);
        }
    }

    var animPanel = document.getElementById('anim-panel');
    if (animPanel && !animPanel.classList.contains('sp-panel--collapsed')) {
        includeRect(animPanel.getBoundingClientRect());
    }

    var guiPanel = document.getElementById('controls-panel');
    if (guiPanel && !guiPanel.classList.contains('sp-panel--collapsed')) {
        includeRect(guiPanel.getBoundingClientRect());
    }

    return safe;
}

function getPresentationAnchorWorldPosition(anchor) {
    if (typeof THREE === 'undefined') return null;

    var a = anchor || {};
    var params = (shader && shader.parameters) ? shader.parameters : null;
    var target = (typeof a.target === 'string') ? a.target.toLowerCase() : '';

    if (target === 'black_hole' || target === 'bh' || target === 'center') {
        return new THREE.Vector3(0.0, 0.0, 0.0);
    }
    if (target === 'disk') {
        var diskR = (params && params.torus && typeof params.torus.r0 === 'number')
            ? params.torus.r0 : 3.4;
        var facing = null;
        if (camera && camera.position) {
            facing = new THREE.Vector3(camera.position.x, camera.position.y, 0.0);
        } else if (typeof observer !== 'undefined' && observer && observer.position) {
            facing = new THREE.Vector3(observer.position.x, observer.position.y, 0.0);
        }
        if (facing && facing.lengthSq() > 1e-6) {
            facing.normalize().multiplyScalar(diskR);
            return facing;
        }
        return new THREE.Vector3(diskR, 0.0, 0.0);
    }
    if (target === 'jet_north') {
        var jetLenN = (params && params.jet && typeof params.jet.length === 'number')
            ? params.jet.length * 0.2 : 6.0;
        return new THREE.Vector3(0.0, 0.0, Math.max(2.0, jetLenN));
    }
    if (target === 'jet_south') {
        var jetLenS = (params && params.jet && typeof params.jet.length === 'number')
            ? params.jet.length * 0.2 : 6.0;
        return new THREE.Vector3(0.0, 0.0, -Math.max(2.0, jetLenS));
    }
    if (target === 'planet') {
        var distance = (params && params.planet && typeof params.planet.distance === 'number')
            ? params.planet.distance : 10.0;
        distance = Math.max(distance, 1.6);

        var tObs = (typeof observer !== 'undefined' &&
            observer && typeof observer.time === 'number')
            ? observer.time : 0.0;

        var orbitalV = 1.0 / Math.sqrt(2.0 * Math.max(distance - 1.0, 0.01));
        var orbitalOmega = -orbitalV *
            Math.sqrt(Math.max(1.0 - 1.0 / distance, 0.0)) / distance;
        var phase = tObs * orbitalOmega;
        return new THREE.Vector3(
            Math.cos(phase) * distance,
            Math.sin(phase) * distance,
            0.0
        );
    }

    var wx = parseFloat(a.x);
    var wy = parseFloat(a.y);
    var wz = parseFloat(a.z);
    return new THREE.Vector3(
        isFinite(wx) ? wx : 0.0,
        isFinite(wy) ? wy : 0.0,
        isFinite(wz) ? wz : 0.0
    );
}

function projectPresentationWorldPoint(worldPoint) {
    if (!worldPoint || typeof THREE === 'undefined' || !camera) return null;
    var projected = worldPoint.clone().project(camera);
    if (!projected ||
        !isFinite(projected.x) ||
        !isFinite(projected.y) ||
        !isFinite(projected.z)) {
        return null;
    }
    var offscreen = projected.z < -1.0 || projected.z > 1.0 ||
        projected.x < -1.0 || projected.x > 1.0 ||
        projected.y < -1.0 || projected.y > 1.0;
    return {
        x: projected.x,
        y: projected.y,
        z: projected.z,
        offscreen: offscreen
    };
}

function getAnnotationAnchorPoint(note, viewWidth, viewHeight) {
    var anchor = note.anchor || {};
    if (anchor.mode === 'world' && typeof THREE !== 'undefined' && camera) {
        var worldPoint = getPresentationAnchorWorldPosition(anchor);
        var projected = projectPresentationWorldPoint(worldPoint);
        if (projected) {
            var ndcX = projected.x;
            var ndcY = projected.y;

            if (projected.offscreen) {
                // If the symbolic target is off-screen, fall back to the black-hole center.
                var centerProjected = projectPresentationWorldPoint(new THREE.Vector3(0.0, 0.0, 0.0));
                if (centerProjected && !centerProjected.offscreen) {
                    ndcX = centerProjected.x;
                    ndcY = centerProjected.y;
                } else {
                    var ndcLen = Math.sqrt(ndcX * ndcX + ndcY * ndcY);
                    if (!isFinite(ndcLen) || ndcLen < 1e-5) {
                        ndcX = 0.0;
                        ndcY = 0.0;
                    } else {
                        ndcX /= ndcLen;
                        ndcY /= ndcLen;
                    }
                    // Keep fallback anchors away from hard viewport edges.
                    ndcX *= 0.72;
                    ndcY *= 0.72;
                }
            } else {
                ndcX = presentationClamp(ndcX, -0.90, 0.90);
                ndcY = presentationClamp(ndcY, -0.90, 0.90);
            }

            var projectedPx = (ndcX * 0.5 + 0.5) * viewWidth;
            var projectedPy = (-ndcY * 0.5 + 0.5) * viewHeight;
            return clampAnnotationAnchorToViewport(projectedPx, projectedPy, viewWidth, viewHeight);
        }
    }

    var sx = (typeof anchor.x === 'number') ? anchor.x : 0.5;
    var sy = (typeof anchor.y === 'number') ? anchor.y : 0.5;
    var screenAnchor = {
        x: presentationClamp(sx, 0.0, 1.0) * viewWidth,
        y: presentationClamp(sy, 0.0, 1.0) * viewHeight
    };
    return clampAnnotationAnchorToViewport(screenAnchor.x, screenAnchor.y, viewWidth, viewHeight);
}

function buildPresentationNoteLayout(ctx, note, viewWidth, viewHeight) {
    var title = (note.title || '').toString();
    var body = (note.text || note.body || '').toString();
    var color = note.color || '#7cc5ff';
    var width = presentationClamp(parseFloat(note.width) || 320, 170, Math.max(180, viewWidth - 40));
    var paddingX = 14;
    var paddingY = 12;
    var lineHeight = 18;
    var bodyLineHeight = 16;

    ctx.font = '600 15px "Segoe UI", "Lucida Grande", sans-serif';
    var bodyMaxWidth = Math.max(120, width - paddingX * 2);
    ctx.font = '13px "Segoe UI", "Lucida Grande", sans-serif';
    var bodyLines = wrapCanvasTextLines(ctx, body, bodyMaxWidth);
    var titleSpace = title ? lineHeight + 2 : 0;
    var bodySpace = bodyLines.length ? (bodyLines.length * bodyLineHeight) : bodyLineHeight;
    var height = paddingY + titleSpace + bodySpace + paddingY;

    var anchor = getAnnotationAnchorPoint(note, viewWidth, viewHeight);

    // ── Direct box position override (set by drag UI) ──
    var hasDirectBox = (typeof note.boxX === 'number' && typeof note.boxY === 'number');
    var x, y;
    if (hasDirectBox) {
        x = presentationClamp(note.boxX * viewWidth, 0, Math.max(0, viewWidth - width));
        y = presentationClamp(note.boxY * viewHeight, 0, Math.max(0, viewHeight - height));
    } else {
        var safeMargins = getPresentationOverlaySafeMargins(viewWidth);
        var placement = (note.placement || 'right').toLowerCase();
        var sideInset = parseFloat(note.sideInset);
        if (!isFinite(sideInset)) sideInset = 26;
        sideInset = presentationClamp(sideInset, 0, Math.max(0, viewWidth * 0.18));
        var sideDockLeft = safeMargins.left + sideInset;
        var sideDockRight = Math.max(sideDockLeft, viewWidth - safeMargins.right - width - sideInset);
        if (placement === 'auto') {
            var usableCenterX = safeMargins.left +
                (viewWidth - safeMargins.left - safeMargins.right) * 0.5;
            placement = (anchor.x >= usableCenterX) ? 'left' : 'right';
        }
        var offset = parseFloat(note.offset);
        if (!isFinite(offset)) offset = 56;
        offset = Math.max(16, offset);
        x = anchor.x;
        y = anchor.y;

        if (placement === 'left') {
            x = sideDockLeft;
            x = Math.min(x, anchor.x - width - offset);
            y -= height * 0.45;
        } else if (placement === 'top') {
            x -= width * 0.5;
            y -= height + offset;
        } else if (placement === 'bottom') {
            x -= width * 0.5;
            y += offset;
        } else {
            x = sideDockRight;
            x = Math.max(x, anchor.x + offset);
            y -= height * 0.45;
        }

        // Keep the center action clear by pushing notes away from the projected black-hole region.
        var bhProj = projectPresentationWorldPoint((typeof THREE !== 'undefined') ? new THREE.Vector3(0.0, 0.0, 0.0) : null);
        if (bhProj && !bhProj.offscreen) {
            var bhX = (bhProj.x * 0.5 + 0.5) * viewWidth;
            var bhY = (-bhProj.y * 0.5 + 0.5) * viewHeight;
            var protectRadius = presentationClamp(Math.min(viewWidth, viewHeight) * 0.24, 140, 340);
            var centerGap = Math.max(18, parseFloat(note.centerGap) || 24);

            if (placement === 'left') {
                var leftMaxX = bhX - protectRadius - width - centerGap;
                if (isFinite(leftMaxX)) x = Math.min(x, leftMaxX);
            } else if (placement === 'right') {
                var rightMinX = bhX + protectRadius + centerGap;
                if (isFinite(rightMinX)) x = Math.max(x, rightMinX);
            } else if (placement === 'top') {
                var topMaxY = bhY - protectRadius - height - centerGap;
                if (isFinite(topMaxY)) y = Math.min(y, topMaxY);
            } else if (placement === 'bottom') {
                var bottomMinY = bhY + protectRadius + centerGap;
                if (isFinite(bottomMinY)) y = Math.max(y, bottomMinY);
            }

            // If still intersecting the protected center, force to the far outer side.
            if (rectIntersectsCircle(x, y, width, height, bhX, bhY, protectRadius)) {
                if (placement === 'left' || placement === 'right') {
                    x = (anchor.x >= bhX) ? sideDockLeft : sideDockRight;
                } else {
                    y = (anchor.y >= bhY)
                        ? Math.max(14, bhY - protectRadius - height - centerGap)
                        : Math.min(viewHeight - height - 14, bhY + protectRadius + centerGap);
                }
            }
        }

        x = presentationClamp(
            x,
            safeMargins.left,
            Math.max(safeMargins.left, viewWidth - width - safeMargins.right)
        );
        y = presentationClamp(y, 14, Math.max(14, viewHeight - height - 14));
    }

    var lineEndX;
    var lineEndY;
    // Compute pointer line end: closest box edge to anchor
    if (hasDirectBox) {
        // With direct positioning, connect to nearest box edge
        var acx = anchor.x, acy = anchor.y;
        var bCenterX = x + width * 0.5, bCenterY = y + height * 0.5;
        var dx = acx - bCenterX, dy = acy - bCenterY;
        if (Math.abs(dx) * height > Math.abs(dy) * width) {
            // Left or right edge
            lineEndX = dx > 0 ? x + width - 8 : x + 8;
            lineEndY = presentationClamp(acy, y + 8, y + height - 8);
        } else {
            // Top or bottom edge
            lineEndX = presentationClamp(acx, x + 8, x + width - 8);
            lineEndY = dy > 0 ? y + height - 8 : y + 8;
        }
    } else {
        var resolvedPlacement = (note.placement || 'right').toLowerCase();
        if (resolvedPlacement === 'auto') {
            resolvedPlacement = (anchor.x >= viewWidth * 0.5) ? 'left' : 'right';
        }
        if (resolvedPlacement === 'left') {
            lineEndX = x + width - 8;
            lineEndY = presentationClamp(anchor.y, y + 8, y + height - 8);
        } else if (resolvedPlacement === 'top') {
            lineEndX = presentationClamp(anchor.x, x + 8, x + width - 8);
            lineEndY = y + height - 8;
        } else if (resolvedPlacement === 'bottom') {
            lineEndX = presentationClamp(anchor.x, x + 8, x + width - 8);
            lineEndY = y + 8;
        } else {
            lineEndX = x + 8;
            lineEndY = presentationClamp(anchor.y, y + 8, y + height - 8);
        }
    }

    return {
        title: title,
        bodyLines: bodyLines,
        color: color,
        anchorX: anchor.x,
        anchorY: anchor.y,
        lineEndX: lineEndX,
        lineEndY: lineEndY,
        x: x,
        y: y,
        width: width,
        height: height,
        paddingX: paddingX,
        paddingY: paddingY,
        lineHeight: lineHeight,
        bodyLineHeight: bodyLineHeight
    };
}

function drawPresentationNote(ctx, layout, alpha) {
    if (!ctx || !layout) return;

    var masterAlpha = (typeof alpha === 'number') ? presentationClamp(alpha, 0, 1) : 1.0;
    if (masterAlpha <= 0) return;

    var accent = layout.color || '#7cc5ff';
    var fill = 'rgba(8, 17, 33, 0.86)';
    var stroke = colorWithAlpha(accent, 0.95, 'rgba(124, 197, 255, 0.95)');
    var line = colorWithAlpha(accent, 0.9, 'rgba(124, 197, 255, 0.9)');
    var glow = colorWithAlpha(accent, 0.3, 'rgba(124, 197, 255, 0.3)');

    ctx.save();
    ctx.globalAlpha = masterAlpha;
    ctx.lineWidth = 2;
    ctx.strokeStyle = line;
    ctx.shadowBlur = 10;
    ctx.shadowColor = glow;
    ctx.beginPath();
    ctx.moveTo(layout.anchorX, layout.anchorY);
    ctx.lineTo(layout.lineEndX, layout.lineEndY);
    ctx.stroke();

    ctx.beginPath();
    ctx.arc(layout.anchorX, layout.anchorY, 4, 0, Math.PI * 2);
    ctx.fillStyle = line;
    ctx.fill();
    ctx.restore();

    ctx.save();
    ctx.globalAlpha = masterAlpha;
    drawRoundedRectPath(ctx, layout.x, layout.y, layout.width, layout.height, 12);
    ctx.fillStyle = fill;
    ctx.fill();
    ctx.lineWidth = 1.4;
    ctx.strokeStyle = stroke;
    ctx.stroke();

    var tx = layout.x + layout.paddingX;
    var ty = layout.y + layout.paddingY;
    if (layout.title) {
        ctx.font = '600 15px "Segoe UI", "Lucida Grande", sans-serif';
        ctx.fillStyle = '#f2f7ff';
        ctx.fillText(layout.title, tx, ty + layout.lineHeight - 3);
        ty += layout.lineHeight + 2;
    }

    ctx.font = '13px "Segoe UI", "Lucida Grande", sans-serif';
    ctx.fillStyle = '#d7e6ff';
    if (!layout.bodyLines.length) layout.bodyLines = [''];
    for (var i = 0; i < layout.bodyLines.length; i++) {
        ctx.fillText(layout.bodyLines[i], tx, ty + layout.bodyLineHeight - 3);
        ty += layout.bodyLineHeight;
    }
    ctx.restore();
}

export function updatePresentationOverlay() {
    var canvas = ensurePresentationAnnotationCanvas();
    var ctx = presentationAnnotationState.ctx;
    if (!canvas || !ctx) return;

    var viewWidth = parseFloat(canvas.style.width) || window.innerWidth || 1;
    var viewHeight = parseFloat(canvas.style.height) || window.innerHeight || 1;
    ctx.clearRect(0, 0, viewWidth, viewHeight);

    if (presentationAnnotationState.enabled) {
        var notes = presentationAnnotationState.notes;
        var channels = Object.keys(notes);
        for (var i = 0; i < channels.length; i++) {
            var ch = channels[i];
            var note = notes[ch];
            if (!note) continue;
            var alpha = getChannelFadeAlpha(ch);
            var layout = buildPresentationNoteLayout(ctx, note, viewWidth, viewHeight);
            drawPresentationNote(ctx, layout, alpha);
        }
    }

    if (presentationParamHudState.enabled && presentationParamHudState.items.length > 0) {
        drawParamHudOnCanvas(ctx, viewWidth, viewHeight);
    }
}

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

function presentationEasing(u, ease) {
    if (ease === 'smooth') {
        return u * u * (3.0 - 2.0 * u);
    }
    if (ease === 'smoother') {
        return u * u * u * (u * (u * 6.0 - 15.0) + 10.0);
    }
    return u;
}

function samplePresentationTrack(track, t) {
    var keys = track.keys;
    if (!keys || !keys.length) return undefined;
    if (keys.length === 1 || t <= keys[0].t) return keys[0].v;
    if (t >= keys[keys.length - 1].t) return keys[keys.length - 1].v;

    for (var i = 0; i < keys.length - 1; i++) {
        var a = keys[i];
        var b = keys[i + 1];
        if (t < a.t || t > b.t) continue;
        var dt = Math.max(b.t - a.t, 1e-8);
        var u = presentationEasing((t - a.t) / dt, b.ease || 'linear');
        if (typeof a.v === 'number' && typeof b.v === 'number') {
            return a.v + (b.v - a.v) * u;
        }
        return (u < 1.0) ? a.v : b.v;
    }

    return keys[keys.length - 1].v;
}

function resolvePresentationPath(path) {
    if (!shader || !path || typeof path !== 'string') return null;

    var clean = path.trim();
    var root = null;
    var parts = [];

    if (clean.indexOf('cameraPan.') === 0) {
        root = cameraPan;
        parts = clean.substring('cameraPan.'.length).split('.');
    } else if (clean.indexOf('camera.') === 0) {
        root = camera;
        parts = clean.substring('camera.'.length).split('.');
    } else if (clean.indexOf('observerState.') === 0) {
        root = observer;
        parts = clean.substring('observerState.'.length).split('.');
    } else if (clean.indexOf('dive.') === 0) {
        root = diveState;
        parts = clean.substring('dive.'.length).split('.');
    } else if (clean.indexOf('hover.') === 0) {
        root = hoverState;
        parts = clean.substring('hover.'.length).split('.');
    } else if (clean.indexOf('params.') === 0) {
        root = shader.parameters;
        parts = clean.substring('params.'.length).split('.');
    } else if (clean.indexOf('shader.parameters.') === 0) {
        root = shader.parameters;
        parts = clean.substring('shader.parameters.'.length).split('.');
    } else {
        root = shader.parameters;
        parts = clean.split('.');
    }

    if (!parts.length) return null;
    return { root: root, parts: parts, originalPath: clean };
}

function presentationPathNeedsCompile(path) {
    if (typeof path !== 'string') return false;
    var clean = path.trim();
    if (!clean) return false;

    if (clean.indexOf('params.') === 0) clean = clean.substring('params.'.length);
    if (clean.indexOf('shader.parameters.') === 0) {
        clean = clean.substring('shader.parameters.'.length);
    }

    switch (clean) {
        case 'kerr_mode':
        case 'accretion_disk':
        case 'accretion_mode':
        case 'disk_self_irradiation':
        case 'jet.enabled':
        case 'jet.mode':
        case 'grmhd.enabled':
        case 'planet.enabled':
        case 'aberration':
        case 'beaming':
        case 'physical_beaming':
        case 'doppler_shift':
        case 'light_travel_time':
        case 'gravitational_time_dilation':
        case 'lorentz_contraction':
        case 'cinematic_tonemap':
        case 'observer.motion':
        case 'quality':
        case 'taa_enabled':
        case 'rk4_integration':
            return true;
        default:
            return false;
    }
}

function refreshPresentationUiBindings() {
    var refreshControllers = (typeof getBlackHoleUiBinding === 'function')
        ? getBlackHoleUiBinding('refreshControllers')
        : refreshAllControllersGlobal;
    if (typeof refreshControllers === 'function') {
        refreshControllers();
    }

    var observerDistanceBinding = (typeof getBlackHoleUiBinding === 'function')
        ? getBlackHoleUiBinding('observerDistance')
        : (typeof distanceController !== 'undefined' ? distanceController : null);
    if (observerDistanceBinding &&
        typeof observerDistanceBinding.updateDisplay === 'function') {
        observerDistanceBinding.updateDisplay();
    }
}

function setPresentationInteractionLock(locked) {
    if (typeof cameraControls !== 'undefined' && cameraControls) {
        cameraControls.enabled = !locked;
    }
}

function setPresentationPathValue(path, value) {
    var cleanPath = (typeof path === 'string') ? path.trim() : '';
    if (cleanPath === 'dive.currentR') {
        var diveRadius = parseFloat(value);
        if (!isFinite(diveRadius) || (!diveState.active && !diveState.reachedSingularity)) {
            return false;
        }
        var diveChanged = Math.abs(diveState.currentR - diveRadius) > 1e-8 ||
            !diveState.timelineDriven || !diveState.paused;
        seekDive(diveRadius, { timelineDriven: true });
        return diveChanged;
    }
    if (cleanPath === 'hover.currentR') {
        var hoverRadius = parseFloat(value);
        if (!isFinite(hoverRadius) || !hoverState.active) {
            return false;
        }
        var hoverChanged = Math.abs(hoverState.currentR - hoverRadius) > 1e-8 ||
            !hoverState.timelineDriven || !hoverState.paused;
        seekHover(hoverRadius, { timelineDriven: true });
        return hoverChanged;
    }

    var resolved = resolvePresentationPath(path);
    if (!resolved) return false;

    var obj = resolved.root;
    var parts = resolved.parts;
    for (var i = 0; i < parts.length - 1; i++) {
        var key = parts[i];
        if (!obj || typeof obj !== 'object' || !(key in obj)) return false;
        obj = obj[key];
    }
    if (!obj || typeof obj !== 'object') return false;

    var leaf = parts[parts.length - 1];
    if (!(leaf in obj)) return false;

    var current = obj[leaf];
    var next = value;
    if (typeof current === 'number') {
        var numeric = parseFloat(next);
        if (!isFinite(numeric)) return false;
        next = numeric;
    } else if (typeof current === 'boolean') {
        next = !!next;
    }

    var changed;
    if (typeof current === 'number' && typeof next === 'number') {
        changed = Math.abs(current - next) > 1e-8;
    } else {
        changed = current !== next;
    }
    if (!changed) return false;

    obj[leaf] = next;

    var isCameraOrObserverPath =
        resolved.originalPath.indexOf('cameraPan.') === 0 ||
        resolved.originalPath.indexOf('camera.') === 0 ||
        resolved.originalPath.indexOf('observerState.') === 0 ||
        resolved.originalPath.indexOf('observer.') === 0 ||
        resolved.originalPath.indexOf('params.observer.') === 0 ||
        resolved.originalPath.indexOf('shader.parameters.observer.') === 0;

    if (resolved.originalPath.indexOf('camera.quaternion.') === 0 &&
        camera && camera.quaternion &&
        typeof camera.quaternion.normalize === 'function') {
        camera.quaternion.normalize();
    }
    if (resolved.originalPath.indexOf('camera.') === 0 &&
        camera && typeof camera.updateMatrixWorld === 'function') {
        camera.updateMatrixWorld(true);
    }

    if (isCameraOrObserverPath && camera && typeof updateCamera === 'function') {
        updateCamera();
    }

    shader.needsUpdate = true;
    return true;
}

export function getPresentationPathValue(path) {
    var resolved = resolvePresentationPath(path);
    if (!resolved) return undefined;

    var obj = resolved.root;
    var parts = resolved.parts;
    for (var i = 0; i < parts.length - 1; i++) {
        var key = parts[i];
        if (!obj || typeof obj !== 'object' || !(key in obj)) return undefined;
        obj = obj[key];
    }
    if (!obj || typeof obj !== 'object') return undefined;

    var leaf = parts[parts.length - 1];
    if (!(leaf in obj)) return undefined;

    var value = obj[leaf];
    if (value && typeof value === 'object') {
        return clonePresentationData(value);
    }
    return value;
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

function choosePresentationMimeType() {
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
    function pad2(v) { return (v < 10 ? '0' : '') + v; }
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

function normalizePresentationRecordingQualityPreset(value) {
    if (typeof value !== 'string') return 'current';
    var clean = value.trim().toLowerCase();
    if (!clean || clean === 'current') return 'current';
    if (typeof QUALITY_PRESETS !== 'undefined' && QUALITY_PRESETS && QUALITY_PRESETS[clean]) {
        return clean;
    }
    return 'current';
}

function normalizePresentationRecordingMode(value) {
    if (typeof value !== 'string') return 'offline';
    var clean = value.trim().toLowerCase();
    if (clean === 'realtime' || clean === 'screen' || clean === 'live') {
        return 'realtime';
    }
    return 'offline';
}

function normalizePresentationRecordingResolutionPreset(value) {
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

function resolvePresentationRecordingResolution(preset) {
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

    // Keep encoder compatibility high: many hardware paths expect even dimensions.
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

function getPresentationRendererRuntimeApi() {
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

function isRendererContextLost() {
    var runtimeApi = getPresentationRendererRuntimeApi();
    return runtimeApi && typeof runtimeApi.isContextLost === 'function' && runtimeApi.isContextLost();
}

function getPresentationWebMMuxerApi() {
    if (typeof window === 'undefined' || !window.WebMMuxer) return null;
    var muxApi = window.WebMMuxer;
    if (typeof muxApi.Muxer !== 'function') return null;
    if (typeof muxApi.ArrayBufferTarget !== 'function') return null;
    return muxApi;
}

function getOfflinePresentationRecordingSupportState() {
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

function isOfflinePresentationRecordingSupported() {
    return getOfflinePresentationRecordingSupportState().supported;
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

function capturePresentationQualitySnapshot() {
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

function restorePresentationQualitySnapshot(snapshot) {
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

function applyPresentationRecordingQualityPreset(presetName) {
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
        // (e.g. recording at 2560×1440 while the window is a different aspect ratio)
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

function stopPresentationCaptureStreamTracks(stream) {
    if (!stream || typeof stream.getTracks !== 'function') return;
    var tracks = stream.getTracks();
    for (var i = 0; i < tracks.length; i++) {
        if (tracks[i] && typeof tracks[i].stop === 'function') {
            tracks[i].stop();
        }
    }
}

function downloadPresentationCaptureDataUrl(dataUrl, mime, filenamePrefix) {
    if (!dataUrl || typeof dataUrl !== 'string') return false;
    var a = document.createElement('a');
    a.href = dataUrl;
    a.download = presentationCaptureFilename(filenamePrefix, mime || 'image/png');
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    return true;
}

function downloadPresentationRecordingBlob(blob, mime, filenamePrefix) {
    if (!blob || blob.size <= 0) return;
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url;
    a.download = presentationCaptureFilename(filenamePrefix, mime);
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function() { URL.revokeObjectURL(url); }, 1500);
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
            console.warn('Offline encoder flush timed out — GPU context may have been lost.');
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
                            // flushes any pending buffered writes) and we're done — no Blob ever
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

        // Always yield via setTimeout — never call encodeNextFrame() directly.
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
