import { THREE } from '../../vendor.js';
import { camera, observer, shader } from '../../core/runtime/runtime-state.js';
var presentationOverlayHooks = {
    onUiConfigChange: null,
    getPathValue: null
};
function clonePresentationData(value) {
    return JSON.parse(JSON.stringify(value));
}
function presentationClamp(v, lo, hi) {
    return Math.max(lo, Math.min(hi, v));
}
export function configurePresentationOverlay(options) {
    presentationOverlayHooks.onUiConfigChange =
        options && typeof options.onUiConfigChange === 'function'
            ? options.onUiConfigChange
            : null;
    presentationOverlayHooks.getPathValue =
        options && typeof options.getPathValue === 'function'
            ? options.getPathValue
            : null;
}
function notifyPresentationUiConfigChanged() {
    if (typeof presentationOverlayHooks.onUiConfigChange === 'function') {
        presentationOverlayHooks.onUiConfigChange();
    }
}

export var presentationAnnotationState = {
    enabled: true,
    includeInRecording: false,
    notes: {},     // { channel: note } keyed by channel number
    fadeMeta: {},  // { channel: { startTime, duration } } â€” active fades
    fadeRafId: 0,  // requestAnimationFrame id (0 = not running)
    canvas: null,
    ctx: null,
    resizeBound: false
};

// â”€â”€ Parameter HUD â€” live numeric / boolean readouts drawn on the overlay canvas â”€â”€
export var presentationParamHudState = {
    enabled: true,
    includeInRecording: false,
    items: [],      // array of { path, label }
    anchorX: 0.0,   // 0â€“1 fractional position (left edge of box)
    anchorY: 1.0,   // 0â€“1 fractional position (bottom edge of box, so 1=bottom)
    fontSize: 11    // px
};

function buildCurrentPresentationAnnotationsConfig() {
    return {
        enabled: !!presentationAnnotationState.enabled,
        includeInRecording: !!presentationAnnotationState.includeInRecording
    };
}

export function normalizePresentationAnnotationsConfig(raw) {
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

export function normalizePresentationParamHudConfig(raw) {
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

export function syncPresentationTimelineUiConfig() {
    notifyPresentationUiConfigChanged();
}


export function applyPresentationTimelineUiConfig(timeline) {
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

export function setPresentationAnnotation(note, channel) {
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

export function clearPresentationAnnotation(channel) {
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

// â”€â”€ Parameter HUD API â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

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
        var val = presentationOverlayHooks.getPathValue ? presentationOverlayHooks.getPathValue(item.path) : undefined;
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

    // â”€â”€ Direct box position override (set by drag UI) â”€â”€
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



