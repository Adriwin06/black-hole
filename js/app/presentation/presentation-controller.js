// Role: Presentation timeline and recording controller.
//       Provides presets, keyframe evaluation, scripted events, and
//       MediaRecorder capture hooks for slideshow-ready sequences.

// ─── Presentation Timeline + Capture ─────────────────────────────────────────
// Timeline system for scripted camera/parameter animations suitable for slides
// and demos. Supports keyframed parameter tracks, timed events (dive/hover),
// and optional canvas recording via MediaRecorder.
var PRESENTATION_PRESET_MANIFEST_PATH = 'js/app/presentation/presets/manifest.json';
var PRESENTATION_PRESETS = {};
var PRESENTATION_PRESET_ORDER = [];

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

function registerPresentationPreset(preset, fallbackName) {
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

function ensurePresentationPresetsLoaded() {
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
    fps: 60,
    bitrateMbps: 20.0,
    filenamePrefix: 'black-hole-presentation',
    autoStopOnPresentationEnd: true,
    mimeType: '',
    includeAnnotationsInRecording: false,
    compositeCanvas: null,
    compositeCtx: null,
    compositeRaf: 0
};

var presentationAnnotationState = {
    enabled: true,
    includeInRecording: false,
    note: null,
    canvas: null,
    ctx: null,
    resizeBound: false
};

var presentationUiRefreshAccumulator = 0.0;

function clonePresentationData(value) {
    return JSON.parse(JSON.stringify(value));
}

function presentationClamp(v, lo, hi) {
    return Math.max(lo, Math.min(hi, v));
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

function setPresentationAnnotation(note) {
    if (!note || typeof note !== 'object') {
        presentationAnnotationState.note = null;
        updatePresentationOverlay(true);
        return false;
    }
    presentationAnnotationState.note = clonePresentationData(note);
    updatePresentationOverlay();
    return true;
}

function clearPresentationAnnotation() {
    presentationAnnotationState.note = null;
    updatePresentationOverlay(true);
}

function setPresentationAnnotationsEnabled(enabled) {
    presentationAnnotationState.enabled = !!enabled;
    updatePresentationOverlay();
    return presentationAnnotationState.enabled;
}

function setPresentationAnnotationsIncludedInRecording(enabled) {
    presentationAnnotationState.includeInRecording = !!enabled;
    return presentationAnnotationState.includeInRecording;
}

function getPresentationAnnotationsState() {
    return {
        enabled: !!presentationAnnotationState.enabled,
        includeInRecording: !!presentationAnnotationState.includeInRecording,
        active: !!presentationAnnotationState.note
    };
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

function getAnnotationAnchorPoint(note, viewWidth, viewHeight) {
    var anchor = note.anchor || {};
    if (anchor.mode === 'world' && typeof THREE !== 'undefined' && camera) {
        var worldPoint = getPresentationAnchorWorldPosition(anchor);
        var projected = worldPoint ? worldPoint.clone().project(camera) : null;
        if (projected &&
            isFinite(projected.x) &&
            isFinite(projected.y) &&
            isFinite(projected.z)) {
            var ndcX = projected.x;
            var ndcY = projected.y;
            var offscreen = projected.z < -1.0 || projected.z > 1.0 ||
                ndcX < -1.0 || ndcX > 1.0 || ndcY < -1.0 || ndcY > 1.0;

            if (offscreen) {
                var ndcLen = Math.sqrt(ndcX * ndcX + ndcY * ndcY);
                if (!isFinite(ndcLen) || ndcLen < 1e-5) {
                    ndcX = 0.0;
                    ndcY = 0.0;
                } else {
                    ndcX /= ndcLen;
                    ndcY /= ndcLen;
                }
                ndcX *= 0.94;
                ndcY *= 0.94;
            } else {
                ndcX = presentationClamp(ndcX, -0.96, 0.96);
                ndcY = presentationClamp(ndcY, -0.96, 0.96);
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
    var placement = (note.placement || 'right').toLowerCase();
    if (placement === 'auto') {
        placement = (anchor.x >= viewWidth * 0.55) ? 'left' : 'right';
    }
    var offset = parseFloat(note.offset);
    if (!isFinite(offset)) offset = 38;
    var x = anchor.x;
    var y = anchor.y;

    if (placement === 'left') {
        x -= width + offset;
        y -= height * 0.45;
    } else if (placement === 'top') {
        x -= width * 0.5;
        y -= height + offset;
    } else if (placement === 'bottom') {
        x -= width * 0.5;
        y += offset;
    } else {
        x += offset;
        y -= height * 0.45;
    }

    var margin = 14;
    x = presentationClamp(x, margin, Math.max(margin, viewWidth - width - margin));
    y = presentationClamp(y, margin, Math.max(margin, viewHeight - height - margin));
    var lineEndX = presentationClamp(anchor.x, x + 6, x + width - 6);
    var lineEndY = presentationClamp(anchor.y, y + 6, y + height - 6);

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

function drawPresentationNote(ctx, layout, clearOnly) {
    var canvas = presentationAnnotationState.canvas;
    if (!canvas || !ctx) return;
    var viewWidth = parseFloat(canvas.style.width) || window.innerWidth || 1;
    var viewHeight = parseFloat(canvas.style.height) || window.innerHeight || 1;
    ctx.clearRect(0, 0, viewWidth, viewHeight);
    if (clearOnly) return;

    var accent = layout.color || '#7cc5ff';
    var fill = 'rgba(8, 17, 33, 0.86)';
    var stroke = colorWithAlpha(accent, 0.95, 'rgba(124, 197, 255, 0.95)');
    var line = colorWithAlpha(accent, 0.9, 'rgba(124, 197, 255, 0.9)');
    var glow = colorWithAlpha(accent, 0.3, 'rgba(124, 197, 255, 0.3)');

    ctx.save();
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

function updatePresentationOverlay(forceClear) {
    var canvas = ensurePresentationAnnotationCanvas();
    var ctx = presentationAnnotationState.ctx;
    if (!canvas || !ctx) return;

    if (forceClear || !presentationAnnotationState.enabled || !presentationAnnotationState.note) {
        drawPresentationNote(ctx, null, true);
        return;
    }

    var viewWidth = parseFloat(canvas.style.width) || window.innerWidth || 1;
    var viewHeight = parseFloat(canvas.style.height) || window.innerHeight || 1;
    var layout = buildPresentationNoteLayout(ctx, presentationAnnotationState.note, viewWidth, viewHeight);
    drawPresentationNote(ctx, layout, false);
}

function normalizePresentationTimeline(timeline) {
    if (!timeline) return null;

    var raw = clonePresentationData(timeline);
    var out = {
        name: raw.name || 'Custom',
        loop: !!raw.loop,
        duration: 0.0,
        tracks: [],
        events: []
    };

    var maxTime = 0.0;
    var tracks = Array.isArray(raw.tracks) ? raw.tracks : [];
    for (var i = 0; i < tracks.length; i++) {
        var track = tracks[i];
        if (!track || typeof track.path !== 'string') continue;

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
        var et = parseFloat(ev.t);
        if (!isFinite(et)) continue;
        et = Math.max(0.0, et);
        if (et > maxTime) maxTime = et;
        out.events.push({
            t: et,
            action: ev.action,
            path: ev.path,
            value: ev.value,
            compile: !!ev.compile,
            note: ev.note
        });
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
    if (typeof refreshAllControllersGlobal === 'function') {
        refreshAllControllersGlobal();
    }
    if (typeof distanceController !== 'undefined' &&
        distanceController &&
        typeof distanceController.updateDisplay === 'function') {
        distanceController.updateDisplay();
    }
}

function setPresentationInteractionLock(locked) {
    if (typeof cameraControls !== 'undefined' && cameraControls) {
        cameraControls.enabled = !locked;
    }
}

function setPresentationPathValue(path, value) {
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
        resolved.originalPath.indexOf('observerState.') === 0 ||
        resolved.originalPath.indexOf('observer.') === 0 ||
        resolved.originalPath.indexOf('params.observer.') === 0 ||
        resolved.originalPath.indexOf('shader.parameters.observer.') === 0;

    if (isCameraOrObserverPath && camera && typeof updateCamera === 'function') {
        updateCamera();
    }

    shader.needsUpdate = true;
    return true;
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
            startDive();
            break;
        case 'resetDive':
            resetDive();
            break;
        case 'pauseDive':
            if (diveState.active) {
                diveState.paused = true;
                updateDiveUI();
            }
            break;
        case 'startHover':
            startHover();
            break;
        case 'resetHover':
            resetHover();
            break;
        case 'pauseHover':
            if (hoverState.active) {
                hoverState.paused = true;
                updateHoverUI();
            }
            break;
        case 'updateShader':
            presentationState.compileRequested = true;
            break;
        case 'annotation':
            setPresentationAnnotation(event.note);
            break;
        case 'clearAnnotation':
            clearPresentationAnnotation();
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

function setPresentationTimeline(timeline) {
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

    applyPresentationTracks(0.0);
    shader.needsUpdate = true;
    return true;
}

function listPresentationPresets() {
    return PRESENTATION_PRESET_ORDER.slice();
}

function loadPresentationPreset(name) {
    if (!PRESENTATION_PRESETS[name]) {
        ensurePresentationPresetsLoaded();
        return false;
    }
    var preset = clonePresentationData(PRESENTATION_PRESETS[name]);
    if (!preset.name) preset.name = name;
    return setPresentationTimeline(preset);
}

function seekPresentation(timeSeconds) {
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

function playPresentation(fromStart) {
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
    } else if (presentationState.time <= 1e-6) {
        presentationState.eventCursor = 0;
        processPresentationEvents(-1.0, 0.0);
    }

    presentationState.active = true;
    presentationState.paused = false;
    setPresentationInteractionLock(true);
    shader.needsUpdate = true;
    return true;
}

function pausePresentation() {
    if (!presentationState.timeline) return false;
    presentationState.active = false;
    presentationState.paused = true;
    setPresentationInteractionLock(false);
    return true;
}

function stopPresentation() {
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

function setPresentationLoop(enabled) {
    presentationState.loop = !!enabled;
    return presentationState.loop;
}

function getPresentationState() {
    return {
        loaded: !!presentationState.timeline,
        name: presentationState.timeline ? presentationState.timeline.name : '',
        presets_loaded: !!presentationPresetLoadState.loaded,
        presets_loading: !!presentationPresetLoadState.loading,
        playing: presentationState.active && !presentationState.paused,
        loop: !!presentationState.loop,
        time: presentationState.time,
        duration: presentationState.duration,
        recording: !!presentationCaptureState.active,
        annotations_enabled: !!presentationAnnotationState.enabled,
        annotations_in_recording: !!presentationAnnotationState.includeInRecording
    };
}

function updatePresentation(dt) {
    if (!presentationState.timeline || !presentationState.active || presentationState.paused) return;

    var previousTime = presentationState.time;
    var nextTime = previousTime + dt;
    var duration = Math.max(presentationState.duration, 0.001);

    if (nextTime < duration) {
        processPresentationEvents(previousTime, nextTime);
        presentationState.time = nextTime;
        applyPresentationTracks(nextTime);
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
    var ext = (mimeType && mimeType.indexOf('mp4') !== -1) ? 'mp4' : 'webm';
    return (prefix || 'black-hole-presentation') + '-' + stamp + '.' + ext;
}

function stopPresentationCompositeCapture() {
    if (presentationCaptureState.compositeRaf) {
        cancelAnimationFrame(presentationCaptureState.compositeRaf);
    }
    presentationCaptureState.compositeRaf = 0;
    presentationCaptureState.compositeCanvas = null;
    presentationCaptureState.compositeCtx = null;
}

function stopPresentationRecording() {
    if (!presentationCaptureState.active || !presentationCaptureState.recorder) return false;

    if (presentationCaptureState.recorder.state !== 'inactive') {
        presentationCaptureState.recorder.stop();
    } else if (presentationCaptureState.stream) {
        var tracks = presentationCaptureState.stream.getTracks();
        for (var i = 0; i < tracks.length; i++) tracks[i].stop();
        presentationCaptureState.active = false;
        presentationCaptureState.recorder = null;
        presentationCaptureState.stream = null;
        presentationCaptureState.chunks = [];
        presentationCaptureState.includeAnnotationsInRecording = false;
        stopPresentationCompositeCapture();
    }
    return true;
}

function startPresentationRecording(options) {
    if (!renderer || !renderer.domElement ||
        typeof renderer.domElement.captureStream !== 'function' ||
        typeof MediaRecorder === 'undefined') {
        return false;
    }
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

    var filenamePrefix = options.filenamePrefix || presentationCaptureState.filenamePrefix;
    var autoStop = (options.autoStopOnPresentationEnd === undefined)
        ? presentationCaptureState.autoStopOnPresentationEnd
        : !!options.autoStopOnPresentationEnd;

    var stream = null;
    var compositeTick = null;
    if (includeAnnotationsInRecording) {
        ensurePresentationAnnotationCanvas();
        updatePresentationOverlay();

        var compositeCanvas = document.createElement('canvas');
        var compositeCtx = compositeCanvas.getContext('2d');
        if (!compositeCtx) return false;

        function syncCompositeSize() {
            compositeCanvas.width = Math.max(1, renderer.domElement.width || 1);
            compositeCanvas.height = Math.max(1, renderer.domElement.height || 1);
        }

        syncCompositeSize();
        stream = compositeCanvas.captureStream(fps);
        presentationCaptureState.compositeCanvas = compositeCanvas;
        presentationCaptureState.compositeCtx = compositeCtx;

        compositeTick = function() {
            if (!presentationCaptureState.active || !presentationCaptureState.compositeCtx) return;

            if (compositeCanvas.width !== renderer.domElement.width ||
                compositeCanvas.height !== renderer.domElement.height) {
                syncCompositeSize();
            }

            if (typeof updatePresentationOverlay === 'function') updatePresentationOverlay();

            var ctx = presentationCaptureState.compositeCtx;
            var w = compositeCanvas.width;
            var h = compositeCanvas.height;
            ctx.clearRect(0, 0, w, h);
            ctx.drawImage(renderer.domElement, 0, 0, w, h);
            if (presentationAnnotationState.enabled && presentationAnnotationState.canvas) {
                ctx.drawImage(presentationAnnotationState.canvas, 0, 0, w, h);
            }

            presentationCaptureState.compositeRaf = requestAnimationFrame(compositeTick);
        };
    } else {
        stream = renderer.domElement.captureStream(fps);
        stopPresentationCompositeCapture();
    }

    var mimeType = choosePresentationMimeType();
    var recorderConfig = {};
    if (mimeType) recorderConfig.mimeType = mimeType;
    recorderConfig.videoBitsPerSecond = Math.round(bitrate * 1000000.0);

    var recorder;
    try {
        recorder = new MediaRecorder(stream, recorderConfig);
    } catch (err) {
        // Retry without explicit config if browser rejects codec/bitrate hints
        try {
            recorder = new MediaRecorder(stream);
            mimeType = recorder.mimeType || '';
        } catch (err2) {
            stream.getTracks().forEach(function(track) { track.stop(); });
            stopPresentationCompositeCapture();
            return false;
        }
    }

    presentationCaptureState.active = true;
    presentationCaptureState.recorder = recorder;
    presentationCaptureState.stream = stream;
    presentationCaptureState.chunks = [];
    presentationCaptureState.fps = fps;
    presentationCaptureState.bitrateMbps = bitrate;
    presentationCaptureState.filenamePrefix = filenamePrefix;
    presentationCaptureState.autoStopOnPresentationEnd = autoStop;
    presentationCaptureState.mimeType = mimeType || recorder.mimeType || 'video/webm';
    presentationCaptureState.includeAnnotationsInRecording = includeAnnotationsInRecording;

    if (includeAnnotationsInRecording && compositeTick) {
        presentationCaptureState.compositeRaf = requestAnimationFrame(compositeTick);
    }

    recorder.ondataavailable = function(event) {
        if (event && event.data && event.data.size > 0) {
            presentationCaptureState.chunks.push(event.data);
        }
    };

    recorder.onstop = function() {
        var mime = presentationCaptureState.mimeType || 'video/webm';
        if (presentationCaptureState.chunks.length > 0) {
            var blob = new Blob(presentationCaptureState.chunks, { type: mime });
            var url = URL.createObjectURL(blob);
            var a = document.createElement('a');
            a.href = url;
            a.download = presentationCaptureFilename(filenamePrefix, mime);
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            setTimeout(function() { URL.revokeObjectURL(url); }, 1500);
        }

        if (presentationCaptureState.stream) {
            var tracks = presentationCaptureState.stream.getTracks();
            for (var i = 0; i < tracks.length; i++) tracks[i].stop();
        }

        presentationCaptureState.active = false;
        presentationCaptureState.recorder = null;
        presentationCaptureState.stream = null;
        presentationCaptureState.chunks = [];
        presentationCaptureState.includeAnnotationsInRecording = false;
        stopPresentationCompositeCapture();
    };

    recorder.onerror = function(err) {
        console.warn('Presentation recorder error:', err);
    };

    recorder.start(250);
    shader.needsUpdate = true;
    return true;
}

ensurePresentationPresetsLoaded();

if (typeof window !== 'undefined') {
    window.blackHolePresentation = {
        listPresets: listPresentationPresets,
        ensurePresetsLoaded: ensurePresentationPresetsLoaded,
        loadPreset: loadPresentationPreset,
        setTimeline: setPresentationTimeline,
        play: playPresentation,
        pause: pausePresentation,
        stop: stopPresentation,
        seek: seekPresentation,
        setLoop: setPresentationLoop,
        state: getPresentationState,
        setAnnotationsEnabled: setPresentationAnnotationsEnabled,
        setAnnotationsIncludedInRecording: setPresentationAnnotationsIncludedInRecording,
        annotationState: getPresentationAnnotationsState,
        showAnnotation: setPresentationAnnotation,
        clearAnnotation: clearPresentationAnnotation,
        startRecording: startPresentationRecording,
        stopRecording: stopPresentationRecording
    };
}
