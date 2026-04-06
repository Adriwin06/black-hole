export function clamp(v, lo, hi) {
    return Math.max(lo, Math.min(hi, v));
}

export function esc(s) {
    return String(s)
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;');
}

export function clonePlain(o) {
    return JSON.parse(JSON.stringify(o));
}

export function parseTime(v, fallback) {
    var n = parseFloat(v);
    return isFinite(n) ? Math.max(0, n) : fallback;
}

export function normalizeEase(e) {
    if (e === 'smooth' || e === 'smoother') return e;
    return 'linear';
}

export function formatValue(v) {
    if (typeof v === 'number') return isFinite(v) ? v.toPrecision(8) : String(v);
    if (typeof v === 'boolean') return v ? 'true' : 'false';
    if (typeof v === 'string') return v;
    if (v === null || typeof v === 'undefined') return '';
    return JSON.stringify(v);
}

export function parseValue(s) {
    if (s === 'true') return true;
    if (s === 'false') return false;
    var n = Number(s);
    if (s !== '' && isFinite(n)) return n;
    return s;
}

export function areValuesEqual(a, b) {
    if (a === b) return true;
    if (typeof a === 'number' && typeof b === 'number') return Math.abs(a - b) < 1e-9;
    return false;
}

export function normalizeHudItems(items) {
    var src = Array.isArray(items) ? items : [];
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

export function normalizeTimelineAnnotationsConfig(raw, fallbackState) {
    var annState = fallbackState || { enabled: true, includeInRecording: false };
    var out = {
        enabled: !!annState.enabled,
        includeInRecording: !!annState.includeInRecording
    };
    if (!raw || typeof raw !== 'object') return out;
    if (raw.enabled !== undefined) out.enabled = !!raw.enabled;
    if (raw.includeInRecording !== undefined) out.includeInRecording = !!raw.includeInRecording;
    return out;
}

export function normalizeTimelineParamHudConfig(raw, fallbackState) {
    var hudState = fallbackState || {
        enabled: true,
        includeInRecording: false,
        anchorX: 0,
        anchorY: 1,
        fontSize: 11,
        items: []
    };
    var out = {
        enabled: !!hudState.enabled,
        includeInRecording: !!hudState.includeInRecording,
        anchorX: isFinite(hudState.anchorX) ? clamp(hudState.anchorX, 0, 1) : 0,
        anchorY: isFinite(hudState.anchorY) ? clamp(hudState.anchorY, 0, 1) : 1,
        fontSize: isFinite(hudState.fontSize) ? clamp(Math.round(hudState.fontSize), 8, 48) : 11,
        items: normalizeHudItems(hudState.items)
    };
    if (!raw || typeof raw !== 'object') return out;
    if (raw.enabled !== undefined) out.enabled = !!raw.enabled;
    if (raw.includeInRecording !== undefined) out.includeInRecording = !!raw.includeInRecording;
    if (typeof raw.anchorX === 'number' && isFinite(raw.anchorX)) out.anchorX = clamp(raw.anchorX, 0, 1);
    if (typeof raw.anchorY === 'number' && isFinite(raw.anchorY)) out.anchorY = clamp(raw.anchorY, 0, 1);
    if (typeof raw.fontSize === 'number' && isFinite(raw.fontSize)) {
        out.fontSize = clamp(Math.round(raw.fontSize), 8, 48);
    }
    out.items = normalizeHudItems(raw.items || out.items);
    return out;
}

export function normalizeTimelineDraft(raw, options) {
    var fallbackAnnotationsState = options && options.annotationsState;
    var fallbackParamHudState = options && options.paramHudState;
    var src = (raw && typeof raw === 'object') ? clonePlain(raw) : {};
    var out = {
        name: src.name || 'Untitled',
        duration: Math.max(0.5, parseTime(src.duration, 12)),
        loop: !!src.loop,
        tracks: [],
        events: [],
        annotationTracks: [],
        annotations: normalizeTimelineAnnotationsConfig(src.annotations, fallbackAnnotationsState),
        paramHud: normalizeTimelineParamHudConfig(src.paramHud, fallbackParamHudState)
    };

    if (Array.isArray(src.annotationTracks)) {
        for (var ati = 0; ati < src.annotationTracks.length; ati++) {
            var at = src.annotationTracks[ati];
            if (at && typeof at === 'object') {
                out.annotationTracks.push({ label: at.label || ('Annotation ' + (ati + 1)) });
            }
        }
    }
    if (!out.annotationTracks.length) out.annotationTracks.push({ label: 'Annotation 1' });

    var tracks = Array.isArray(src.tracks) ? src.tracks : [];
    for (var i = 0; i < tracks.length; i++) {
        var tr = tracks[i];
        if (!tr || !tr.path) continue;
        var keys = [];
        var srcK = Array.isArray(tr.keys) ? tr.keys : [];
        for (var k = 0; k < srcK.length; k++) {
            var t = parseTime(srcK[k] && srcK[k].t, NaN);
            if (!isFinite(t)) continue;
            keys.push({ t: t, v: srcK[k].v, ease: normalizeEase(srcK[k].ease) });
        }
        if (!keys.length) continue;
        keys.sort(function(a, b) { return a.t - b.t; });
        out.tracks.push({ path: tr.path.trim(), compile: !!tr.compile, keys: keys });
    }
    out.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });

    var events = Array.isArray(src.events) ? src.events : [];
    for (var j = 0; j < events.length; j++) {
        var ev = events[j];
        if (!ev || !ev.action) continue;
        var et = parseTime(ev.t, NaN);
        if (!isFinite(et)) continue;
        out.events.push(clonePlain(ev));
    }
    out.events.sort(function(a, b) {
        return parseTime(a.t, 0) - parseTime(b.t, 0);
    });

    return out;
}
