export function presentationEasing(u, ease) {
    if (ease === 'smooth') {
        return u * u * (3.0 - 2.0 * u);
    }
    if (ease === 'smoother') {
        return u * u * u * (u * (u * 6.0 - 15.0) + 10.0);
    }
    return u;
}

export function samplePresentationTrack(track, t) {
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
