import {
    lookAtOriginQuat,
    sampleTrackDraft
} from './motion-utils.js';

export function setupTimelineMotionModal(options) {
    options = options || {};

    var panel = options.panel || null;
    var motionBtn = options.motionBtn || null;
    var motionModal = options.motionModal || null;
    var motionTypeEl = options.motionTypeEl || null;
    var motionParamsEl = options.motionParamsEl || null;
    var motionApplyBtn = options.motionApplyBtn || null;
    var motionCloseBtn = options.motionCloseBtn || null;
    var esc = options.esc || function(value) { return value; };
    var currentTime = options.currentTime || function() { return 0; };
    var getDraft = options.getDraft || function() { return null; };
    var pushUndo = options.pushUndo || function() {};
    var getTrackByPath = options.getTrackByPath || function() { return null; };
    var upsertKey = options.upsertKey || function() {};
    var normalizeQuatSigns = options.normalizeQuatSigns || function() {};
    var setSelectedTrack = options.setSelectedTrack || function() {};
    var setStatus = options.setStatus || function() {};
    var setDraftDuration = options.setDraftDuration || function() {};
    var sortDraftTracks = options.sortDraftTracks || function() {};
    var applyDraft = options.applyDraft || function() {};
    var rebuildAll = options.rebuildAll || function() {};
    var getPresentationPathValue = options.getPresentationPathValue || function() { return undefined; };

    var MOTION_TYPES = {
        sky_reveal: {
            params: [
                { id: 'start', label: 'Start time (s)', type: 'number', min: 0, step: 0.1, defaultFn: function() { return currentTime().toFixed(2); } },
                { id: 'duration', label: 'Duration (s)', type: 'number', min: 0.1, step: 1, def: 12 },
                { id: 'from', label: 'BH enters from', type: 'select', options: [['bottom', 'Bottom - see sky above'], ['top', 'Top - see sky below'], ['left', 'Left - see sky to the right'], ['right', 'Right - see sky to the left']] },
                { id: 'ease', label: 'Ease', type: 'select', options: [['smoother', 'smoother'], ['smooth', 'smooth'], ['linear', 'linear']] }
            ]
        },
        orbit: {
            params: [
                { id: 'start', label: 'Start time (s)', type: 'number', min: 0, step: 0.1, defaultFn: function() { return currentTime().toFixed(2); } },
                { id: 'duration', label: 'Duration (s)', type: 'number', min: 0.1, step: 1, def: 10 },
                { id: 'orbits', label: 'Number of orbits', type: 'number', min: 0.1, step: 0.5, def: 1 },
                { id: 'dir', label: 'Direction', type: 'select', options: [['ccw', 'Counter-clockwise (CCW)'], ['cw', 'Clockwise (CW)']] }
            ]
        },
        zoom: {
            params: [
                { id: 'start', label: 'Start time (s)', type: 'number', min: 0, step: 0.1, defaultFn: function() { return currentTime().toFixed(2); } },
                { id: 'duration', label: 'Duration (s)', type: 'number', min: 0.1, step: 1, def: 8 },
                { id: 'from', label: 'From distance', type: 'number', min: 1, step: 0.5, defaultFn: function() { var value = getPresentationPathValue('observer.distance'); return (typeof value === 'number' ? value : 15).toFixed(2); } },
                { id: 'to', label: 'To distance', type: 'number', min: 1, step: 0.5, def: 7 },
                { id: 'ease', label: 'Ease', type: 'select', options: [['smooth', 'smooth'], ['smoother', 'smoother'], ['linear', 'linear']] }
            ]
        },
        exposure: {
            params: [
                { id: 'start', label: 'Start time (s)', type: 'number', min: 0, step: 0.1, defaultFn: function() { return currentTime().toFixed(2); } },
                { id: 'duration', label: 'Duration (s)', type: 'number', min: 0.1, step: 1, def: 6 },
                { id: 'from', label: 'From exposure', type: 'number', min: 0, step: 0.05, defaultFn: function() { var value = getPresentationPathValue('look.exposure'); return (typeof value === 'number' ? value : 1.0).toFixed(3); } },
                { id: 'to', label: 'To exposure', type: 'number', min: 0, step: 0.05, def: 1.5 },
                { id: 'ease', label: 'Ease', type: 'select', options: [['smooth', 'smooth'], ['smoother', 'smoother'], ['linear', 'linear']] }
            ]
        },
        inclination: {
            params: [
                { id: 'start', label: 'Start time (s)', type: 'number', min: 0, step: 0.1, defaultFn: function() { return currentTime().toFixed(2); } },
                { id: 'duration', label: 'Duration (s)', type: 'number', min: 0.1, step: 1, def: 8 },
                { id: 'from', label: 'From angle (deg)', type: 'number', min: -90, step: 5, defaultFn: function() { var value = getPresentationPathValue('observer.orbital_inclination'); return (typeof value === 'number' ? value : 0).toFixed(1); } },
                { id: 'to', label: 'To angle (deg)', type: 'number', min: -90, step: 5, def: 30 },
                { id: 'ease', label: 'Ease', type: 'select', options: [['smooth', 'smooth'], ['smoother', 'smoother'], ['linear', 'linear']] }
            ]
        }
    };

    function renderMotionParams() {
        var type = motionTypeEl.value;
        var def = MOTION_TYPES[type];
        if (!def) {
            motionParamsEl.innerHTML = '';
            return;
        }
        var html = '';
        for (var i = 0; i < def.params.length; i++) {
            var p = def.params[i];
            html += '<div class="tl-motion-row"><label>' + esc(p.label) + '</label>';
            if (p.type === 'select') {
                html += '<select data-pid="' + esc(p.id) + '" class="tl-motion-input">';
                for (var j = 0; j < p.options.length; j++) {
                    html += '<option value="' + esc(p.options[j][0]) + '">' + esc(p.options[j][1]) + '</option>';
                }
                html += '</select>';
            } else {
                var defVal = p.defaultFn ? p.defaultFn() : (p.def !== undefined ? p.def : '');
                html += '<input type="number" data-pid="' + esc(p.id) + '" class="tl-motion-input"' +
                    ' step="' + (p.step || 'any') + '"' +
                    ' min="' + (p.min !== undefined ? p.min : '') + '"' +
                    ' value="' + esc(String(defVal)) + '">';
            }
            html += '</div>';
        }
        motionParamsEl.innerHTML = html;
    }

    function getMotionParam(id) {
        var el = motionParamsEl.querySelector('[data-pid="' + CSS.escape(id) + '"]');
        if (!el) return undefined;
        return el.tagName === 'SELECT' ? el.value : parseFloat(el.value);
    }

    function applyMotionFn() {
        var draft = getDraft();
        if (!draft) {
            setStatus('Load or create a timeline first.', 'tl-status--warn');
            return;
        }
        var type = motionTypeEl.value;
        var start = getMotionParam('start');
        if (!isFinite(start)) start = 0;
        var duration = getMotionParam('duration');
        if (!isFinite(duration) || duration < 0.01) duration = 5;
        var ease = getMotionParam('ease') || 'smooth';

        pushUndo();

        if (type === 'orbit') {
            var numOrbits = getMotionParam('orbits');
            if (!isFinite(numOrbits) || numOrbits < 0.01) numOrbits = 1;
            var dirSign = getMotionParam('dir') === 'cw' ? -1 : 1;

            var cpx = 0;
            var cpy = 0;
            var cpz = 11;
            var vx = getPresentationPathValue('camera.position.x');
            var vy = getPresentationPathValue('camera.position.y');
            var vz = getPresentationPathValue('camera.position.z');
            if (typeof vx === 'number') cpx = vx;
            if (typeof vy === 'number') cpy = vy;
            if (typeof vz === 'number') cpz = vz;

            var txd = getTrackByPath('camera.position.x');
            var tyd = getTrackByPath('camera.position.y');
            var tzd = getTrackByPath('camera.position.z');
            if (txd && txd.keys.length) cpx = sampleTrackDraft(txd, start);
            if (tyd && tyd.keys.length) cpy = sampleTrackDraft(tyd, start);
            if (tzd && tzd.keys.length) cpz = sampleTrackDraft(tzd, start);

            var dist = Math.sqrt(cpx * cpx + cpy * cpy + cpz * cpz);
            if (dist < 1e-6) {
                dist = 11;
                cpx = 0;
                cpy = 0;
                cpz = 11;
            }

            var theta = Math.asin(Math.max(-1, Math.min(1, cpy / dist)));
            var phi = Math.atan2(cpz, cpx);
            var cosTheta = Math.cos(theta);
            var totalAngle = dirSign * numOrbits * 2 * Math.PI;
            var numSteps = Math.max(8, Math.round(numOrbits * 32));

            for (var si = 0; si <= numSteps; si++) {
                var frac = si / numSteps;
                var angle = phi + frac * totalAngle;
                var keyTime = start + frac * duration;
                var kx = dist * cosTheta * Math.cos(angle);
                var ky = dist * Math.sin(theta);
                var kz = dist * cosTheta * Math.sin(angle);
                var q = lookAtOriginQuat(kx, ky, kz);
                upsertKey('camera.position.x', keyTime, kx, 'linear');
                upsertKey('camera.position.y', keyTime, ky, 'linear');
                upsertKey('camera.position.z', keyTime, kz, 'linear');
                upsertKey('camera.quaternion.x', keyTime, q.x, 'linear');
                upsertKey('camera.quaternion.y', keyTime, q.y, 'linear');
                upsertKey('camera.quaternion.z', keyTime, q.z, 'linear');
                upsertKey('camera.quaternion.w', keyTime, q.w, 'linear');
            }
            normalizeQuatSigns();
            setSelectedTrack('camera.position.x');
            setStatus('Orbit: ' + numOrbits + ' orbit(s) over ' + duration + 's starting at t=' + start.toFixed(2) + 's.', '');
        } else if (type === 'zoom') {
            var zoomFrom = getMotionParam('from');
            if (!isFinite(zoomFrom)) zoomFrom = 15;
            var zoomTo = getMotionParam('to');
            if (!isFinite(zoomTo)) zoomTo = 7;
            upsertKey('observer.distance', start, zoomFrom, 'linear');
            upsertKey('observer.distance', start + duration, zoomTo, ease);
            setSelectedTrack('observer.distance');
            setStatus('Zoom: distance ' + zoomFrom + ' -> ' + zoomTo + ' over ' + duration + 's.', '');
        } else if (type === 'exposure') {
            var exposureFrom = getMotionParam('from');
            if (!isFinite(exposureFrom)) exposureFrom = 1.0;
            var exposureTo = getMotionParam('to');
            if (!isFinite(exposureTo)) exposureTo = 1.5;
            upsertKey('look.exposure', start, exposureFrom, 'linear');
            upsertKey('look.exposure', start + duration, exposureTo, ease);
            setSelectedTrack('look.exposure');
            setStatus('Exposure: ' + exposureFrom + ' -> ' + exposureTo + ' over ' + duration + 's.', '');
        } else if (type === 'inclination') {
            var inclFrom = getMotionParam('from');
            if (!isFinite(inclFrom)) inclFrom = 0;
            var inclTo = getMotionParam('to');
            if (!isFinite(inclTo)) inclTo = 30;
            upsertKey('observer.orbital_inclination', start, inclFrom, 'linear');
            upsertKey('observer.orbital_inclination', start + duration, inclTo, ease);
            setSelectedTrack('observer.orbital_inclination');
            setStatus('Inclination: ' + inclFrom + ' deg -> ' + inclTo + ' deg over ' + duration + 's.', '');
        } else if (type === 'sky_reveal') {
            var fromDir = getMotionParam('from') || 'bottom';
            var startPanX = 0;
            var startPanY = 0;
            if (fromDir === 'bottom') startPanY = 1.0;
            else if (fromDir === 'top') startPanY = -1.0;
            else if (fromDir === 'left') startPanX = 2.0;
            else if (fromDir === 'right') startPanX = -2.0;

            var endPanX = 0;
            var endPanY = 0;
            var epx = getPresentationPathValue('cameraPan.x');
            var epy = getPresentationPathValue('cameraPan.y');
            if (typeof epx === 'number') endPanX = epx;
            if (typeof epy === 'number') endPanY = epy;

            upsertKey('cameraPan.x', start, endPanX + startPanX, 'linear');
            upsertKey('cameraPan.x', start + duration, endPanX, ease);
            upsertKey('cameraPan.y', start, endPanY + startPanY, 'linear');
            upsertKey('cameraPan.y', start + duration, endPanY, ease);
            setSelectedTrack('cameraPan.x');
            setStatus('Sky reveal: BH enters from ' + fromDir + ' over ' + duration + 's at t=' + start.toFixed(2) + 's.', '');
        }

        setDraftDuration(Math.max(draft.duration, start + duration));
        sortDraftTracks();
        applyDraft();
        rebuildAll();
    }

    function openMotionModal() {
        renderMotionParams();
        motionModal.classList.add('is-open');
        motionBtn.classList.add('is-active');
    }

    function closeMotionModal() {
        motionModal.classList.remove('is-open');
        motionBtn.classList.remove('is-active');
    }

    if (motionBtn) {
        motionBtn.addEventListener('click', function(e) {
            e.stopPropagation();
            if (motionModal.classList.contains('is-open')) closeMotionModal();
            else openMotionModal();
        });
    }
    if (motionCloseBtn) {
        motionCloseBtn.addEventListener('click', function(e) {
            e.stopPropagation();
            closeMotionModal();
        });
    }
    if (motionTypeEl) motionTypeEl.addEventListener('change', renderMotionParams);
    if (motionApplyBtn) {
        motionApplyBtn.addEventListener('click', function(e) {
            e.stopPropagation();
            applyMotionFn();
        });
    }
    if (panel) {
        panel.addEventListener('click', function(e) {
            if (motionModal.classList.contains('is-open') &&
                !motionModal.contains(e.target) &&
                e.target !== motionBtn) {
                closeMotionModal();
            }
        }, true);
    }

    return {
        open: openMotionModal,
        close: closeMotionModal,
        renderParams: renderMotionParams
    };
}
