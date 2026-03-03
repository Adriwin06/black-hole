// Role: Interactive presentation timeline editor embedded in the
//       presentation panel. Supports new/apply, keyframes, events, and JSON.

var PRESENTATION_EDITOR_COMMON_PATHS = [
    'observer.distance',
    'observer.orbital_inclination',
    'observer.azimuth',
    'observer.motion',
    'look.exposure',
    'look.disk_gain',
    'look.glow',
    'black_hole.spin',
    'black_hole.spin_enabled',
    'accretion_mode',
    'kerr_mode',
    'jet.enabled',
    'jet.mode',
    'grmhd.enabled',
    'turbulence_loop_enabled',
    'turbulence_loop_seconds',
    'cameraPan.x',
    'cameraPan.y',
    'camera.position.x',
    'camera.position.y',
    'camera.position.z',
    'camera.quaternion.x',
    'camera.quaternion.y',
    'camera.quaternion.z',
    'camera.quaternion.w'
];

function createPresentationEditorHtml() {
    return '' +
        '<div id="presentation-editor-root" class="presentation-editor-root is-hidden">' +
            '<button id="presentation-editor-toggle" class="presentation-editor-toggle" type="button" aria-expanded="false" aria-controls="presentation-editor-body">' +
                '<span class="presentation-editor-toggle-arrow">&#9654;</span> TIMELINE EDITOR' +
            '</button>' +
            '<div id="presentation-editor-body" class="presentation-editor-body">' +

            '<div class="presentation-row">' +
                '<label for="presentation-editor-name">Name</label>' +
                '<input id="presentation-editor-name" class="presentation-input-wide" type="text" placeholder="New Presentation">' +
            '</div>' +

            '<div class="presentation-row">' +
                '<label for="presentation-editor-duration">Duration</label>' +
                '<input id="presentation-editor-duration" type="number" min="0.5" max="1200" step="0.1" value="12">' +
                '<label class="presentation-inline-toggle">' +
                    '<input id="presentation-editor-loop" type="checkbox"> Loop' +
                '</label>' +
            '</div>' +

            '<div class="presentation-btn-row presentation-btn-row-tight">' +
                '<button id="presentation-editor-new-btn" class="presentation-btn" type="button">NEW</button>' +
                '<button id="presentation-editor-from-loaded-btn" class="presentation-btn" type="button">FROM LOADED</button>' +
                '<button id="presentation-editor-apply-btn" class="presentation-btn presentation-btn-primary" type="button">APPLY</button>' +
            '</div>' +

            '<details class="presentation-editor-section">' +
                '<summary class="presentation-editor-section-title">Simulation loop</summary>' +
                '<div class="presentation-editor-section-body">' +
                    '<div class="presentation-row presentation-toggles">' +
                        '<label><input id="presentation-editor-turbulence-loop-enabled" type="checkbox"> Loop turbulence</label>' +
                    '</div>' +
                    '<div class="presentation-row">' +
                        '<label for="presentation-editor-turbulence-loop-seconds">Period</label>' +
                        '<input id="presentation-editor-turbulence-loop-seconds" type="number" min="0.25" max="240" step="0.25" value="20">' +
                        '<label class="presentation-inline-label">s</label>' +
                    '</div>' +
                '</div>' +
            '</details>' +

            '<details class="presentation-editor-section" open>' +
                '<summary class="presentation-editor-section-title">Tracks / keyframes</summary>' +
                '<div class="presentation-editor-section-body">' +
                    '<div class="presentation-row">' +
                        '<label for="presentation-editor-track-path">Path</label>' +
                        '<input id="presentation-editor-track-path" class="presentation-input-wide" list="presentation-editor-path-list" type="text" placeholder="observer.distance">' +
                    '</div>' +

                    '<div class="presentation-row">' +
                        '<label for="presentation-editor-key-time">Key</label>' +
                        '<input id="presentation-editor-key-time" type="number" min="0" step="0.01" value="0">' +
                        '<select id="presentation-editor-key-ease" class="presentation-select presentation-editor-compact-select">' +
                            '<option value="linear">linear</option>' +
                            '<option value="smooth">smooth</option>' +
                            '<option value="smoother">smoother</option>' +
                        '</select>' +
                    '</div>' +

                    '<div class="presentation-row">' +
                        '<label for="presentation-editor-key-value">Value</label>' +
                        '<input id="presentation-editor-key-value" class="presentation-input-wide" type="text" placeholder="11 or true or thin">' +
                    '</div>' +

                    '<div class="presentation-btn-row presentation-btn-row-tight">' +
                        '<button id="presentation-editor-key-time-now-btn" class="presentation-btn" type="button">USE TIME</button>' +
                        '<button id="presentation-editor-key-capture-btn" class="presentation-btn" type="button">CAPTURE VALUE</button>' +
                        '<button id="presentation-editor-key-add-btn" class="presentation-btn presentation-btn-primary" type="button">ADD/UPDATE</button>' +
                        '<button id="presentation-editor-key-remove-btn" class="presentation-btn" type="button">REMOVE</button>' +
                    '</div>' +
                    '<div class="presentation-btn-row presentation-btn-row-tight">' +
                        '<button id="presentation-editor-auto-keyframe-btn" class="presentation-btn presentation-btn-primary" type="button">AUTO KEYFRAME (FROM CONTROLS)</button>' +
                    '</div>' +
                    '<div class="presentation-editor-list-toolbar">' +
                        '<button id="presentation-editor-track-list-expand-btn" class="presentation-mini-btn" type="button" aria-pressed="false">EXPAND LIST</button>' +
                        '<span class="presentation-editor-list-hint">Drag bottom edge to resize</span>' +
                    '</div>' +

                    '<div id="presentation-editor-track-list" class="presentation-editor-list"></div>' +
                '</div>' +
            '</details>' +

            '<details class="presentation-editor-section">' +
                '<summary class="presentation-editor-section-title">Events</summary>' +
                '<div class="presentation-editor-section-body">' +
                    '<div class="presentation-row">' +
                        '<label for="presentation-editor-event-time">Event</label>' +
                        '<input id="presentation-editor-event-time" type="number" min="0" step="0.01" value="0">' +
                        '<select id="presentation-editor-event-action" class="presentation-select presentation-editor-action-select">' +
                            '<option value="set">set</option>' +
                            '<option value="updateShader">updateShader</option>' +
                            '<option value="startDive">startDive</option>' +
                            '<option value="pauseDive">pauseDive</option>' +
                            '<option value="resetDive">resetDive</option>' +
                            '<option value="startHover">startHover</option>' +
                            '<option value="pauseHover">pauseHover</option>' +
                            '<option value="resetHover">resetHover</option>' +
                            '<option value="annotation">annotation</option>' +
                            '<option value="clearAnnotation">clearAnnotation</option>' +
                        '</select>' +
                    '</div>' +

                    '<div id="presentation-editor-event-set-fields">' +
                        '<div class="presentation-row">' +
                            '<label for="presentation-editor-event-path">Path</label>' +
                            '<input id="presentation-editor-event-path" class="presentation-input-wide" list="presentation-editor-path-list" type="text" placeholder="accretion_mode">' +
                        '</div>' +
                        '<div class="presentation-row">' +
                            '<label for="presentation-editor-event-value">Value</label>' +
                            '<input id="presentation-editor-event-value" class="presentation-input-wide" type="text" placeholder="thin">' +
                        '</div>' +
                        '<div class="presentation-row presentation-toggles">' +
                            '<label><input id="presentation-editor-event-compile" type="checkbox"> Force compile</label>' +
                        '</div>' +
                    '</div>' +

                    '<div id="presentation-editor-event-annotation-fields" class="presentation-editor-note-fields">' +
                        '<div class="presentation-row">' +
                            '<label for="presentation-editor-note-title">Title</label>' +
                            '<input id="presentation-editor-note-title" class="presentation-input-wide" type="text" placeholder="Feature">' +
                        '</div>' +
                        '<div class="presentation-row">' +
                            '<label for="presentation-editor-note-text">Text</label>' +
                            '<textarea id="presentation-editor-note-text" class="presentation-editor-note-text" rows="2" placeholder="Explain what changed"></textarea>' +
                        '</div>' +
                        '<div class="presentation-row">' +
                            '<label for="presentation-editor-note-target">Anchor</label>' +
                            '<select id="presentation-editor-note-target" class="presentation-select">' +
                                '<option value="black_hole">black_hole</option>' +
                                '<option value="disk">disk</option>' +
                                '<option value="jet_north">jet_north</option>' +
                                '<option value="jet_south">jet_south</option>' +
                                '<option value="planet">planet</option>' +
                            '</select>' +
                        '</div>' +
                        '<div class="presentation-row">' +
                            '<label for="presentation-editor-note-placement">Placement</label>' +
                            '<select id="presentation-editor-note-placement" class="presentation-select">' +
                                '<option value="auto">auto</option>' +
                                '<option value="left">left</option>' +
                                '<option value="right">right</option>' +
                                '<option value="top">top</option>' +
                                '<option value="bottom">bottom</option>' +
                            '</select>' +
                        '</div>' +
                    '</div>' +

                    '<div class="presentation-btn-row presentation-btn-row-tight">' +
                        '<button id="presentation-editor-event-time-now-btn" class="presentation-btn" type="button">USE TIME</button>' +
                        '<button id="presentation-editor-event-add-btn" class="presentation-btn presentation-btn-primary" type="button">ADD EVENT</button>' +
                    '</div>' +

                    '<div id="presentation-editor-event-list" class="presentation-editor-list"></div>' +
                '</div>' +
            '</details>' +

            '<details class="presentation-editor-section">' +
                '<summary class="presentation-editor-section-title">JSON</summary>' +
                '<div class="presentation-editor-section-body">' +
                    '<textarea id="presentation-editor-json" class="presentation-editor-json" rows="5" spellcheck="false"></textarea>' +
                    '<div class="presentation-btn-row presentation-btn-row-tight">' +
                        '<button id="presentation-editor-json-export-btn" class="presentation-btn" type="button">EXPORT</button>' +
                        '<button id="presentation-editor-json-import-file-btn" class="presentation-btn" type="button">IMPORT FILE</button>' +
                        '<button id="presentation-editor-json-import-btn" class="presentation-btn" type="button">IMPORT TEXT</button>' +
                        '<button id="presentation-editor-json-download-btn" class="presentation-btn" type="button">DOWNLOAD</button>' +
                    '</div>' +
                    '<input id="presentation-editor-json-file-input" type="file" accept=".json,application/json" hidden>' +
                '</div>' +
            '</details>' +

            '<div id="presentation-editor-status" class="presentation-editor-status">No draft loaded.</div>' +
            '<datalist id="presentation-editor-path-list"></datalist>' +
            '</div>' +
        '</div>';
}

function bindPresentationTimelineEditor(section, hooks) {
    if (!section) return null;
    hooks = hooks || {};

    var root = section.querySelector('#presentation-editor-root');
    if (!root) return null;
    var toggleBtn = section.querySelector('#presentation-editor-toggle');
    var body = section.querySelector('#presentation-editor-body');

    var nameInput = section.querySelector('#presentation-editor-name');
    var durationInput = section.querySelector('#presentation-editor-duration');
    var loopInput = section.querySelector('#presentation-editor-loop');
    var newBtn = section.querySelector('#presentation-editor-new-btn');
    var fromLoadedBtn = section.querySelector('#presentation-editor-from-loaded-btn');
    var applyBtn = section.querySelector('#presentation-editor-apply-btn');
    var turbulenceLoopEnabledInput = section.querySelector('#presentation-editor-turbulence-loop-enabled');
    var turbulenceLoopSecondsInput = section.querySelector('#presentation-editor-turbulence-loop-seconds');
    var trackPathInput = section.querySelector('#presentation-editor-track-path');
    var keyTimeInput = section.querySelector('#presentation-editor-key-time');
    var keyEaseSelect = section.querySelector('#presentation-editor-key-ease');
    var keyValueInput = section.querySelector('#presentation-editor-key-value');
    var keyTimeNowBtn = section.querySelector('#presentation-editor-key-time-now-btn');
    var keyCaptureBtn = section.querySelector('#presentation-editor-key-capture-btn');
    var keyAddBtn = section.querySelector('#presentation-editor-key-add-btn');
    var keyRemoveBtn = section.querySelector('#presentation-editor-key-remove-btn');
    var autoKeyframeBtn = section.querySelector('#presentation-editor-auto-keyframe-btn');
    var trackListExpandBtn = section.querySelector('#presentation-editor-track-list-expand-btn');
    var trackList = section.querySelector('#presentation-editor-track-list');
    var eventTimeInput = section.querySelector('#presentation-editor-event-time');
    var eventActionSelect = section.querySelector('#presentation-editor-event-action');
    var eventSetFields = section.querySelector('#presentation-editor-event-set-fields');
    var eventPathInput = section.querySelector('#presentation-editor-event-path');
    var eventValueInput = section.querySelector('#presentation-editor-event-value');
    var eventCompileInput = section.querySelector('#presentation-editor-event-compile');
    var eventAnnotationFields = section.querySelector('#presentation-editor-event-annotation-fields');
    var noteTitleInput = section.querySelector('#presentation-editor-note-title');
    var noteTextInput = section.querySelector('#presentation-editor-note-text');
    var noteTargetSelect = section.querySelector('#presentation-editor-note-target');
    var notePlacementSelect = section.querySelector('#presentation-editor-note-placement');
    var eventTimeNowBtn = section.querySelector('#presentation-editor-event-time-now-btn');
    var eventAddBtn = section.querySelector('#presentation-editor-event-add-btn');
    var eventList = section.querySelector('#presentation-editor-event-list');
    var jsonTextarea = section.querySelector('#presentation-editor-json');
    var jsonExportBtn = section.querySelector('#presentation-editor-json-export-btn');
    var jsonImportFileBtn = section.querySelector('#presentation-editor-json-import-file-btn');
    var jsonImportBtn = section.querySelector('#presentation-editor-json-import-btn');
    var jsonDownloadBtn = section.querySelector('#presentation-editor-json-download-btn');
    var jsonFileInput = section.querySelector('#presentation-editor-json-file-input');
    var statusEl = section.querySelector('#presentation-editor-status');
    var pathList = section.querySelector('#presentation-editor-path-list');

    var draft = null;
    var draftDirty = false;
    var trackListExpanded = false;
    var autoKeyframeState = {
        lastSnapshot: null
    };

    function clonePlain(v) { return JSON.parse(JSON.stringify(v)); }
    function clamp(v, lo, hi) { return Math.max(lo, Math.min(hi, v)); }
    function parseTime(v, fallback) {
        var t = parseFloat(v);
        if (!isFinite(t)) t = isFinite(fallback) ? fallback : 0;
        return Math.max(0, t);
    }
    function parseValue(text) {
        var raw = (text || '').toString().trim();
        if (!raw) return '';
        if (raw === 'true') return true;
        if (raw === 'false') return false;
        if (/^[-+]?(\d+(\.\d+)?|\.\d+)(e[-+]?\d+)?$/i.test(raw)) return parseFloat(raw);
        if ((raw[0] === '{' && raw[raw.length - 1] === '}') ||
            (raw[0] === '[' && raw[raw.length - 1] === ']') ||
            (raw[0] === '"' && raw[raw.length - 1] === '"')) {
            try { return JSON.parse(raw); } catch (err) {}
        }
        return raw;
    }
    function areValuesEqual(a, b) {
        if (typeof a === 'number' && typeof b === 'number') {
            return Math.abs(a - b) <= 1e-8;
        }
        return a === b;
    }
    function formatValue(v) {
        if (typeof v === 'string') return v;
        try { return JSON.stringify(v); } catch (err) { return String(v); }
    }
    function normalizeEase(v) {
        if (v === 'smooth' || v === 'smoother') return v;
        return 'linear';
    }
    function esc(text) {
        return String(text)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#39;');
    }
    function setEditorStatus(text, mode) {
        if (!statusEl) return;
        statusEl.textContent = text;
        statusEl.classList.remove('is-warning', 'is-error', 'is-dirty', 'is-clean');
        if (mode) statusEl.classList.add(mode);
    }
    function setMainStatus(text, mode) {
        if (hooks && typeof hooks.setStatus === 'function') hooks.setStatus(text, mode);
    }
    function setTrackListExpanded(expanded) {
        trackListExpanded = !!expanded;
        if (trackList) {
            trackList.classList.toggle('is-expanded', trackListExpanded);
        }
        if (trackListExpandBtn) {
            trackListExpandBtn.textContent = trackListExpanded ? 'COLLAPSE LIST' : 'EXPAND LIST';
            trackListExpandBtn.setAttribute('aria-pressed', trackListExpanded ? 'true' : 'false');
            trackListExpandBtn.title = trackListExpanded
                ? 'Collapse keyframe list height'
                : 'Expand keyframe list height';
        }
    }
    function importJsonTextIntoDraft(rawText, sourceLabel) {
        var text = (rawText || '').trim();
        if (!text) {
            setEditorStatus('JSON text is empty.', 'is-warning');
            return false;
        }
        try {
            setDraft(JSON.parse(text), true);
            var src = sourceLabel ? (' from ' + sourceLabel) : '';
            setEditorStatus('JSON imported' + src + '. Press APPLY to use it.', 'is-dirty');
            return true;
        } catch (err) {
            setEditorStatus('Invalid JSON: ' + err.message, 'is-error');
            return false;
        }
    }
    function setOpen(open) {
        if (!root || !toggleBtn || !body) return;
        var isOpen = !!open;
        root.classList.toggle('is-open', isOpen);
        toggleBtn.setAttribute('aria-expanded', isOpen ? 'true' : 'false');
    }
    function setVisible(visible) {
        if (!root) return;
        root.classList.toggle('is-hidden', !visible);
    }
    function startNewPresetDraft() {
        if (draftDirty) {
            var proceed = window.confirm('Discard draft changes and start a new preset?');
            if (!proceed) return false;
        }
        setDraft(createEmptyTimeline(), false);
        autoKeyframeState.lastSnapshot = null;
        setVisible(true);
        setOpen(true);
        setEditorStatus('New preset draft ready. Edit then press APPLY.', 'is-clean');
        return true;
    }
    function getRuntimeTimeline() {
        if (typeof getPresentationTimeline === 'function') return getPresentationTimeline();
        if (window.blackHolePresentation && typeof window.blackHolePresentation.getTimeline === 'function') {
            return window.blackHolePresentation.getTimeline();
        }
        return null;
    }
    function getPathValue(path) {
        if (typeof getPresentationPathValue === 'function') return getPresentationPathValue(path);
        if (window.blackHolePresentation && typeof window.blackHolePresentation.getPathValue === 'function') {
            return window.blackHolePresentation.getPathValue(path);
        }
        return undefined;
    }
    function setPathValue(path, value) {
        if (typeof setPresentationPathValue === 'function') return !!setPresentationPathValue(path, value);
        if (typeof shader === 'undefined' || !shader || !shader.parameters) return false;

        if (path === 'turbulence_loop_enabled') {
            var nextEnabled = !!value;
            var enabledChanged = shader.parameters.turbulence_loop_enabled !== nextEnabled;
            shader.parameters.turbulence_loop_enabled = nextEnabled;
            if (enabledChanged) shader.needsUpdate = true;
            return enabledChanged;
        }
        if (path === 'turbulence_loop_seconds') {
            var nextSeconds = clamp(parseFloat(value) || 20.0, 0.25, 240.0);
            var prevSeconds = parseFloat(shader.parameters.turbulence_loop_seconds) || 20.0;
            var secondsChanged = Math.abs(prevSeconds - nextSeconds) > 1e-8;
            shader.parameters.turbulence_loop_seconds = nextSeconds;
            if (secondsChanged) shader.needsUpdate = true;
            return secondsChanged;
        }
        return false;
    }
    function syncTurbulenceLoopInputs() {
        if (!turbulenceLoopEnabledInput || !turbulenceLoopSecondsInput) return;
        var enabled = !!getPathValue('turbulence_loop_enabled');
        var period = parseFloat(getPathValue('turbulence_loop_seconds'));
        if (!isFinite(period)) period = 20.0;
        period = clamp(period, 0.25, 240.0);

        turbulenceLoopEnabledInput.checked = enabled;
        turbulenceLoopSecondsInput.value = period.toFixed(2);
        turbulenceLoopSecondsInput.disabled = !enabled;
    }
    function normalizeTimeline(raw) {
        var src = (raw && typeof raw === 'object') ? clonePlain(raw) : {};
        var out = {
            name: (typeof src.name === 'string' && src.name.trim()) ? src.name.trim() : 'Untitled',
            duration: Math.max(0.5, parseTime(src.duration, 12)),
            loop: !!src.loop,
            tracks: [],
            events: []
        };
        var maxT = 0;
        var tracks = Array.isArray(src.tracks) ? src.tracks : [];
        for (var i = 0; i < tracks.length; i++) {
            var tr = tracks[i];
            if (!tr || typeof tr.path !== 'string' || !tr.path.trim()) continue;
            var keys = [];
            var srcKeys = Array.isArray(tr.keys) ? tr.keys : [];
            for (var k = 0; k < srcKeys.length; k++) {
                var key = srcKeys[k];
                if (!key) continue;
                var t = parseTime(key.t, NaN);
                if (!isFinite(t)) continue;
                maxT = Math.max(maxT, t);
                keys.push({ t: t, v: key.v, ease: normalizeEase(key.ease) });
            }
            if (!keys.length) continue;
            keys.sort(function(a, b) { return a.t - b.t; });
            out.tracks.push({ path: tr.path.trim(), compile: !!tr.compile, keys: keys });
        }
        out.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });

        var events = Array.isArray(src.events) ? src.events : [];
        for (var j = 0; j < events.length; j++) {
            var ev = events[j];
            if (!ev || typeof ev.action !== 'string' || !ev.action.trim()) continue;
            var et = parseTime(ev.t, NaN);
            if (!isFinite(et)) continue;
            maxT = Math.max(maxT, et);
            var evOut = clonePlain(ev);
            evOut.t = et;
            evOut.action = ev.action.trim();
            out.events.push(evOut);
        }
        out.events.sort(function(a, b) { return a.t - b.t; });
        out.duration = Math.max(out.duration, maxT, 0.5);
        return out;
    }
    function createEmptyTimeline() {
        return normalizeTimeline({ name: 'New Presentation', duration: 12, loop: false, tracks: [], events: [] });
    }
    function currentTime() {
        if (typeof getPresentationState === 'function') {
            var st = getPresentationState();
            if (st && isFinite(st.time)) return Math.max(0, st.time);
        }
        return 0;
    }
    function goToTime(t) {
        if (typeof seekPresentation === 'function') seekPresentation(t);
        if (hooks && typeof hooks.syncFromState === 'function') hooks.syncFromState();
    }
    function markDirty(text) {
        draftDirty = true;
        setEditorStatus(text || 'Draft changed. Press APPLY to use it.', 'is-dirty');
    }
    function setDraft(timeline, dirty) {
        draft = normalizeTimeline(timeline);
        draftDirty = !!dirty;
        autoKeyframeState.lastSnapshot = null;
        nameInput.value = draft.name;
        durationInput.value = draft.duration.toFixed(2);
        loopInput.checked = !!draft.loop;
        updatePathOptions();
        renderTracks();
        renderEvents();
        setEditorStatus(dirty ? 'Draft loaded.' : 'Editor synced.', dirty ? 'is-dirty' : 'is-clean');
    }
    function readMetaIntoDraft() {
        if (!draft) return;
        draft.name = (nameInput.value || '').trim() || 'Untitled';
        draft.duration = Math.max(0.5, parseTime(durationInput.value, draft.duration));
        draft.loop = !!loopInput.checked;
        durationInput.value = draft.duration.toFixed(2);
    }
    function upsertTrackKey(path, timeSeconds, value, ease) {
        if (!draft) return;
        var track = null;
        for (var i = 0; i < draft.tracks.length; i++) {
            if (draft.tracks[i].path === path) {
                track = draft.tracks[i];
                break;
            }
        }
        if (!track) {
            track = { path: path, compile: false, keys: [] };
            draft.tracks.push(track);
        }
        var t = Math.max(0, parseFloat(timeSeconds) || 0);
        var replaced = false;
        for (var k = 0; k < track.keys.length; k++) {
            if (Math.abs(track.keys[k].t - t) <= 1e-4) {
                track.keys[k] = { t: t, v: value, ease: normalizeEase(ease) };
                replaced = true;
                break;
            }
        }
        if (!replaced) {
            track.keys.push({ t: t, v: value, ease: normalizeEase(ease) });
        }
        track.keys.sort(function(a, b) { return a.t - b.t; });
    }
    function captureControlsSnapshot() {
        var out = {};

        function addValue(path, value) {
            if (typeof value === 'number') {
                if (!isFinite(value)) return;
                out[path] = value;
                return;
            }
            if (typeof value === 'boolean' || typeof value === 'string') {
                out[path] = value;
            }
        }

        function walk(prefix, node, depth) {
            if (!node || typeof node !== 'object' || depth > 8) return;
            if (Array.isArray(node)) return;
            var keys = Object.keys(node);
            for (var i = 0; i < keys.length; i++) {
                var key = keys[i];
                var value = node[key];
                if (typeof value === 'function') continue;
                var path = prefix ? (prefix + '.' + key) : key;
                if (value && typeof value === 'object') {
                    if (Array.isArray(value)) continue;
                    walk(path, value, depth + 1);
                } else {
                    addValue(path, value);
                }
            }
        }

        if (typeof shader !== 'undefined' && shader && shader.parameters) {
            walk('', shader.parameters, 0);
        }
        if (typeof cameraPan !== 'undefined' && cameraPan) {
            addValue('cameraPan.x', cameraPan.x);
            addValue('cameraPan.y', cameraPan.y);
        }
        if (typeof camera !== 'undefined' && camera) {
            if (camera.position) {
                addValue('camera.position.x', camera.position.x);
                addValue('camera.position.y', camera.position.y);
                addValue('camera.position.z', camera.position.z);
            }
            if (camera.quaternion) {
                addValue('camera.quaternion.x', camera.quaternion.x);
                addValue('camera.quaternion.y', camera.quaternion.y);
                addValue('camera.quaternion.z', camera.quaternion.z);
                addValue('camera.quaternion.w', camera.quaternion.w);
            }
        }

        return out;
    }
    function addAutoKeyframeFromControls() {
        if (!draft) return;
        readMetaIntoDraft();

        var keyTime = parseTime(keyTimeInput.value, currentTime());
        keyTimeInput.value = keyTime.toFixed(2);
        var currentSnapshot = captureControlsSnapshot();
        var currentPaths = Object.keys(currentSnapshot);
        if (!currentPaths.length) {
            setEditorStatus('No controllable values were found for auto keyframing.', 'is-warning');
            return;
        }

        var previous = autoKeyframeState.lastSnapshot;
        if (!previous) {
            autoKeyframeState.lastSnapshot = {
                time: keyTime,
                values: clonePlain(currentSnapshot)
            };
            setEditorStatus(
                'Auto-key baseline captured at ' + keyTime.toFixed(2) +
                's. Change controls, then press AUTO KEYFRAME again.',
                'is-clean'
            );
            return;
        }

        var changedCount = 0;
        for (var i = 0; i < currentPaths.length; i++) {
            var path = currentPaths[i];
            if (!Object.prototype.hasOwnProperty.call(previous.values, path)) continue;
            var prevValue = previous.values[path];
            var currValue = currentSnapshot[path];
            if (areValuesEqual(prevValue, currValue)) continue;
            upsertTrackKey(path, previous.time, prevValue, 'linear');
            upsertTrackKey(path, keyTime, currValue, 'linear');
            changedCount += 1;
        }

        // If some paths disappeared from snapshot, ignore them (runtime shape changes).
        draft.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });
        draft.duration = Math.max(draft.duration, keyTime, previous.time);
        durationInput.value = draft.duration.toFixed(2);

        autoKeyframeState.lastSnapshot = {
            time: keyTime,
            values: clonePlain(currentSnapshot)
        };

        updatePathOptions();
        renderTracks();

        if (changedCount <= 0) {
            setEditorStatus(
                'No control changes detected since last auto keyframe. Baseline moved to ' +
                keyTime.toFixed(2) + 's.',
                'is-warning'
            );
            return;
        }

        markDirty(
            'Auto keyframe captured ' + changedCount +
            ' changed setting' + (changedCount === 1 ? '' : 's') +
            ' at ' + keyTime.toFixed(2) + 's.'
        );
    }
    function updatePathOptions() {
        if (!pathList) return;
        var map = {};
        var html = '';
        function addPath(path) {
            var clean = (path || '').trim();
            if (!clean || map[clean]) return;
            map[clean] = true;
            html += '<option value="' + clean.replace(/"/g, '&quot;') + '"></option>';
        }
        for (var i = 0; i < PRESENTATION_EDITOR_COMMON_PATHS.length; i++) addPath(PRESENTATION_EDITOR_COMMON_PATHS[i]);
        if (draft && Array.isArray(draft.tracks)) {
            for (var t = 0; t < draft.tracks.length; t++) addPath(draft.tracks[t].path);
        }
        pathList.innerHTML = html;
    }
    function summarizeEvent(ev) {
        if (ev.action === 'set') return (ev.path || '?') + ' = ' + formatValue(ev.value);
        if (ev.action === 'annotation') return 'annotation: ' + ((ev.note && (ev.note.title || ev.note.text || ev.note.body)) || 'note');
        return ev.action;
    }
    function renderTracks() {
        if (!trackList || !draft) return;
        if (!draft.tracks.length) {
            trackList.innerHTML = '<div class="presentation-editor-empty">No tracks yet.</div>';
            return;
        }
        var html = '';
        for (var i = 0; i < draft.tracks.length; i++) {
            var tr = draft.tracks[i];
            html += '<div class="presentation-editor-item"><div class="presentation-editor-item-head">' +
                '<span class="presentation-editor-item-title">' + esc(tr.path) + '</span>' +
                '<button type="button" class="presentation-mini-btn" data-action="del-track" data-ti="' + i + '">DEL TRACK</button>' +
                '</div>';
            for (var k = 0; k < tr.keys.length; k++) {
                var key = tr.keys[k];
                html += '<div class="presentation-editor-entry">' +
                    '<span class="presentation-editor-entry-time">' + key.t.toFixed(2) + 's</span>' +
                    '<span class="presentation-editor-entry-value">' + esc(formatValue(key.v)) + '</span>' +
                    '<span class="presentation-editor-entry-ease">' + esc(key.ease) + '</span>' +
                    '<button type="button" class="presentation-mini-btn" data-action="go-key" data-ti="' + i + '" data-ki="' + k + '">GO</button>' +
                    '<button type="button" class="presentation-mini-btn" data-action="edit-key" data-ti="' + i + '" data-ki="' + k + '">EDIT</button>' +
                    '<button type="button" class="presentation-mini-btn" data-action="del-key" data-ti="' + i + '" data-ki="' + k + '">DEL</button>' +
                    '</div>';
            }
            html += '</div>';
        }
        trackList.innerHTML = html;
    }
    function renderEvents() {
        if (!eventList || !draft) return;
        if (!draft.events.length) {
            eventList.innerHTML = '<div class="presentation-editor-empty">No events yet.</div>';
            return;
        }
        var html = '';
        for (var i = 0; i < draft.events.length; i++) {
            var ev = draft.events[i];
            html += '<div class="presentation-editor-entry">' +
                '<span class="presentation-editor-entry-time">' + ev.t.toFixed(2) + 's</span>' +
                '<span class="presentation-editor-entry-value">' + esc(summarizeEvent(ev)) + '</span>' +
                '<button type="button" class="presentation-mini-btn" data-action="go-event" data-ei="' + i + '">GO</button>' +
                '<button type="button" class="presentation-mini-btn" data-action="edit-event" data-ei="' + i + '">EDIT</button>' +
                '<button type="button" class="presentation-mini-btn" data-action="del-event" data-ei="' + i + '">DEL</button>' +
                '</div>';
        }
        eventList.innerHTML = html;
    }
    function updateActionVisibility() {
        var action = eventActionSelect.value;
        eventSetFields.style.display = (action === 'set') ? '' : 'none';
        eventAnnotationFields.style.display = (action === 'annotation') ? '' : 'none';
    }

    if (toggleBtn) {
        toggleBtn.addEventListener('click', function() {
            var willOpen = !root.classList.contains('is-open');
            setOpen(willOpen);
        });
    }

    newBtn.addEventListener('click', function() {
        if (draftDirty && !window.confirm('Discard draft changes?')) return;
        setDraft(createEmptyTimeline(), false);
        setVisible(true);
        setOpen(true);
    });
    fromLoadedBtn.addEventListener('click', function() {
        if (draftDirty && !window.confirm('Replace draft with loaded timeline?')) return;
        if (hooks && typeof hooks.loadSelectedPreset === 'function' && !getRuntimeTimeline()) {
            hooks.loadSelectedPreset();
        }
        var rt = getRuntimeTimeline();
        if (!rt) {
            setEditorStatus('No loaded timeline found.', 'is-warning');
            return;
        }
        setDraft(rt, false);
        setVisible(true);
        setOpen(true);
    });
    applyBtn.addEventListener('click', function() {
        if (!draft) return;
        readMetaIntoDraft();
        var normalized = normalizeTimeline(draft);
        if (typeof setPresentationTimeline !== 'function' || !setPresentationTimeline(normalized)) {
            setEditorStatus('Failed to apply timeline.', 'is-error');
            return;
        }
        draft = normalized;
        draftDirty = false;
        autoKeyframeState.lastSnapshot = null;
        setEditorStatus('Draft applied to runtime timeline.', 'is-clean');
        setMainStatus('Loaded: ' + (normalized.name || 'timeline'), '');
        if (hooks && typeof hooks.syncFromState === 'function') hooks.syncFromState();
        if (typeof CustomEvent === 'function') {
            section.dispatchEvent(new CustomEvent('presentation:timeline-loaded', { detail: { source: 'editor' } }));
        }
    });
    if (turbulenceLoopEnabledInput) {
        turbulenceLoopEnabledInput.addEventListener('change', function() {
            var enabled = !!turbulenceLoopEnabledInput.checked;
            setPathValue('turbulence_loop_enabled', enabled);
            if (turbulenceLoopSecondsInput) {
                turbulenceLoopSecondsInput.disabled = !enabled;
            }
            if (typeof refreshAllControllersGlobal === 'function') {
                refreshAllControllersGlobal();
            }
        });
    }
    if (turbulenceLoopSecondsInput) {
        turbulenceLoopSecondsInput.addEventListener('change', function() {
            var seconds = clamp(parseFloat(turbulenceLoopSecondsInput.value) || 20.0, 0.25, 240.0);
            turbulenceLoopSecondsInput.value = seconds.toFixed(2);
            setPathValue('turbulence_loop_seconds', seconds);
            if (typeof refreshAllControllersGlobal === 'function') {
                refreshAllControllersGlobal();
            }
        });
    }

    nameInput.addEventListener('input', function() {
        if (!draft) return;
        readMetaIntoDraft();
        markDirty();
    });
    durationInput.addEventListener('change', function() {
        if (!draft) return;
        readMetaIntoDraft();
        markDirty();
    });
    loopInput.addEventListener('change', function() {
        if (!draft) return;
        readMetaIntoDraft();
        markDirty();
    });
    keyTimeNowBtn.addEventListener('click', function() {
        keyTimeInput.value = currentTime().toFixed(2);
    });
    eventTimeNowBtn.addEventListener('click', function() {
        eventTimeInput.value = currentTime().toFixed(2);
    });
    keyCaptureBtn.addEventListener('click', function() {
        var path = (trackPathInput.value || '').trim();
        if (!path) {
            setEditorStatus('Path is required.', 'is-warning');
            return;
        }
        var v = getPathValue(path);
        if (typeof v === 'undefined') {
            setEditorStatus('Unknown path: ' + path, 'is-warning');
            return;
        }
        keyValueInput.value = formatValue(v);
        setEditorStatus('Captured value for ' + path + '.', 'is-clean');
    });
    if (autoKeyframeBtn) {
        autoKeyframeBtn.addEventListener('click', function() {
            addAutoKeyframeFromControls();
        });
    }
    if (trackListExpandBtn) {
        trackListExpandBtn.addEventListener('click', function() {
            setTrackListExpanded(!trackListExpanded);
        });
    }
    keyAddBtn.addEventListener('click', function() {
        if (!draft) return;
        readMetaIntoDraft();
        var path = (trackPathInput.value || '').trim();
        if (!path) {
            setEditorStatus('Path is required.', 'is-warning');
            return;
        }
        var t = parseTime(keyTimeInput.value, currentTime());
        var ease = normalizeEase(keyEaseSelect.value);
        var value = parseValue(keyValueInput.value);
        var track = null;
        for (var i = 0; i < draft.tracks.length; i++) {
            if (draft.tracks[i].path === path) track = draft.tracks[i];
        }
        if (!track) {
            track = { path: path, compile: false, keys: [] };
            draft.tracks.push(track);
        }
        var updated = false;
        for (var k = 0; k < track.keys.length; k++) {
            if (Math.abs(track.keys[k].t - t) <= 1e-4) {
                track.keys[k] = { t: t, v: value, ease: ease };
                updated = true;
            }
        }
        if (!updated) track.keys.push({ t: t, v: value, ease: ease });
        track.keys.sort(function(a, b) { return a.t - b.t; });
        draft.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });
        draft.duration = Math.max(draft.duration, t);
        durationInput.value = draft.duration.toFixed(2);
        updatePathOptions();
        renderTracks();
        markDirty(updated ? 'Keyframe updated in draft.' : 'Keyframe added to draft.');
    });
    keyRemoveBtn.addEventListener('click', function() {
        if (!draft) return;
        var path = (trackPathInput.value || '').trim();
        if (!path) {
            setEditorStatus('Path is required.', 'is-warning');
            return;
        }
        var t = parseTime(keyTimeInput.value, currentTime());
        var ti = -1;
        for (var i = 0; i < draft.tracks.length; i++) {
            if (draft.tracks[i].path === path) ti = i;
        }
        if (ti < 0) {
            setEditorStatus('Track not found.', 'is-warning');
            return;
        }
        var keys = draft.tracks[ti].keys;
        var best = -1;
        var bestDiff = Infinity;
        for (var k = 0; k < keys.length; k++) {
            var diff = Math.abs(keys[k].t - t);
            if (diff < bestDiff) {
                bestDiff = diff;
                best = k;
            }
        }
        if (best < 0 || bestDiff > 0.1) {
            setEditorStatus('No keyframe near selected time.', 'is-warning');
            return;
        }
        keys.splice(best, 1);
        if (!keys.length) draft.tracks.splice(ti, 1);
        updatePathOptions();
        renderTracks();
        markDirty('Keyframe removed from draft.');
    });
    eventActionSelect.addEventListener('change', updateActionVisibility);
    eventAddBtn.addEventListener('click', function() {
        if (!draft) return;
        readMetaIntoDraft();
        var t = parseTime(eventTimeInput.value, currentTime());
        var action = eventActionSelect.value;
        var ev = { t: t, action: action };
        if (action === 'set') {
            var path = (eventPathInput.value || '').trim();
            if (!path) {
                setEditorStatus('Set-event path is required.', 'is-warning');
                return;
            }
            ev.path = path;
            ev.value = parseValue(eventValueInput.value);
            ev.compile = !!eventCompileInput.checked;
        } else if (action === 'annotation') {
            ev.note = {
                title: (noteTitleInput.value || '').trim() || 'Note',
                text: (noteTextInput.value || '').trim() || '...',
                placement: notePlacementSelect.value || 'auto',
                anchor: { mode: 'world', target: noteTargetSelect.value || 'black_hole' }
            };
        }
        draft.events.push(ev);
        draft.events.sort(function(a, b) { return a.t - b.t; });
        draft.duration = Math.max(draft.duration, t);
        durationInput.value = draft.duration.toFixed(2);
        updatePathOptions();
        renderEvents();
        markDirty('Event added to draft.');
    });

    trackList.addEventListener('click', function(event) {
        var btn = event.target;
        if (!btn || btn.tagName !== 'BUTTON') return;
        if (!draft) return;
        var action = btn.getAttribute('data-action');
        var ti = parseInt(btn.getAttribute('data-ti'), 10);
        var ki = parseInt(btn.getAttribute('data-ki'), 10);
        var track = isFinite(ti) ? draft.tracks[ti] : null;
        var key = track && isFinite(ki) ? track.keys[ki] : null;
        if (action === 'go-key' && key) goToTime(key.t);
        if (action === 'edit-key' && key) {
            trackPathInput.value = track.path;
            keyTimeInput.value = key.t.toFixed(2);
            keyEaseSelect.value = normalizeEase(key.ease);
            keyValueInput.value = formatValue(key.v);
        }
        if (action === 'del-key' && key) {
            track.keys.splice(ki, 1);
            if (!track.keys.length) draft.tracks.splice(ti, 1);
            renderTracks();
            markDirty('Keyframe removed from draft.');
        }
        if (action === 'del-track' && track) {
            draft.tracks.splice(ti, 1);
            renderTracks();
            markDirty('Track removed from draft.');
        }
    });
    eventList.addEventListener('click', function(event) {
        var btn = event.target;
        if (!btn || btn.tagName !== 'BUTTON') return;
        if (!draft) return;
        var action = btn.getAttribute('data-action');
        var ei = parseInt(btn.getAttribute('data-ei'), 10);
        var ev = isFinite(ei) ? draft.events[ei] : null;
        if (!ev) return;
        if (action === 'go-event') goToTime(ev.t);
        if (action === 'edit-event') {
            eventTimeInput.value = ev.t.toFixed(2);
            eventActionSelect.value = ev.action || 'set';
            updateActionVisibility();
            if (ev.action === 'set') {
                eventPathInput.value = ev.path || '';
                eventValueInput.value = formatValue(ev.value);
                eventCompileInput.checked = !!ev.compile;
            } else if (ev.action === 'annotation') {
                var note = ev.note || {};
                noteTitleInput.value = note.title || '';
                noteTextInput.value = note.text || note.body || '';
                noteTargetSelect.value = (note.anchor && note.anchor.target) || 'black_hole';
                notePlacementSelect.value = note.placement || 'auto';
            }
        }
        if (action === 'del-event') {
            draft.events.splice(ei, 1);
            renderEvents();
            markDirty('Event removed from draft.');
        }
    });

    jsonExportBtn.addEventListener('click', function() {
        if (!draft) return;
        jsonTextarea.value = JSON.stringify(normalizeTimeline(draft), null, 2);
        setEditorStatus('Draft exported to JSON text area.', 'is-clean');
    });
    if (jsonImportFileBtn) {
        jsonImportFileBtn.addEventListener('click', function() {
            if (!jsonFileInput) {
                setEditorStatus('File import is unavailable. Use IMPORT TEXT.', 'is-warning');
                return;
            }
            // Re-selecting the same file should still trigger change.
            jsonFileInput.value = '';
            jsonFileInput.click();
        });
    }
    if (jsonFileInput) {
        jsonFileInput.addEventListener('change', function() {
            var file = (jsonFileInput.files && jsonFileInput.files[0]) ? jsonFileInput.files[0] : null;
            if (!file) return;

            if (typeof file.text === 'function') {
                file.text()
                    .then(function(text) {
                        importJsonTextIntoDraft(text, 'file "' + file.name + '"');
                    })
                    .catch(function(err) {
                        setEditorStatus('Failed to read file: ' + (err && err.message ? err.message : err), 'is-error');
                    });
                return;
            }

            if (typeof FileReader === 'undefined') {
                setEditorStatus('FileReader API is unavailable. Use IMPORT TEXT instead.', 'is-warning');
                return;
            }

            var reader = new FileReader();
            reader.onload = function() {
                importJsonTextIntoDraft(String(reader.result || ''), 'file "' + file.name + '"');
            };
            reader.onerror = function() {
                setEditorStatus('Failed to read file "' + file.name + '".', 'is-error');
            };
            reader.readAsText(file);
        });
    }
    jsonImportBtn.addEventListener('click', function() {
        importJsonTextIntoDraft(jsonTextarea.value, 'text area');
    });
    jsonDownloadBtn.addEventListener('click', function() {
        if (!draft) return;
        var normalized = normalizeTimeline(draft);
        var filename = ((normalized.name || 'presentation')
            .toLowerCase()
            .replace(/[^a-z0-9\-_.]+/g, '-') || 'presentation') + '.json';
        var blob = new Blob([JSON.stringify(normalized, null, 2)], { type: 'application/json' });
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        setTimeout(function() { URL.revokeObjectURL(url); }, 1000);
        setEditorStatus('Draft downloaded as ' + filename + '.', 'is-clean');
    });

    section.addEventListener('presentation:timeline-loaded', function() {
        if (draftDirty) return;
        var runtime = getRuntimeTimeline();
        if (runtime) setDraft(runtime, false);
    });

    function syncFromRuntime(force) {
        syncTurbulenceLoopInputs();
        if (draftDirty && !force) return;
        var runtime = getRuntimeTimeline();
        if (runtime) {
            setDraft(runtime, false);
            return;
        }
        setDraft(createEmptyTimeline(), false);
        setEditorStatus('No loaded timeline. Create a draft or load a preset.', 'is-warning');
    }

    updateActionVisibility();
    keyTimeInput.value = currentTime().toFixed(2);
    eventTimeInput.value = currentTime().toFixed(2);
    trackPathInput.value = 'observer.distance';
    eventPathInput.value = 'accretion_mode';
    eventValueInput.value = 'thin';
    syncTurbulenceLoopInputs();
    setTrackListExpanded(false);
    setOpen(false);
    setVisible(false);

    return {
        syncFromRuntime: syncFromRuntime,
        setVisible: setVisible,
        setOpen: setOpen,
        startNewPresetDraft: startNewPresetDraft
    };
}
