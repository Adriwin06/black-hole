// ═══════════════════════════════════════════════════════════════════════════════
// Timeline Panel — bottom-docked dopesheet for animation editing
// ═══════════════════════════════════════════════════════════════════════════════
// Provides a full-width bottom panel similar to After Effects / Blender dopesheet
// with transport controls, track list, keyframe lanes, ruler, and key inspector.
// Public API exposed via global `timelinePanelBinding`.

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

var timelinePanelBinding = null;

function buildTimelinePanel() {
    'use strict';

    // ── Utility helpers ─────────────────────────────────────────────────────
    function clamp(v, lo, hi) { return Math.max(lo, Math.min(hi, v)); }
    function esc(s) { return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
    function clonePlain(o) { return JSON.parse(JSON.stringify(o)); }
    function parseTime(v, fallback) { var n = parseFloat(v); return isFinite(n) ? Math.max(0, n) : fallback; }
    function normalizeEase(e) {
        if (e === 'smooth' || e === 'smoother') return e;
        return 'linear';
    }
    function formatValue(v) {
        if (typeof v === 'number') return isFinite(v) ? v.toPrecision(8) : String(v);
        if (typeof v === 'boolean') return v ? 'true' : 'false';
        if (typeof v === 'string') return v;
        if (v === null || typeof v === 'undefined') return '';
        return JSON.stringify(v);
    }
    function parseValue(s) {
        if (s === 'true') return true;
        if (s === 'false') return false;
        var n = Number(s);
        if (s !== '' && isFinite(n)) return n;
        return s;
    }
    function areValuesEqual(a, b) {
        if (a === b) return true;
        if (typeof a === 'number' && typeof b === 'number') return Math.abs(a - b) < 1e-9;
        return false;
    }

    // ── Build DOM ───────────────────────────────────────────────────────────
    var panel = document.createElement('div');
    panel.id = 'tl-panel';
    panel.className = 'tl-panel tl-panel--collapsed';
    panel.innerHTML =
        '<div id="tl-resize-handle" class="tl-resize-handle"></div>' +

        // ── Transport bar ──
        '<div class="tl-transport">' +
            '<div class="tl-transport-left">' +
                '<button id="tl-btn-play" class="tl-btn tl-btn--accent" type="button" title="Play">&#9654;</button>' +
                '<button id="tl-btn-pause" class="tl-btn" type="button" title="Pause">&#10074;&#10074;</button>' +
                '<button id="tl-btn-stop" class="tl-btn" type="button" title="Stop">&#9632;</button>' +
                '<span class="tl-transport-sep"></span>' +
                '<span class="tl-time-display">' +
                    '<input id="tl-time-current" class="tl-time-input" type="number" min="0" step="0.01" value="0" title="Current time (editable)">' +
                    '<span class="tl-time-sep">/</span>' +
                    '<input id="tl-time-duration" class="tl-time-input" type="number" min="0.5" step="0.1" value="12" title="Duration (editable)">' +
                '</span>' +
            '</div>' +
            '<div class="tl-transport-center">' +
                '<div id="tl-scrubber" class="tl-scrubber">' +
                    '<div id="tl-scrub-fill" class="tl-scrub-fill"></div>' +
                    '<div id="tl-scrub-head" class="tl-scrub-head"></div>' +
                '</div>' +
            '</div>' +
            '<div class="tl-transport-right">' +
                '<label for="tl-preset-select" class="tl-preset-label">Preset</label>' +
                '<select id="tl-preset-select" class="tl-preset-select"></select>' +
                '<span class="tl-transport-sep"></span>' +
                '<button id="tl-btn-motion" class="tl-btn tl-btn--motion" type="button" title="Insert a predefined motion function">⊕&nbsp;FX</button>' +
                '<button id="tl-btn-rec" class="tl-btn tl-btn--rec" type="button" title="Recording settings">&#9679;&nbsp;REC</button>' +
                '<span class="tl-transport-sep"></span>' +
                '<button id="tl-btn-auto-key" class="tl-btn tl-btn--warn" type="button" title="Auto Keyframe: capture changes">AUTO KEY</button>' +
                '<button id="tl-btn-add-track" class="tl-btn" type="button" title="Add a new track">+ TRACK</button>' +
                '<span class="tl-transport-sep"></span>' +
                '<button id="tl-btn-import" class="tl-btn" type="button" title="Import timeline from JSON file">&#8593; IMPORT</button>' +
                '<button id="tl-btn-export" class="tl-btn" type="button" title="Export timeline as JSON file">&#8595; EXPORT</button>' +
                '<span class="tl-transport-sep"></span>' +
                '<button id="tl-btn-close" class="tl-btn" type="button" title="Close panel">&times;</button>' +
            '</div>' +
        '</div>' +

        // ── 3-column body ──
        '<div class="tl-body">' +

            // Left: track list
            '<div class="tl-tracks">' +
                '<div class="tl-tracks-head">' +
                    '<span class="tl-tracks-title">TRACKS</span>' +
                '</div>' +
                '<div id="tl-track-list" class="tl-track-list"></div>' +
            '</div>' +

            // Center: dopesheet lanes
            '<div class="tl-dopesheet">' +
                '<div id="tl-ruler" class="tl-ruler"></div>' +
                '<div id="tl-lanes" class="tl-lanes"></div>' +
            '</div>' +

            // Right: key inspector
            '<div class="tl-inspector">' +
                '<div class="tl-inspector-title">KEY INSPECTOR</div>' +
                '<div id="tl-insp-summary" class="tl-insp-summary">No keyframe selected.</div>' +
                '<div class="tl-insp-row">' +
                    '<label>Path</label>' +
                    '<input id="tl-insp-path" list="tl-path-datalist" type="text" placeholder="observer.distance">' +
                '</div>' +
                '<div class="tl-insp-row">' +
                    '<label>Time</label>' +
                    '<input id="tl-insp-time" type="number" min="0" step="0.01" value="0">' +
                    '<select id="tl-insp-ease">' +
                        '<option value="linear">linear</option>' +
                        '<option value="smooth">smooth</option>' +
                        '<option value="smoother">smoother</option>' +
                    '</select>' +
                '</div>' +
                '<div class="tl-insp-row">' +
                    '<label>Value</label>' +
                    '<input id="tl-insp-value" type="text" placeholder="11 or true">' +
                '</div>' +
                '<div class="tl-insp-actions">' +
                    '<button id="tl-insp-use-time" class="tl-mini-btn" type="button">USE TIME</button>' +
                    '<button id="tl-insp-capture" class="tl-mini-btn" type="button">CAPTURE</button>' +
                '</div>' +
                '<div class="tl-insp-actions">' +
                    '<button id="tl-insp-set" class="tl-mini-btn tl-mini-btn--accent" type="button">SET KEY</button>' +
                    '<button id="tl-insp-del" class="tl-mini-btn tl-mini-btn--danger" type="button">DELETE KEY</button>' +
                '</div>' +
                '<datalist id="tl-path-datalist"></datalist>' +
            '</div>' +

        '</div>' +

        // ── Motion Functions modal (floats above panel) ──
        '<div id="tl-motion-modal" class="tl-motion-modal">' +
            '<div class="tl-motion-hdr">' +
                '<span class="tl-motion-title">MOTION FUNCTION</span>' +
                '<button id="tl-motion-close" class="tl-btn" type="button" title="Close">&times;</button>' +
            '</div>' +
            '<div class="tl-motion-row tl-motion-type-row">' +
                '<label>Type</label>' +
                '<select id="tl-motion-type" class="tl-motion-input">' +
                    '<option value="orbit">Orbit around&nbsp;BH</option>' +
                    '<option value="zoom">Zoom in / out</option>' +
                    '<option value="exposure">Exposure fade</option>' +
                    '<option value="inclination">Inclination sweep</option>' +
                '</select>' +
            '</div>' +
            '<div id="tl-motion-params" class="tl-motion-params"></div>' +
            '<div class="tl-motion-footer">' +
                '<button id="tl-motion-apply" class="tl-mini-btn tl-mini-btn--accent" type="button">APPLY</button>' +
            '</div>' +
        '</div>' +

        // ── REC modal (floats above panel) ──
        '<div id="tl-rec-modal" class="tl-rec-modal">' +
            '<div class="tl-rec-hdr">' +
                '<span class="tl-rec-title">RECORDING</span>' +
                '<button id="tl-rec-close" class="tl-btn" type="button" title="Close">&times;</button>' +
            '</div>' +
            '<div class="tl-rec-body">' +
                '<div class="tl-rec-checks">' +
                    '<label class="tl-rec-check"><input type="checkbox" id="tl-rec-loop"> Loop timeline</label>' +
                    '<label class="tl-rec-check"><input type="checkbox" id="tl-rec-annotations" checked> Show explanations</label>' +
                    '<label class="tl-rec-check"><input type="checkbox" id="tl-rec-annot-record"> Include text in recording</label>' +
                    '<label class="tl-rec-check"><input type="checkbox" id="tl-rec-reset-sim" checked> Reset simulation on rec start</label>' +
                '</div>' +
                '<div class="tl-rec-row">' +
                    '<label>Quality</label>' +
                    '<select id="tl-rec-quality" class="tl-rec-select"></select>' +
                '</div>' +
                '<div class="tl-rec-row">' +
                    '<label>Mode</label>' +
                    '<select id="tl-rec-mode" class="tl-rec-select"></select>' +
                '</div>' +
                '<div class="tl-rec-row">' +
                    '<label>Resolution</label>' +
                    '<select id="tl-rec-resolution" class="tl-rec-select"></select>' +
                '</div>' +
                '<div class="tl-rec-row">' +
                    '<label>FPS</label>' +
                    '<input id="tl-rec-fps" class="tl-rec-num" type="number" min="24" max="120" step="1" value="60">' +
                    '<label style="margin-left:8px">Mbps</label>' +
                    '<input id="tl-rec-bitrate" class="tl-rec-num" type="number" min="4" max="80" step="1" value="20">' +
                '</div>' +
            '</div>' +
            '<div class="tl-rec-actions">' +
                '<button id="tl-rec-start" class="tl-mini-btn tl-mini-btn--accent" type="button">&#9679; START REC</button>' +
                '<button id="tl-rec-stop" class="tl-mini-btn tl-mini-btn--danger" type="button">&#9632; STOP REC</button>' +
            '</div>' +
            '<div id="tl-rec-status" class="tl-rec-status">Idle</div>' +
        '</div>' +

        // ── Status bar ──
        '<div id="tl-status" class="tl-status"></div>';

    document.body.appendChild(panel);

    // ── Element refs ──
    var playBtn      = panel.querySelector('#tl-btn-play');
    var pauseBtn     = panel.querySelector('#tl-btn-pause');
    var stopBtn      = panel.querySelector('#tl-btn-stop');
    var autoKeyBtn   = panel.querySelector('#tl-btn-auto-key');
    var addTrackBtn  = panel.querySelector('#tl-btn-add-track');
    var importBtn    = panel.querySelector('#tl-btn-import');
    var exportBtn    = panel.querySelector('#tl-btn-export');
    var closeBtn     = panel.querySelector('#tl-btn-close');
    var timeCurrent  = panel.querySelector('#tl-time-current');
    var timeDuration = panel.querySelector('#tl-time-duration');
    var scrubber    = panel.querySelector('#tl-scrubber');
    var scrubFill   = panel.querySelector('#tl-scrub-fill');
    var scrubHead   = panel.querySelector('#tl-scrub-head');
    var presetSelect= panel.querySelector('#tl-preset-select');
    var trackListEl = panel.querySelector('#tl-track-list');
    var rulerEl     = panel.querySelector('#tl-ruler');
    var lanesEl     = panel.querySelector('#tl-lanes');
    var inspPath    = panel.querySelector('#tl-insp-path');
    var inspTime    = panel.querySelector('#tl-insp-time');
    var inspEase    = panel.querySelector('#tl-insp-ease');
    var inspValue   = panel.querySelector('#tl-insp-value');
    var inspSummary = panel.querySelector('#tl-insp-summary');
    var inspUseTime = panel.querySelector('#tl-insp-use-time');
    var inspCapture = panel.querySelector('#tl-insp-capture');
    var inspSet     = panel.querySelector('#tl-insp-set');
    var inspDel     = panel.querySelector('#tl-insp-del');
    var statusEl    = panel.querySelector('#tl-status');
    var pathDatalist= panel.querySelector('#tl-path-datalist');
    var resizeHandle= panel.querySelector('#tl-resize-handle');
    var motionBtn      = panel.querySelector('#tl-btn-motion');
    var motionModal    = panel.querySelector('#tl-motion-modal');
    var motionCloseBtn = panel.querySelector('#tl-motion-close');
    var motionTypeEl   = panel.querySelector('#tl-motion-type');
    var motionParamsEl = panel.querySelector('#tl-motion-params');
    var motionApplyBtn = panel.querySelector('#tl-motion-apply');
    var recBtn           = panel.querySelector('#tl-btn-rec');
    var recModal         = panel.querySelector('#tl-rec-modal');
    var recCloseBtn      = panel.querySelector('#tl-rec-close');
    var recLoopCb        = panel.querySelector('#tl-rec-loop');
    var recAnnotCb       = panel.querySelector('#tl-rec-annotations');
    var recAnnotRecordCb = panel.querySelector('#tl-rec-annot-record');
    var recResetSimCb    = panel.querySelector('#tl-rec-reset-sim');
    var recQualitySelect = panel.querySelector('#tl-rec-quality');
    var recModeSelect    = panel.querySelector('#tl-rec-mode');
    var recResSelect     = panel.querySelector('#tl-rec-resolution');
    var recFpsInput      = panel.querySelector('#tl-rec-fps');
    var recBitrateInput  = panel.querySelector('#tl-rec-bitrate');
    var recStartBtn      = panel.querySelector('#tl-rec-start');
    var recStopBtn       = panel.querySelector('#tl-rec-stop');
    var recStatusEl      = panel.querySelector('#tl-rec-status');

    // ── State ───────────────────────────────────────────────────────────────
    var draft          = null;
    var selectedTrack  = '';
    var selectedKeyT   = NaN;
    // Multi-select: array of { path: string, t: number }
    var selectedKeys   = [];
    var panelOpen      = false;
    var autoKeySnapshot= null;
    var syncTimer      = null;
    var PANEL_HEIGHT_KEY  = 'black-hole.tl-panel.height';
    var PANEL_STATE_KEY   = 'black-hole.tl-panel.state';
    var PANEL_DEFAULT_H   = 320;

    // ── Clipboard (copy/paste keyframes) ───────────────────────────────────
    // Array of { path, t, v, ease } with anchorT = min(t) in the set
    var clipboard = null;

    // ── Key drag state ──────────────────────────────────────────────────────
    // Set when the user starts dragging a diamond; reset on pointerup/cancel
    var keyDragState = null;

    // ── Undo/Redo history ───────────────────────────────────────────────────
    var undoStack = [];
    var redoStack = [];
    var UNDO_MAX  = 40;
    function pushUndo() {
        if (!draft) return;
        undoStack.push(clonePlain(draft));
        if (undoStack.length > UNDO_MAX) undoStack.shift();
        redoStack = [];
    }
    function undo() {
        if (!undoStack.length) { setStatus('Nothing to undo.', 'tl-status--warn'); return; }
        redoStack.push(clonePlain(draft));
        draft = undoStack.pop();
        applyDraft();
        rebuildAll();
        setStatus('Undo.', '');
    }
    function redo() {
        if (!redoStack.length) { setStatus('Nothing to redo.', 'tl-status--warn'); return; }
        undoStack.push(clonePlain(draft));
        draft = redoStack.pop();
        applyDraft();
        rebuildAll();
        setStatus('Redo.', '');
    }
    var PANEL_MIN_H      = 180;
    var PANEL_MAX_H      = 700;

    // ── Persistence ───────────────────────────────────────────────────────
    function saveState() {
        try {
            var state = {
                preset: presetSelect.value || '',
                draft: draft ? clonePlain(draft) : null,
                selectedTrack: selectedTrack,
                selectedKeys: selectedKeys.slice(),
                wasOpen: panelOpen
            };
            sessionStorage.setItem(PANEL_STATE_KEY, JSON.stringify(state));
        } catch(e) {}
    }
    function loadState() {
        try {
            var raw = sessionStorage.getItem(PANEL_STATE_KEY);
            if (!raw) return null;
            return JSON.parse(raw);
        } catch(e) { return null; }
    }

    // ── Panel open/close ──────────────────────────────────────────────────────
    function updatePushedOffset() {
        var h = panel.getBoundingClientRect().height;
        document.body.style.setProperty('--tl-h', Math.round(h + 10) + 'px');
    }

    function setPanelOpen(open) {
        panelOpen = !!open;
        panel.classList.toggle('tl-panel--collapsed', !panelOpen);
        document.body.classList.toggle('has-timeline-panel', panelOpen);
        if (panelOpen) {
            updatePushedOffset();
            var saved = loadState();
            populatePresets();
            if (saved) {
                // Restore persisted state
                if (saved.preset && presetSelect.querySelector('option[value="' + CSS.escape(saved.preset) + '"]')) {
                    presetSelect.value = saved.preset;
                }
                if (saved.draft) {
                    draft = normalizeTL(saved.draft);
                    rebuildAll();
                } else {
                    syncFromRuntime();
                }
                if (saved.selectedTrack) selectedTrack = saved.selectedTrack;
                if (saved.selectedKeys && saved.selectedKeys.length) selectedKeys = saved.selectedKeys;
            } else {
                // No saved state — default to a fresh empty timeline
                presetSelect.value = '';
                loadPresetByName('');
            }
            updateTimeInputs();
            startSync();
        } else {
            saveState();
            stopSync();
            document.body.style.removeProperty('--tl-h');
        }
    }

    function toggle() { setPanelOpen(!panelOpen); }

    // ── Panel resize ────────────────────────────────────────────────────────
    var storedH = null;
    try { storedH = parseFloat(localStorage.getItem(PANEL_HEIGHT_KEY)); } catch(e) {}
    if (isFinite(storedH)) panel.style.height = clamp(storedH, PANEL_MIN_H, PANEL_MAX_H) + 'px';
    else panel.style.height = PANEL_DEFAULT_H + 'px';

    (function initResize() {
        var dragging = false, startY = 0, startH = 0;
        resizeHandle.addEventListener('pointerdown', function(e) {
            if (e.button !== 0) return;
            dragging = true; startY = e.clientY;
            startH = panel.getBoundingClientRect().height;
            resizeHandle.setPointerCapture(e.pointerId);
            document.body.classList.add('tl-resizing');
            e.preventDefault();
        });
        document.addEventListener('pointermove', function(e) {
            if (!dragging) return;
            var h = clamp(startH - (e.clientY - startY), PANEL_MIN_H, PANEL_MAX_H);
            panel.style.height = h + 'px';
            updatePushedOffset();
            e.preventDefault();
        });
        function end(e) {
            if (!dragging) return;
            dragging = false;
            document.body.classList.remove('tl-resizing');
            try { localStorage.setItem(PANEL_HEIGHT_KEY, String(Math.round(panel.getBoundingClientRect().height))); } catch(ex) {}
        }
        document.addEventListener('pointerup', end);
        document.addEventListener('pointercancel', end);
    })();

    // ── Scrubber ────────────────────────────────────────────────────────────
    function currentTime() {
        if (typeof getPresentationState === 'function') {
            var st = getPresentationState();
            if (st && isFinite(st.time)) return Math.max(0, st.time);
        }
        return 0;
    }
    function getDuration() {
        if (draft && draft.duration > 0) return draft.duration;
        if (typeof getPresentationState === 'function') {
            var st = getPresentationState();
            if (st && st.duration > 0) return st.duration;
        }
        return 12;
    }
    function updateTimeInputs() {
        var t = currentTime(), d = getDuration();
        // Only update if not focused (don't interrupt user typing)
        if (document.activeElement !== timeCurrent)  timeCurrent.value = t.toFixed(2);
        if (document.activeElement !== timeDuration) timeDuration.value = d.toFixed(2);
    }
    function updateScrubber() {
        var d = getDuration(), t = currentTime();
        var pct = clamp(t / Math.max(d, 0.001) * 100, 0, 100);
        scrubFill.style.width = pct + '%';
        scrubHead.style.left = pct + '%';
    }

    (function initScrub() {
        var dragging = false;
        function seek(clientX) {
            var rect = scrubber.getBoundingClientRect();
            if (!rect.width) return;
            var ratio = clamp((clientX - rect.left) / rect.width, 0, 1);
            var t = ratio * getDuration();
            if (typeof seekPresentation === 'function') seekPresentation(t);
            updateTimeInputs(); updateScrubber();
        }
        scrubber.addEventListener('pointerdown', function(e) {
            if (e.button !== 0) return;
            dragging = true;
            scrubber.setPointerCapture(e.pointerId);
            seek(e.clientX); e.preventDefault();
        });
        scrubber.addEventListener('pointermove', function(e) {
            if (dragging) { seek(e.clientX); e.preventDefault(); }
        });
        function endScrub() { dragging = false; }
        scrubber.addEventListener('pointerup', endScrub);
        scrubber.addEventListener('pointercancel', endScrub);
    })();

    // ── Transport ───────────────────────────────────────────────────────────
    playBtn.addEventListener('click', function() {
        // If nothing loaded, load the selected preset first
        if (typeof getPresentationState === 'function') {
            var st = getPresentationState();
            if (!st.loaded && presetSelect.value) {
                loadPresetByName(presetSelect.value);
            }
        }
        if (typeof playPresentation === 'function') playPresentation(false);
    });
    pauseBtn.addEventListener('click', function() {
        if (typeof pausePresentation === 'function') pausePresentation();
    });
    stopBtn.addEventListener('click', function() {
        if (typeof stopPresentation === 'function') stopPresentation();
    });
    closeBtn.addEventListener('click', function() { setPanelOpen(false); });

    // ── Editable time / duration inputs ───────────────────────────────────
    timeCurrent.addEventListener('change', function() {
        var t = parseTime(timeCurrent.value, 0);
        t = clamp(t, 0, getDuration());
        if (typeof seekPresentation === 'function') seekPresentation(t);
        updateTimeInputs(); updateScrubber(); updatePlayheads();
    });
    timeCurrent.addEventListener('keydown', function(e) {
        if (e.key === 'Enter') { timeCurrent.blur(); }
    });
    timeDuration.addEventListener('change', function() {
        var d = parseTime(timeDuration.value, 12);
        d = Math.max(0.5, d);
        if (draft) {
            pushUndo();
            draft.duration = d;
            applyDraft();
            rebuildAll();
            setStatus('Duration set to ' + d.toFixed(2) + 's.', '');
        }
    });
    timeDuration.addEventListener('keydown', function(e) {
        if (e.key === 'Enter') { timeDuration.blur(); }
    });

    // ── Preset selector ─────────────────────────────────────────────────────
    function populatePresets() {
        var names = (typeof listPresentationPresets === 'function')
            ? listPresentationPresets() : [];
        var currentVal = presetSelect.value;
        presetSelect.innerHTML = '';

        // "New (empty)" option
        var emptyOpt = document.createElement('option');
        emptyOpt.value = '';
        emptyOpt.textContent = '— new empty —';
        presetSelect.appendChild(emptyOpt);

        for (var i = 0; i < names.length; i++) {
            var opt = document.createElement('option');
            opt.value = names[i];
            opt.textContent = names[i];
            presetSelect.appendChild(opt);
        }

        // Restore previous selection if possible
        if (currentVal && presetSelect.querySelector('option[value="' + CSS.escape(currentVal) + '"]')) {
            presetSelect.value = currentVal;
        } else if (names.length) {
            // Pick current loaded timeline name or first preset
            var st = typeof getPresentationState === 'function' ? getPresentationState() : null;
            var loadedName = st && st.name ? st.name : '';
            if (loadedName && names.indexOf(loadedName) !== -1) {
                presetSelect.value = loadedName;
            } else {
                var defaultName = (names.indexOf('Full Feature Tour') !== -1) ? 'Full Feature Tour' : names[0];
                presetSelect.value = defaultName;
            }
        }
    }

    function loadPresetByName(name) {
        if (!name) {
            // Empty = start fresh
            draft = normalizeTL({ name: 'Untitled', duration: 12, tracks: [], events: [] });
            if (typeof setPresentationTimeline === 'function') {
                setPresentationTimeline(clonePlain(draft));
            }
            rebuildAll();
            setStatus('New empty timeline created.', '');
            return;
        }
        if (typeof loadPresentationPreset === 'function') {
            if (loadPresentationPreset(name)) {
                syncFromRuntime();
                setStatus('Loaded preset: ' + name, '');
            } else {
                setStatus('Failed to load preset: ' + name, 'tl-status--warn');
            }
        }
    }

    presetSelect.addEventListener('change', function() {
        loadPresetByName(presetSelect.value);
    });

    // ── Sync from runtime ───────────────────────────────────────────────────
    function syncFromRuntime() {
        var rt = null;
        if (typeof getPresentationTimeline === 'function') rt = getPresentationTimeline();
        if (rt) {
            draft = normalizeTL(rt);
            rebuildAll();
        }
        updateTimeInputs();
        updateScrubber();
    }

    function rebuildAll() {
        rebuildTrackList();
        rebuildLanes();
        updateRuler();
        updateInspector();
        updatePathDatalist();
    }

    function startSync() {
        if (syncTimer) return;
        syncTimer = setInterval(function() {
            updateTimeInputs();
            updateScrubber();
            updatePlayheads();
            syncRecModal();
        }, 80);
    }
    function stopSync() { clearInterval(syncTimer); syncTimer = null; }

    // ── Timeline data helpers ───────────────────────────────────────────────
    function normalizeTL(raw) {
        var src = (raw && typeof raw === 'object') ? clonePlain(raw) : {};
        var out = { name: src.name || 'Untitled', duration: Math.max(0.5, parseTime(src.duration, 12)),
                    loop: !!src.loop, tracks: [], events: [] };
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
        out.events.sort(function(a, b) { return a.t - b.t; });
        return out;
    }

    function getTrackByPath(path) {
        if (!draft) return null;
        for (var i = 0; i < draft.tracks.length; i++)
            if (draft.tracks[i].path === path) return draft.tracks[i];
        return null;
    }

    function getKeyAt(track, t) {
        if (!track) return -1;
        for (var i = 0; i < track.keys.length; i++)
            if (Math.abs(track.keys[i].t - t) <= 1e-4) return i;
        return -1;
    }

    function applyDraft() {
        if (!draft) return;
        if (typeof setPresentationTimeline === 'function') {
            setPresentationTimeline(clonePlain(draft));
            setStatus('Timeline applied.', '');
        }
    }

    // ── Path datalist ───────────────────────────────────────────────────────
    function updatePathDatalist() {
        var map = {}, html = '';
        function add(p) { if (!p || map[p]) return; map[p] = 1; html += '<option value="' + esc(p) + '">'; }
        if (typeof PRESENTATION_EDITOR_COMMON_PATHS !== 'undefined') {
            for (var i = 0; i < PRESENTATION_EDITOR_COMMON_PATHS.length; i++) add(PRESENTATION_EDITOR_COMMON_PATHS[i]);
        }
        if (draft) for (var j = 0; j < draft.tracks.length; j++) add(draft.tracks[j].path);
        pathDatalist.innerHTML = html;
    }

    function setStatus(text, cls) {
        statusEl.textContent = text || '';
        statusEl.className = 'tl-status' + (cls ? ' ' + cls : '');
    }

    // ── Track list (left column) ────────────────────────────────────────────
    function rebuildTrackList() {
        if (!draft || !draft.tracks.length) {
            trackListEl.innerHTML = '<div class="tl-empty">No tracks. Load a preset or add a track.</div>';
            return;
        }
        var html = '';
        for (var i = 0; i < draft.tracks.length; i++) {
            var tr = draft.tracks[i];
            var sel = (tr.path === selectedTrack);
            html += '<button type="button" class="tl-track-item' + (sel ? ' is-sel' : '') +
                '" data-path="' + esc(tr.path) + '">' +
                '<span class="tl-track-path">' + esc(tr.path) + '</span>' +
                '<span class="tl-track-keycount">' + tr.keys.length + ' key' + (tr.keys.length === 1 ? '' : 's') + '</span>' +
                '</button>';
        }
        trackListEl.innerHTML = html;
    }

    trackListEl.addEventListener('click', function(e) {
        var btn = e.target.closest('[data-path]');
        if (!btn) return;
        selectedTrack = btn.getAttribute('data-path');
        selectedKeyT = NaN;
        clearMultiSelect();
        rebuildTrackList();
        rebuildLanes();
        updateInspector();
    });

    // ── Ruler (time header in dopesheet) ────────────────────────────────────
    function updateRuler() {
        var d = getDuration();
        var tickCount = d >= 30 ? 12 : (d >= 10 ? 10 : 8);
        var html = '';
        for (var i = 0; i <= tickCount; i++) {
            var pct = (i / tickCount) * 100;
            var t = (i / tickCount) * d;
            html += '<span class="tl-ruler-tick" style="left:' + pct.toFixed(2) + '%">' +
                '<span class="tl-ruler-label">' + t.toFixed(1) + 's</span></span>';
        }
        // playhead
        var phPct = clamp(currentTime() / Math.max(d, 0.001) * 100, 0, 100);
        html += '<span class="tl-ruler-playhead" id="tl-ruler-playhead" style="left:' + phPct.toFixed(2) + '%"></span>';
        rulerEl.innerHTML = html;
    }

    // Click ruler to seek
    rulerEl.addEventListener('click', function(e) {
        var rect = rulerEl.getBoundingClientRect();
        if (!rect.width) return;
        var ratio = clamp((e.clientX - rect.left) / rect.width, 0, 1);
        var t = ratio * getDuration();
        if (typeof seekPresentation === 'function') seekPresentation(t);
        updateTimeInputs(); updateScrubber(); updatePlayheads();
    });

    // ── Multi-select helpers ────────────────────────────────────────────────
    function isKeyMultiSelected(path, t) {
        for (var i = 0; i < selectedKeys.length; i++) {
            if (selectedKeys[i].path === path && Math.abs(selectedKeys[i].t - t) <= 1e-4) return true;
        }
        return false;
    }
    function addToMultiSelect(path, t) {
        if (!isKeyMultiSelected(path, t)) selectedKeys.push({ path: path, t: t });
    }
    function removeFromMultiSelect(path, t) {
        selectedKeys = selectedKeys.filter(function(s) {
            return !(s.path === path && Math.abs(s.t - t) <= 1e-4);
        });
    }
    function clearMultiSelect() { selectedKeys = []; }
    function selectionCount() { return selectedKeys.length; }
    function selectAllKeysAtTime(t) {
        clearMultiSelect();
        if (!draft) return;
        for (var i = 0; i < draft.tracks.length; i++) {
            var tr = draft.tracks[i];
            for (var j = 0; j < tr.keys.length; j++) {
                if (Math.abs(tr.keys[j].t - t) <= 1e-4) addToMultiSelect(tr.path, tr.keys[j].t);
            }
        }
        rebuildTrackList();
        rebuildLanes();
        if (selectedKeys.length === 1) {
            selectedTrack = selectedKeys[0].path;
            selectedKeyT = selectedKeys[0].t;
            var tk = getTrackByPath(selectedTrack);
            if (tk) fillInspector(tk, getKeyAt(tk, selectedKeyT));
        } else {
            inspSummary.textContent = selectedKeys.length + ' keyframes selected at t=' + t.toFixed(2) + '.';
        }
    }

    // ── Lanes (center dopesheet) ────────────────────────────────────────────
    function rebuildLanes() {
        if (!draft || !draft.tracks.length) {
            lanesEl.innerHTML = '<div class="tl-empty">No tracks.</div>';
            return;
        }
        var d = getDuration();
        var phPct = clamp(currentTime() / Math.max(d, 0.001) * 100, 0, 100);
        var html = '';
        for (var i = 0; i < draft.tracks.length; i++) {
            var tr = draft.tracks[i];
            var isSel = (tr.path === selectedTrack);
            html += '<div class="tl-lane' + (isSel ? ' is-sel' : '') + '" data-path="' + esc(tr.path) + '">';
            html += '<span class="tl-lane-playhead" style="left:' + phPct.toFixed(2) + '%"></span>';
            for (var k = 0; k < tr.keys.length; k++) {
                var key = tr.keys[k];
                var kPct = clamp(key.t / Math.max(d, 0.001) * 100, 0, 100);
                var isSelKey = isKeyMultiSelected(tr.path, key.t);
                var tip = key.t.toFixed(2) + 's | ' + normalizeEase(key.ease) + ' | ' + formatValue(key.v);
                html += '<button type="button" class="tl-diamond' + (isSelKey ? ' is-sel' : '') +
                    '" style="left:' + kPct.toFixed(3) + '%" title="' + esc(tip) +
                    '" data-path="' + esc(tr.path) + '" data-ki="' + k + '"></button>';
            }
            html += '</div>';
        }
        lanesEl.innerHTML = html;
    }

    lanesEl.addEventListener('dblclick', function(e) {
        var diamond = e.target.closest('.tl-diamond');
        if (diamond) {
            var path = diamond.getAttribute('data-path');
            var ki = parseInt(diamond.getAttribute('data-ki'), 10);
            var track = getTrackByPath(path);
            if (track && track.keys[ki]) {
                e.preventDefault();
                selectAllKeysAtTime(track.keys[ki].t);
            }
        }
    });

    var dragJustCommitted = false;

    lanesEl.addEventListener('click', function(e) {
        // Ignore the synthetic click that fires after a drag commit
        if (dragJustCommitted) { dragJustCommitted = false; return; }
        var diamond = e.target.closest('.tl-diamond');
        if (diamond) {
            var path = diamond.getAttribute('data-path');
            var ki = parseInt(diamond.getAttribute('data-ki'), 10);
            var track = getTrackByPath(path);
            if (track && track.keys[ki]) {
                var t = track.keys[ki].t;
                var isCtrl = e.ctrlKey || e.metaKey;
                if (isCtrl) {
                    // Toggle this key in multi-selection
                    if (isKeyMultiSelected(path, t)) {
                        removeFromMultiSelect(path, t);
                    } else {
                        addToMultiSelect(path, t);
                    }
                } else {
                    // Replace selection with just this key
                    clearMultiSelect();
                    addToMultiSelect(path, t);
                }
                selectedTrack = path;
                selectedKeyT = t;
                rebuildTrackList();
                rebuildLanes();
                if (selectedKeys.length === 1) {
                    fillInspector(track, ki);
                } else {
                    inspSummary.textContent = selectedKeys.length + ' keyframes selected.';
                }
            }
            return;
        }
        // Click on lane background = select track + seek to clicked time
        var lane = e.target.closest('.tl-lane');
        if (lane) {
            selectedTrack = lane.getAttribute('data-path');
            selectedKeyT = NaN;
            clearMultiSelect();
            // Seek to the clicked time position
            var rect = lane.getBoundingClientRect();
            if (rect.width > 0) {
                var ratio = clamp((e.clientX - rect.left) / rect.width, 0, 1);
                var t2 = ratio * getDuration();
                if (typeof seekPresentation === 'function') seekPresentation(t2);
                updateTimeInputs(); updateScrubber(); updatePlayheads();
            }
            rebuildTrackList();
            rebuildLanes();
            updateInspector();
        }
    });

    // ── Keyframe drag (move diamonds by dragging) ────────────────────────────
    lanesEl.addEventListener('pointerdown', function(e) {
        if (e.button !== 0) return;
        var diamond = e.target.closest('.tl-diamond');
        if (!diamond) return;
        var path = diamond.getAttribute('data-path');
        var ki   = parseInt(diamond.getAttribute('data-ki'), 10);
        var track = getTrackByPath(path);
        if (!track || !track.keys[ki]) return;

        var clickedT = track.keys[ki].t;
        // If the clicked diamond isn't in the current multi-select, replace selection
        if (!isKeyMultiSelected(path, clickedT)) {
            clearMultiSelect();
            addToMultiSelect(path, clickedT);
            selectedTrack = path;
            selectedKeyT  = clickedT;
            rebuildTrackList();
            rebuildLanes();
            if (selectedKeys.length === 1) fillInspector(track, ki);
        }

        // Collect all selected keys with their original positions
        var rect = lanesEl.getBoundingClientRect();
        var dur  = getDuration();
        var dragKeys = [];
        for (var s = 0; s < selectedKeys.length; s++) {
            var sk = selectedKeys[s];
            var tr = getTrackByPath(sk.path);
            if (!tr) continue;
            var trKi = getKeyAt(tr, sk.t);
            if (trKi < 0) continue;
            dragKeys.push({ path: sk.path, origT: sk.t, ki: trKi });
        }
        keyDragState = {
            pointerId : e.pointerId,
            startX    : e.clientX,
            rect      : rect,
            duration  : dur,
            keys      : dragKeys,
            dragging  : false
        };
        lanesEl.setPointerCapture(e.pointerId);
        e.preventDefault();
    });

    lanesEl.addEventListener('pointermove', function(e) {
        if (!keyDragState) return;
        var dx = e.clientX - keyDragState.startX;
        if (!keyDragState.dragging && Math.abs(dx) > 4) {
            keyDragState.dragging = true;
            lanesEl.classList.add('tl-lanes--dragging');
        }
        if (!keyDragState.dragging) return;
        var dt = (dx / keyDragState.rect.width) * keyDragState.duration;
        var dur = keyDragState.duration;
        for (var i = 0; i < keyDragState.keys.length; i++) {
            var dk = keyDragState.keys[i];
            var newT = clamp(dk.origT + dt, 0, dur);
            var pct  = (newT / Math.max(dur, 0.001) * 100).toFixed(3) + '%';
            var escaped = dk.path.replace(/\\/g, '\\\\').replace(/"/g, '\\"');
            var el = lanesEl.querySelector('.tl-diamond[data-path="' + escaped + '"][data-ki="' + dk.ki + '"]');
            if (el) el.style.left = pct;
        }
        e.preventDefault();
    });

    lanesEl.addEventListener('pointerup', function(e) {
        if (!keyDragState || keyDragState.pointerId !== e.pointerId) return;
        var wasDragging = keyDragState.dragging;
        var drKeys = keyDragState.keys;
        var dx = e.clientX - keyDragState.startX;
        var dt = (dx / keyDragState.rect.width) * keyDragState.duration;
        var dur = keyDragState.duration;
        lanesEl.classList.remove('tl-lanes--dragging');
        keyDragState = null;

        if (wasDragging && draft) {
            pushUndo();
            clearMultiSelect();
            var affectedPaths = {};
            for (var i = 0; i < drKeys.length; i++) {
                var dk = drKeys[i];
                var tr  = getTrackByPath(dk.path);
                if (!tr || !tr.keys[dk.ki]) continue;
                var newT = clamp(dk.origT + dt, 0, dur);
                tr.keys[dk.ki].t = newT;
                addToMultiSelect(dk.path, newT);
                affectedPaths[dk.path] = 1;
            }
            for (var ap in affectedPaths) {
                var tr2 = getTrackByPath(ap);
                if (tr2) tr2.keys.sort(function(a, b) { return a.t - b.t; });
            }
            applyDraft();
            rebuildAll();
            dragJustCommitted = true;
            var count = drKeys.length;
            setStatus(count + ' key' + (count > 1 ? 's' : '') + ' moved.', '');
        }
    });

    lanesEl.addEventListener('pointercancel', function(e) {
        if (!keyDragState) return;
        lanesEl.classList.remove('tl-lanes--dragging');
        keyDragState = null;
        rebuildLanes(); // restore visual positions
    });

    // ── Playhead live-update ────────────────────────────────────────────────
    function updatePlayheads() {
        var d = getDuration();
        var phPct = clamp(currentTime() / Math.max(d, 0.001) * 100, 0, 100).toFixed(2) + '%';
        var rulerPH = panel.querySelector('#tl-ruler-playhead');
        if (rulerPH) rulerPH.style.left = phPct;
        var lanePHs = panel.querySelectorAll('.tl-lane-playhead');
        for (var i = 0; i < lanePHs.length; i++) lanePHs[i].style.left = phPct;
    }

    // ── Inspector (right column) ────────────────────────────────────────────
    function updateInspector() {
        // Multi-selection summary
        if (selectedKeys.length > 1) {
            var pathSet = {};
            for (var s = 0; s < selectedKeys.length; s++) pathSet[selectedKeys[s].path] = 1;
            var pathCount = Object.keys(pathSet).length;
            inspSummary.textContent = selectedKeys.length + ' keyframes selected (' + pathCount + ' track' + (pathCount > 1 ? 's' : '') + ')';
            inspPath.value = '';
            inspTime.value = '';
            inspValue.value = '';
            return;
        }
        var track = getTrackByPath(selectedTrack);
        if (!track) {
            inspSummary.textContent = selectedKeys.length === 0 ? 'No keyframe selected.' : 'No track.';
            inspPath.value = selectedTrack || '';
            inspTime.value = currentTime().toFixed(2);
            inspValue.value = '';
            return;
        }
        if (selectedKeys.length === 1) {
            var sk = selectedKeys[0];
            var sTrack = getTrackByPath(sk.path);
            if (sTrack) {
                var ki = getKeyAt(sTrack, sk.t);
                if (ki >= 0) return fillInspector(sTrack, ki);
            }
        }
        if (isFinite(selectedKeyT)) {
            var ki2 = getKeyAt(track, selectedKeyT);
            if (ki2 >= 0) return fillInspector(track, ki2);
        }
        inspSummary.textContent = 'Track: ' + track.path + ' (' + track.keys.length + ' keys)';
        inspPath.value = track.path;
        inspTime.value = currentTime().toFixed(2);
        inspValue.value = '';
    }

    function fillInspector(track, ki) {
        var key = track.keys[ki];
        if (!key) return;
        inspPath.value = track.path;
        inspTime.value = key.t.toFixed(2);
        inspEase.value = normalizeEase(key.ease);
        inspValue.value = formatValue(key.v);

        var prev = ki > 0 ? track.keys[ki - 1] : null;
        var deltaText = '';
        if (prev && typeof prev.v === 'number' && typeof key.v === 'number') {
            var d = key.v - prev.v;
            deltaText = ' | \u0394 ' + (d >= 0 ? '+' : '') + d.toFixed(4);
        }
        inspSummary.textContent = track.path + ' @ ' + key.t.toFixed(2) + 's' + deltaText;
    }

    // Inspector actions
    inspUseTime.addEventListener('click', function() {
        inspTime.value = currentTime().toFixed(2);
    });
    inspCapture.addEventListener('click', function() {
        var path = (inspPath.value || '').trim();
        if (!path) return;
        if (typeof getPresentationPathValue === 'function') {
            var v = getPresentationPathValue(path);
            if (typeof v !== 'undefined') inspValue.value = formatValue(v);
        }
    });
    inspSet.addEventListener('click', function() {
        if (!draft) return;
        var path = (inspPath.value || '').trim();
        if (!path) { setStatus('Path is required.', 'tl-status--warn'); return; }
        var t = parseTime(inspTime.value, currentTime());
        var ease = normalizeEase(inspEase.value);
        var value = parseValue(inspValue.value);

        pushUndo();
        var track = getTrackByPath(path);
        if (!track) {
            track = { path: path, compile: false, keys: [] };
            draft.tracks.push(track);
            draft.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });
        }
        var existing = getKeyAt(track, t);
        if (existing >= 0) {
            track.keys[existing] = { t: t, v: value, ease: ease };
        } else {
            track.keys.push({ t: t, v: value, ease: ease });
            track.keys.sort(function(a, b) { return a.t - b.t; });
        }
        draft.duration = Math.max(draft.duration, t);
        selectedTrack = path;
        selectedKeyT = t;
        clearMultiSelect();
        addToMultiSelect(path, t);
        normalizeQuatSigns();
        applyDraft();
        rebuildAll();
        setStatus('Key set: ' + path + ' @ ' + t.toFixed(2) + 's', '');
    });
    inspDel.addEventListener('click', function() {
        if (!draft) return;
        // If we have a multi-selection, delete all selected keys
        if (selectedKeys.length > 0) {
            pushUndo();
            var count = 0;
            for (var s = 0; s < selectedKeys.length; s++) {
                var sk = selectedKeys[s];
                var tr = getTrackByPath(sk.path);
                if (!tr) continue;
                var ki = getKeyAt(tr, sk.t);
                if (ki >= 0) { tr.keys.splice(ki, 1); count++; }
            }
            // Remove empty tracks
            draft.tracks = draft.tracks.filter(function(tr) { return tr.keys.length > 0; });
            selectedKeyT = NaN;
            clearMultiSelect();
            applyDraft();
            rebuildAll();
            if (count) setStatus(count + ' key' + (count > 1 ? 's' : '') + ' deleted.', '');
            else setStatus('No keys deleted.', 'tl-status--warn');
            return;
        }
        // Fallback: delete from inspector fields
        var path = (inspPath.value || '').trim();
        var t = parseTime(inspTime.value, NaN);
        var track = getTrackByPath(path);
        if (!track) { setStatus('No track found for path.', 'tl-status--warn'); return; }
        var ki2 = getKeyAt(track, t);
        if (ki2 < 0) { setStatus('No key at that time. Select a keyframe first.', 'tl-status--warn'); return; }
        pushUndo();
        track.keys.splice(ki2, 1);
        if (!track.keys.length) {
            draft.tracks = draft.tracks.filter(function(tr2) { return tr2.path !== path; });
        }
        selectedKeyT = NaN;
        clearMultiSelect();
        applyDraft();
        rebuildAll();
        setStatus('Key deleted.', '');
    });

    // ── Add track ───────────────────────────────────────────────────────────
    addTrackBtn.addEventListener('click', function() {
        var path = prompt('Track path (e.g. observer.distance):');
        if (!path || !path.trim()) return;
        path = path.trim();
        if (!draft) {
            draft = normalizeTL({ name: 'New', duration: 12, tracks: [], events: [] });
        }
        if (getTrackByPath(path)) { setStatus('Track already exists.', 'tl-status--warn'); return; }
        pushUndo();
        var val = '';
        if (typeof getPresentationPathValue === 'function') {
            var v = getPresentationPathValue(path);
            if (typeof v !== 'undefined') val = v;
        }
        draft.tracks.push({ path: path, compile: false, keys: [{ t: currentTime(), v: val, ease: 'linear' }] });
        draft.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });
        selectedTrack = path;
        selectedKeyT = currentTime();
        applyDraft();
        rebuildAll();
        setStatus('Track added: ' + path, '');
    });

    // ── Auto Keyframe ───────────────────────────────────────────────────────
    function captureSnapshot() {
        var out = {};
        function add(p, v) {
            if (typeof v === 'number' && isFinite(v)) { out[p] = v; return; }
            if (typeof v === 'boolean' || typeof v === 'string') out[p] = v;
        }
        function walk(prefix, node, depth) {
            if (!node || typeof node !== 'object' || depth > 8 || Array.isArray(node)) return;
            var keys = Object.keys(node);
            for (var i = 0; i < keys.length; i++) {
                var k = keys[i], val = node[k];
                if (typeof val === 'function') continue;
                var path = prefix ? prefix + '.' + k : k;
                if (val && typeof val === 'object') walk(path, val, depth + 1);
                else add(path, val);
            }
        }
        if (typeof shader !== 'undefined' && shader && shader.parameters) walk('', shader.parameters, 0);
        if (typeof cameraPan !== 'undefined' && cameraPan) { add('cameraPan.x', cameraPan.x); add('cameraPan.y', cameraPan.y); }
        if (typeof camera !== 'undefined' && camera) {
            if (camera.position) { add('camera.position.x', camera.position.x); add('camera.position.y', camera.position.y); add('camera.position.z', camera.position.z); }
            if (camera.quaternion) { add('camera.quaternion.x', camera.quaternion.x); add('camera.quaternion.y', camera.quaternion.y); add('camera.quaternion.z', camera.quaternion.z); add('camera.quaternion.w', camera.quaternion.w); }
        }
        return out;
    }
    autoKeyBtn.addEventListener('click', function() {
        if (!draft) { setStatus('Load a timeline first.', 'tl-status--warn'); return; }
        var t = currentTime();
        var snap = captureSnapshot();
        var paths = Object.keys(snap);
        if (!paths.length) { setStatus('Nothing to capture.', 'tl-status--warn'); return; }

        if (!autoKeySnapshot) {
            autoKeySnapshot = { time: t, values: clonePlain(snap) };
            setStatus('Baseline captured at ' + t.toFixed(2) + 's. Change controls, press AUTO KEY again.', 'tl-status--info');
            return;
        }

        pushUndo();
        var changed = 0;
        for (var i = 0; i < paths.length; i++) {
            var p = paths[i];
            if (!autoKeySnapshot.values.hasOwnProperty(p)) continue;
            if (areValuesEqual(autoKeySnapshot.values[p], snap[p])) continue;
            // upsert both baseline and current
            upsertKey(p, autoKeySnapshot.time, autoKeySnapshot.values[p], 'linear');
            upsertKey(p, t, snap[p], 'linear');
            if (!selectedTrack) selectedTrack = p;
            changed++;
        }
        draft.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });
        draft.duration = Math.max(draft.duration, t, autoKeySnapshot.time);
        autoKeySnapshot = { time: t, values: clonePlain(snap) };

        if (changed) {
            selectedKeyT = t;
            normalizeQuatSigns();
            applyDraft();
            rebuildAll();
            setStatus('Auto-keyed ' + changed + ' changed value' + (changed > 1 ? 's' : '') + ' at ' + t.toFixed(2) + 's.', '');
        } else {
            setStatus('No changes detected. Baseline moved to ' + t.toFixed(2) + 's.', 'tl-status--warn');
        }
    });

    function upsertKey(path, t, value, ease) {
        if (!draft) return;
        var track = getTrackByPath(path);
        if (!track) {
            track = { path: path, compile: false, keys: [] };
            draft.tracks.push(track);
        }
        var idx = getKeyAt(track, t);
        if (idx >= 0) track.keys[idx] = { t: t, v: value, ease: normalizeEase(ease) };
        else { track.keys.push({ t: t, v: value, ease: normalizeEase(ease) }); track.keys.sort(function(a,b){return a.t-b.t;}); }
    }

    // ── Quaternion sign normalisation ───────────────────────────────────────
    // Quaternions q and -q represent the same rotation. When the four components
    // (x,y,z,w) are stored as independent numeric tracks and interpolated
    // component-by-component, adjacent keys in opposite hemispheres (dot < 0)
    // will interpolate the long way around.  This function walks every group of
    // tracks whose paths share a common prefix and end in .x/.y/.z/.w, and
    // negates any key whose quaternion dot-product with the previous normalised
    // key is negative, ensuring the shortest arc is always taken.
    function normalizeQuatSigns() {
        if (!draft) return;
        // Collect unique prefixes for xyzw groups
        var groups = {};
        for (var i = 0; i < draft.tracks.length; i++) {
            var m = draft.tracks[i].path.match(/^(.+)\.(x|y|z|w)$/);
            if (m) groups[m[1]] = true;
        }
        var prefixes = Object.keys(groups);
        for (var pi = 0; pi < prefixes.length; pi++) {
            var prefix = prefixes[pi];
            var tx = getTrackByPath(prefix + '.x');
            var ty = getTrackByPath(prefix + '.y');
            var tz = getTrackByPath(prefix + '.z');
            var tw = getTrackByPath(prefix + '.w');
            if (!tx || !ty || !tz || !tw) continue;
            // Iterate over x-track key times in order; all four components must
            // share the same set of times (AutoKey guarantees this).
            for (var k = 1; k < tx.keys.length; k++) {
                var t0 = tx.keys[k - 1].t;
                var t1 = tx.keys[k].t;
                // Previous quat (already normalised in earlier iterations)
                var px = quatCompAt(tx, t0), py = quatCompAt(ty, t0);
                var pz = quatCompAt(tz, t0), pw = quatCompAt(tw, t0);
                // Current quat
                var cx = quatCompAt(tx, t1), cy = quatCompAt(ty, t1);
                var cz = quatCompAt(tz, t1), cw = quatCompAt(tw, t1);
                if (px * cx + py * cy + pz * cz + pw * cw < 0) {
                    setQuatCompAt(tx, t1, -cx); setQuatCompAt(ty, t1, -cy);
                    setQuatCompAt(tz, t1, -cz); setQuatCompAt(tw, t1, -cw);
                }
            }
        }
    }
    function quatCompAt(track, t) {
        for (var i = 0; i < track.keys.length; i++)
            if (Math.abs(track.keys[i].t - t) <= 1e-4) return +track.keys[i].v;
        return 0;
    }
    function setQuatCompAt(track, t, v) {
        for (var i = 0; i < track.keys.length; i++)
            if (Math.abs(track.keys[i].t - t) <= 1e-4) { track.keys[i].v = v; return; }
    }

    importBtn.addEventListener('click', function() {
        var input = document.createElement('input');
        input.type = 'file';
        input.accept = '.json,application/json';
        input.style.display = 'none';
        input.addEventListener('change', function() {
            var file = input.files && input.files[0];
            if (!file) return;
            var reader = new FileReader();
            reader.onload = function() {
                try {
                    var obj = JSON.parse(reader.result);
                    if (typeof setPresentationTimeline === 'function') {
                        if (setPresentationTimeline(obj)) {
                            syncFromRuntime();
                            setStatus('Imported: ' + file.name, '');
                        } else {
                            setStatus('Invalid timeline data in ' + file.name, 'tl-status--error');
                        }
                    }
                } catch (err) {
                    setStatus('JSON parse error: ' + err.message, 'tl-status--error');
                }
            };
            reader.onerror = function() {
                setStatus('Failed to read file.', 'tl-status--error');
            };
            reader.readAsText(file);
        });
        document.body.appendChild(input);
        input.click();
        input.remove();
    });

    exportBtn.addEventListener('click', function() {
        var tl = typeof getPresentationTimeline === 'function' ? getPresentationTimeline() : null;
        if (!tl && draft) tl = clonePlain(draft);
        if (!tl) { setStatus('No timeline to export.', 'tl-status--warn'); return; }
        var json = JSON.stringify(tl, null, 2);
        var blob = new Blob([json], { type: 'application/json' });
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = (tl.name ? tl.name.replace(/[^a-zA-Z0-9_-]/g, '_') : 'timeline') + '.json';
        document.body.appendChild(a);
        a.click();
        a.remove();
        URL.revokeObjectURL(url);
        setStatus('Exported: ' + a.download, '');
    });

    // ── Listen for external timeline loads ──────────────────────────────────
    window.addEventListener('presentation:timeline-panel-sync', function() {
        if (panelOpen) syncFromRuntime();
    });

    // ══ Motion Functions ═════════════════════════════════════════════════════
    // Each entry: params[] with { id, label, type, min, step, def, defaultFn, options }
    var MOTION_TYPES = {
        orbit: {
            params: [
                { id: 'start',    label: 'Start time (s)',     type: 'number', min: 0,    step: 0.1,  defaultFn: function() { return currentTime().toFixed(2); } },
                { id: 'duration', label: 'Duration (s)',       type: 'number', min: 0.1,  step: 1,    def: 10 },
                { id: 'orbits',   label: 'Number of orbits',   type: 'number', min: 0.1,  step: 0.5,  def: 1 },
                { id: 'dir',      label: 'Direction',          type: 'select', options: [['ccw','Counter-clockwise (CCW)'],['cw','Clockwise (CW)']] }
            ]
        },
        zoom: {
            params: [
                { id: 'start',    label: 'Start time (s)',     type: 'number', min: 0,    step: 0.1,  defaultFn: function() { return currentTime().toFixed(2); } },
                { id: 'duration', label: 'Duration (s)',       type: 'number', min: 0.1,  step: 1,    def: 8 },
                { id: 'from',     label: 'From distance',      type: 'number', min: 1,    step: 0.5,  defaultFn: function() { var v = typeof getPresentationPathValue === 'function' ? getPresentationPathValue('observer.distance') : undefined; return (typeof v === 'number' ? v : 15).toFixed(2); } },
                { id: 'to',       label: 'To distance',        type: 'number', min: 1,    step: 0.5,  def: 7 },
                { id: 'ease',     label: 'Ease',               type: 'select', options: [['smooth','smooth'],['smoother','smoother'],['linear','linear']] }
            ]
        },
        exposure: {
            params: [
                { id: 'start',    label: 'Start time (s)',     type: 'number', min: 0,    step: 0.1,  defaultFn: function() { return currentTime().toFixed(2); } },
                { id: 'duration', label: 'Duration (s)',       type: 'number', min: 0.1,  step: 1,    def: 6 },
                { id: 'from',     label: 'From exposure',      type: 'number', min: 0,    step: 0.05, defaultFn: function() { var v = typeof getPresentationPathValue === 'function' ? getPresentationPathValue('look.exposure') : undefined; return (typeof v === 'number' ? v : 1.0).toFixed(3); } },
                { id: 'to',       label: 'To exposure',        type: 'number', min: 0,    step: 0.05, def: 1.5 },
                { id: 'ease',     label: 'Ease',               type: 'select', options: [['smooth','smooth'],['smoother','smoother'],['linear','linear']] }
            ]
        },
        inclination: {
            params: [
                { id: 'start',    label: 'Start time (s)',     type: 'number', min: 0,    step: 0.1,  defaultFn: function() { return currentTime().toFixed(2); } },
                { id: 'duration', label: 'Duration (s)',       type: 'number', min: 0.1,  step: 1,    def: 8 },
                { id: 'from',     label: 'From angle (\u00b0)',       type: 'number', min: -90,  step: 5,    defaultFn: function() { var v = typeof getPresentationPathValue === 'function' ? getPresentationPathValue('observer.orbital_inclination') : undefined; return (typeof v === 'number' ? v : 0).toFixed(1); } },
                { id: 'to',       label: 'To angle (\u00b0)',         type: 'number', min: -90,  step: 5,    def: 30 },
                { id: 'ease',     label: 'Ease',               type: 'select', options: [['smooth','smooth'],['smoother','smoother'],['linear','linear']] }
            ]
        }
    };

    function renderMotionParams() {
        var type = motionTypeEl.value;
        var def = MOTION_TYPES[type];
        if (!def) { motionParamsEl.innerHTML = ''; return; }
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

    // Three.js quaternion: camera at (px, py, pz) looking toward world origin.
    // The camera's +Z axis (which points AWAY from the look target) = normalize(P).
    function lookAtOriginQuat(px, py, pz) {
        var r = Math.sqrt(px*px + py*py + pz*pz);
        if (r < 1e-8) return { x: 0, y: 0, z: 0, w: 1 };
        var zx = px/r, zy = py/r, zz = pz/r; // cam +Z in world
        // world up — switch to world Z when camera points straight up/down
        var upx = 0, upy = 1, upz = 0;
        if (Math.abs(zy) > 0.999) { upx = 0; upy = 0; upz = 1; }
        // cam +X = right = normalize(worldUp × camZ)
        var rx = upy*zz - upz*zy, ry = upz*zx - upx*zz, rz = upx*zy - upy*zx;
        var rlen = Math.sqrt(rx*rx + ry*ry + rz*rz);
        if (rlen < 1e-8) return { x: 0, y: 0, z: 0, w: 1 };
        rx /= rlen; ry /= rlen; rz /= rlen;
        // cam +Y = camZ × camX
        var ux = zy*rz - zz*ry, uy = zz*rx - zx*rz, uz = zx*ry - zy*rx;
        // Rotation matrix columns: [camX, camY, camZ]
        var m00 = rx,  m01 = ux,  m02 = zx;
        var m10 = ry,  m11 = uy,  m12 = zy;
        var m20 = rz,  m21 = uz,  m22 = zz;
        var trace = m00 + m11 + m22;
        var qx, qy, qz, qw;
        if (trace > 0) {
            var s = 0.5 / Math.sqrt(trace + 1.0);
            qw = 0.25/s; qx = (m21-m12)*s; qy = (m02-m20)*s; qz = (m10-m01)*s;
        } else if (m00 > m11 && m00 > m22) {
            var s2 = 2.0*Math.sqrt(1.0+m00-m11-m22);
            qw = (m21-m12)/s2; qx = 0.25*s2; qy = (m01+m10)/s2; qz = (m02+m20)/s2;
        } else if (m11 > m22) {
            var s3 = 2.0*Math.sqrt(1.0+m11-m00-m22);
            qw = (m02-m20)/s3; qx = (m01+m10)/s3; qy = 0.25*s3; qz = (m12+m21)/s3;
        } else {
            var s4 = 2.0*Math.sqrt(1.0+m22-m00-m11);
            qw = (m10-m01)/s4; qx = (m02+m20)/s4; qy = (m12+m21)/s4; qz = 0.25*s4;
        }
        return { x: qx, y: qy, z: qz, w: qw };
    }

    // Linear-interpolate a draft track at a given time (used to seed orbit start)
    function sampleTrackDraft(track, t) {
        var keys = track.keys;
        if (!keys || !keys.length) return 0;
        if (t <= keys[0].t) return +keys[0].v;
        if (t >= keys[keys.length-1].t) return +keys[keys.length-1].v;
        for (var i = 0; i < keys.length-1; i++) {
            if (t >= keys[i].t && t <= keys[i+1].t) {
                var dt = keys[i+1].t - keys[i].t;
                if (dt < 1e-8) return +keys[i].v;
                return +keys[i].v + (+keys[i+1].v - +keys[i].v) * ((t - keys[i].t) / dt);
            }
        }
        return +keys[keys.length-1].v;
    }

    function applyMotionFn() {
        if (!draft) { setStatus('Load or create a timeline first.', 'tl-status--warn'); return; }
        var type = motionTypeEl.value;
        var start    = getMotionParam('start');    if (!isFinite(start))    start    = 0;
        var duration = getMotionParam('duration'); if (!isFinite(duration) || duration < 0.01) duration = 5;
        var ease     = getMotionParam('ease')    || 'smooth';

        pushUndo();

        if (type === 'orbit') {
            var numOrbits = getMotionParam('orbits'); if (!isFinite(numOrbits) || numOrbits < 0.01) numOrbits = 1;
            var dirSign   = getMotionParam('dir') === 'cw' ? -1 : 1;

            // Seed camera position — prefer draft track value at start time, then runtime
            var cpx = 0, cpy = 0, cpz = 11;
            if (typeof getPresentationPathValue === 'function') {
                var vx = getPresentationPathValue('camera.position.x');
                var vy = getPresentationPathValue('camera.position.y');
                var vz = getPresentationPathValue('camera.position.z');
                if (typeof vx === 'number') cpx = vx;
                if (typeof vy === 'number') cpy = vy;
                if (typeof vz === 'number') cpz = vz;
            }
            var txd = getTrackByPath('camera.position.x');
            var tyd = getTrackByPath('camera.position.y');
            var tzd = getTrackByPath('camera.position.z');
            if (txd && txd.keys.length) cpx = sampleTrackDraft(txd, start);
            if (tyd && tyd.keys.length) cpy = sampleTrackDraft(tyd, start);
            if (tzd && tzd.keys.length) cpz = sampleTrackDraft(tzd, start);

            var dist = Math.sqrt(cpx*cpx + cpy*cpy + cpz*cpz);
            if (dist < 1e-6) { dist = 11; cpz = 11; cpy = 0; cpx = 0; }

            // Spherical coords (Y-up): keep elevation constant, sweep azimuth
            var theta    = Math.asin(Math.max(-1, Math.min(1, cpy / dist))); // elevation
            var phi      = Math.atan2(cpz, cpx);                             // azimuth in XZ-plane
            var cosTheta = Math.cos(theta);
            var totalAngle = dirSign * numOrbits * 2 * Math.PI;
            // Use 32 steps per orbit — all linear so velocity is perfectly constant.
            // Each step spans 11.25 degrees; chord/arc deviation is < 0.05 % at that density.
            var numSteps   = Math.max(8, Math.round(numOrbits * 32));

            for (var si = 0; si <= numSteps; si++) {
                var frac  = si / numSteps;
                var angle = phi + frac * totalAngle;
                var ktime = start + frac * duration;
                var kx = dist * cosTheta * Math.cos(angle);
                var ky = dist * Math.sin(theta);
                var kz = dist * cosTheta * Math.sin(angle);
                var q  = lookAtOriginQuat(kx, ky, kz);
                // Always linear — smooth/smoother ease would decelerate at every
                // intermediate keyframe, creating stops every 11.25 degrees.
                upsertKey('camera.position.x',   ktime, kx,  'linear');
                upsertKey('camera.position.y',   ktime, ky,  'linear');
                upsertKey('camera.position.z',   ktime, kz,  'linear');
                upsertKey('camera.quaternion.x', ktime, q.x, 'linear');
                upsertKey('camera.quaternion.y', ktime, q.y, 'linear');
                upsertKey('camera.quaternion.z', ktime, q.z, 'linear');
                upsertKey('camera.quaternion.w', ktime, q.w, 'linear');
            }
            normalizeQuatSigns();
            selectedTrack = 'camera.position.x';
            setStatus('Orbit: ' + numOrbits + ' orbit(s) over ' + duration + 's starting at t=' + start.toFixed(2) + 's.', '');

        } else if (type === 'zoom') {
            var fromV = getMotionParam('from'); if (!isFinite(fromV)) fromV = 15;
            var toV   = getMotionParam('to');   if (!isFinite(toV))   toV   = 7;
            upsertKey('observer.distance', start,            fromV, 'linear');
            upsertKey('observer.distance', start + duration, toV,   ease);
            selectedTrack = 'observer.distance';
            setStatus('Zoom: distance ' + fromV + ' \u2192 ' + toV + ' over ' + duration + 's.', '');

        } else if (type === 'exposure') {
            var fromV = getMotionParam('from'); if (!isFinite(fromV)) fromV = 1.0;
            var toV   = getMotionParam('to');   if (!isFinite(toV))   toV   = 1.5;
            upsertKey('look.exposure', start,            fromV, 'linear');
            upsertKey('look.exposure', start + duration, toV,   ease);
            selectedTrack = 'look.exposure';
            setStatus('Exposure: ' + fromV + ' \u2192 ' + toV + ' over ' + duration + 's.', '');

        } else if (type === 'inclination') {
            var fromV = getMotionParam('from'); if (!isFinite(fromV)) fromV = 0;
            var toV   = getMotionParam('to');   if (!isFinite(toV))   toV   = 30;
            upsertKey('observer.orbital_inclination', start,            fromV, 'linear');
            upsertKey('observer.orbital_inclination', start + duration, toV,   ease);
            selectedTrack = 'observer.orbital_inclination';
            setStatus('Inclination: ' + fromV + '\u00b0 \u2192 ' + toV + '\u00b0 over ' + duration + 's.', '');
        }

        draft.duration = Math.max(draft.duration, start + duration);
        draft.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });
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

    motionBtn.addEventListener('click', function(e) {
        e.stopPropagation();
        if (motionModal.classList.contains('is-open')) closeMotionModal();
        else openMotionModal();
    });
    motionCloseBtn.addEventListener('click', function(e) {
        e.stopPropagation();
        closeMotionModal();
    });
    motionTypeEl.addEventListener('change', renderMotionParams);
    motionApplyBtn.addEventListener('click', function(e) {
        e.stopPropagation();
        applyMotionFn();
    });
    // Close modal on click outside
    panel.addEventListener('click', function(e) {
        if (motionModal.classList.contains('is-open') &&
            !motionModal.contains(e.target) &&
            e.target !== motionBtn) {
            closeMotionModal();
        }
    }, true);

    // ── Recording modal ─────────────────────────────────────────────────────
    function clampNum(v, lo, hi) { return Math.min(hi, Math.max(lo, v)); }

    function formatRecDuration(s) {
        if (isNaN(s) || s < 0) return '--:--';
        var m = Math.floor(s / 60), sec = Math.floor(s % 60);
        return (m < 10 ? '0' : '') + m + ':' + (sec < 10 ? '0' : '') + sec;
    }

    function setRecStatus(text, cls) {
        recStatusEl.textContent = text;
        recStatusEl.className   = 'tl-rec-status' + (cls ? ' ' + cls : '');
    }

    function populateRecQuality() {
        if (!recQualitySelect) return;
        var cur = recQualitySelect.value;
        recQualitySelect.innerHTML = '';
        var presets = (typeof QUALITY_PRESETS !== 'undefined' && QUALITY_PRESETS)
            ? Object.keys(QUALITY_PRESETS) : [];
        var keys = presets.length ? presets : ['optimal', 'high', 'ultra', 'cinematic', 'medium', 'mobile'];
        for (var i = 0; i < keys.length; i++) {
            var o = document.createElement('option');
            o.value = keys[i];
            o.textContent = keys[i].charAt(0).toUpperCase() + keys[i].slice(1);
            recQualitySelect.appendChild(o);
        }
        if (cur) recQualitySelect.value = cur;
        if (!recQualitySelect.value && recQualitySelect.options.length) {
            recQualitySelect.value = recQualitySelect.options[0].value;
        }
    }

    function populateRecMode() {
        if (!recModeSelect) return;
        var cur = recModeSelect.value;
        recModeSelect.innerHTML = '';
        var state = (typeof getPresentationState === 'function') ? getPresentationState() : {};
        var modes = [
            { v: 'offline',  l: 'Offline (fixed FPS)' },
            { v: 'realtime', l: 'Realtime (screen capture)' }
        ];
        for (var i = 0; i < modes.length; i++) {
            var o = document.createElement('option');
            o.value = modes[i].v;
            o.textContent = modes[i].l;
            recModeSelect.appendChild(o);
        }
        if (cur) recModeSelect.value = cur;
        if (!recModeSelect.value) recModeSelect.value = 'offline';
    }

    function refreshRecResolutionLabel() {
        if (!recResSelect) return;
        var state = (typeof getPresentationState === 'function') ? getPresentationState() : {};
        for (var i = 0; i < recResSelect.options.length; i++) {
            var opt = recResSelect.options[i];
            if (opt.value === 'current') {
                var w = (state.recording_output_width || window.innerWidth);
                var h = (state.recording_output_height || window.innerHeight);
                opt.textContent = 'Current viewport (' + w + '\xd7' + h + ')';
                break;
            }
        }
    }

    function populateRecResolution() {
        if (!recResSelect) return;
        var cur = recResSelect.value;
        recResSelect.innerHTML = '';
        var opts = [
            { v: 'current',   l: 'Current viewport' },
            { v: '1280x720',  l: '1280\xd7720 (HD)' },
            { v: '1920x1080', l: '1920\xd71080 (Full HD)' },
            { v: '2560x1440', l: '2560\xd71440 (2K)' },
            { v: '3840x2160', l: '3840\xd72160 (4K)' }
        ];
        for (var i = 0; i < opts.length; i++) {
            var o = document.createElement('option');
            o.value = opts[i].v;
            o.textContent = opts[i].l;
            recResSelect.appendChild(o);
        }
        if (cur) recResSelect.value = cur;
        if (!recResSelect.value) recResSelect.value = 'current';
        refreshRecResolutionLabel();
    }

    function syncRecModal() {
        if (!recModal || !recModal.classList.contains('is-open')) {
            if (recBtn) {
                var st = (typeof getPresentationState === 'function') ? getPresentationState() : {};
                recBtn.classList.toggle('is-recording', !!st.recording);
            }
            return;
        }
        if (typeof getPresentationState !== 'function') return;
        var s = getPresentationState();

        // Sync loop / annotations checkboxes from state
        if (recLoopCb)        recLoopCb.checked        = !!s.loop;
        if (recAnnotCb)       recAnnotCb.checked       = !!s.annotations_enabled;
        if (recAnnotRecordCb) recAnnotRecordCb.checked = !!s.annotations_in_recording;

        // Refresh resolution label with live viewport size
        refreshRecResolutionLabel();

        // Keep REC button indicator in sync
        if (recBtn) recBtn.classList.toggle('is-recording', !!s.recording);

        if (!s.recording) {
            recStartBtn.disabled = false;
            recStopBtn.disabled  = true;
            setRecStatus('Idle', '');
            return;
        }

        recStartBtn.disabled = true;
        recStopBtn.disabled  = false;

        if (s.recording_mode === 'realtime') {
            if (s.recording_background_throttle_detected) {
                setRecStatus('\u26a0 Background throttle detected \u2014 keep window focused!', 'is-warning');
            } else {
                setRecStatus('Recording\u2026 (realtime)', 'is-recording');
            }
            return;
        }

        // Offline mode
        var phase = s.recording_offline_phase || '';
        if (phase === 'rendering') {
            var pct  = Math.round((s.recording_offline_progress || 0) * 100);
            var done = s.recording_offline_frames_done || 0;
            var total= s.recording_offline_frames_total || 0;
            var fps  = s.recording_offline_render_fps ? (' @ ' + s.recording_offline_render_fps.toFixed(1) + ' fps') : '';
            var eta  = (s.recording_offline_eta_s != null && s.recording_offline_eta_s >= 0)
                ? (' ETA ' + formatRecDuration(s.recording_offline_eta_s)) : '';
            setRecStatus('Rendering ' + pct + '% (' + done + '/' + total + ')' + fps + eta, 'is-recording');
        } else if (phase === 'finalizing') {
            var fp = s.recording_offline_finalizing_progress;
            var sub = s.recording_offline_finalizing_sub || '';
            var label = sub === 'encode' ? 'Encoding' : sub === 'mux' ? 'Muxing' : sub === 'download' ? 'Downloading' : 'Finalizing';
            var pctStr = (fp != null && fp >= 0) ? (' ' + Math.round(fp * 100) + '%') : '\u2026';
            setRecStatus(label + pctStr, 'is-recording');
        } else {
            setRecStatus('Recording\u2026', 'is-recording');
        }
    }

    function openRecModal() {
        populateRecQuality();
        populateRecMode();
        populateRecResolution();
        syncRecModal();
        recModal.classList.add('is-open');
        recBtn.classList.add('is-active');
    }
    function closeRecModal() {
        recModal.classList.remove('is-open');
        recBtn.classList.remove('is-active');
    }

    recBtn.addEventListener('click', function(e) {
        e.stopPropagation();
        if (recModal.classList.contains('is-open')) closeRecModal();
        else openRecModal();
    });
    recCloseBtn.addEventListener('click', function(e) {
        e.stopPropagation();
        closeRecModal();
    });

    recLoopCb.addEventListener('change', function() {
        if (typeof setPresentationLoop === 'function') setPresentationLoop(recLoopCb.checked);
    });
    recAnnotCb.addEventListener('change', function() {
        if (typeof setPresentationAnnotationsEnabled === 'function') setPresentationAnnotationsEnabled(recAnnotCb.checked);
    });
    recAnnotRecordCb.addEventListener('change', function() {
        if (typeof setPresentationAnnotationsIncludedInRecording === 'function') setPresentationAnnotationsIncludedInRecording(recAnnotRecordCb.checked);
    });

    recStartBtn.addEventListener('click', function() {
        if (typeof startPresentationRecording !== 'function') return;
        if (recResetSimCb && recResetSimCb.checked) {
            if (typeof observer !== 'undefined' && observer) observer.time = 0.0;
            if (typeof shader !== 'undefined' && shader) shader.needsUpdate = true;
        }
        var fps     = clampNum(parseFloat(recFpsInput.value) || 60, 24, 120);
        var bitrate = clampNum(parseFloat(recBitrateInput.value) || 20, 4, 80);
        recFpsInput.value     = Math.round(fps);
        recBitrateInput.value = Math.round(bitrate);
        var started = startPresentationRecording({
            fps: fps,
            bitrateMbps: bitrate,
            autoStopOnPresentationEnd: true,
            recordingMode:       recModeSelect    ? recModeSelect.value    : 'offline',
            recordingResolution: recResSelect     ? recResSelect.value     : 'current',
            qualityPreset:       recQualitySelect ? recQualitySelect.value : 'optimal',
            includeAnnotationsInRecording: !!recAnnotRecordCb.checked
        });
        if (!started) {
            var st = (typeof getPresentationState === 'function') ? getPresentationState() : {};
            setRecStatus(st.recording_offline_unavailable_reason || 'Failed to start recording.', 'is-warning');
        }
        syncRecModal();
    });
    recStopBtn.addEventListener('click', function() {
        if (typeof stopPresentationRecording === 'function') stopPresentationRecording();
        syncRecModal();
    });

    // Close REC modal on click outside
    panel.addEventListener('click', function(e) {
        if (recModal.classList.contains('is-open') &&
            !recModal.contains(e.target) &&
            e.target !== recBtn) {
            closeRecModal();
        }
    }, true);

    // ── Keyboard shortcuts ──────────────────────────────────────────────────
    function deleteSelectedKeys() {
        if (!draft || !selectedKeys.length) return;
        pushUndo();
        var count = 0;
        for (var s = 0; s < selectedKeys.length; s++) {
            var sk = selectedKeys[s];
            var tr = getTrackByPath(sk.path);
            if (!tr) continue;
            var ki = getKeyAt(tr, sk.t);
            if (ki >= 0) { tr.keys.splice(ki, 1); count++; }
        }
        draft.tracks = draft.tracks.filter(function(tr) { return tr.keys.length > 0; });
        selectedKeyT = NaN;
        clearMultiSelect();
        applyDraft();
        rebuildAll();
        if (count) setStatus(count + ' key' + (count > 1 ? 's' : '') + ' deleted.', '');
    }

    function selectAllKeysOnTrack() {
        if (!draft || !selectedTrack) return;
        var track = getTrackByPath(selectedTrack);
        if (!track) return;
        clearMultiSelect();
        for (var k = 0; k < track.keys.length; k++) {
            addToMultiSelect(track.path, track.keys[k].t);
        }
        rebuildLanes();
        inspSummary.textContent = selectedKeys.length + ' keyframes selected on ' + selectedTrack;
    }

    function selectAllKeys() {
        if (!draft) return;
        clearMultiSelect();
        for (var i = 0; i < draft.tracks.length; i++) {
            var tr = draft.tracks[i];
            for (var k = 0; k < tr.keys.length; k++) {
                addToMultiSelect(tr.path, tr.keys[k].t);
            }
        }
        rebuildLanes();
        inspSummary.textContent = selectedKeys.length + ' keyframes selected (all).';
    }

    // ── Copy / Paste keyframes ──────────────────────────────────────────────
    function copySelectedKeys() {
        if (!draft || !selectedKeys.length) return;
        var anchorT = Infinity;
        for (var s = 0; s < selectedKeys.length; s++) {
            if (selectedKeys[s].t < anchorT) anchorT = selectedKeys[s].t;
        }
        var entries = [];
        for (var s = 0; s < selectedKeys.length; s++) {
            var sk = selectedKeys[s];
            var tr = getTrackByPath(sk.path);
            if (!tr) continue;
            var ki = getKeyAt(tr, sk.t);
            if (ki < 0) continue;
            var key = tr.keys[ki];
            entries.push({ path: sk.path, relT: sk.t - anchorT, v: clonePlain(key.v), ease: key.ease });
        }
        if (!entries.length) return;
        clipboard = { anchorT: anchorT, entries: entries };
        setStatus(entries.length + ' key' + (entries.length > 1 ? 's' : '') + ' copied.', '');
    }

    function pasteKeys() {
        if (!draft || !clipboard) return;
        var pasteT = currentTime();
        pushUndo();
        clearMultiSelect();
        for (var i = 0; i < clipboard.entries.length; i++) {
            var entry = clipboard.entries[i];
            var t = clamp(pasteT + entry.relT, 0, getDuration());
            var track = getTrackByPath(entry.path);
            if (!track) {
                track = { path: entry.path, compile: false, keys: [] };
                draft.tracks.push(track);
                draft.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });
            }
            var existing = getKeyAt(track, t);
            if (existing >= 0) {
                track.keys[existing] = { t: t, v: clonePlain(entry.v), ease: entry.ease };
            } else {
                track.keys.push({ t: t, v: clonePlain(entry.v), ease: entry.ease });
                track.keys.sort(function(a, b) { return a.t - b.t; });
            }
            draft.duration = Math.max(draft.duration, t);
            addToMultiSelect(entry.path, t);
        }
        applyDraft();
        rebuildAll();
        setStatus(clipboard.entries.length + ' key' + (clipboard.entries.length > 1 ? 's' : '') + ' pasted.', '');
    }

    window.addEventListener('keydown', function(e) {
        if (!panelOpen) return;
        // Don't capture when typing in input/textarea (except for specific shortcuts)
        var tag = (e.target.tagName || '').toLowerCase();
        var inInput = (tag === 'input' || tag === 'textarea' || tag === 'select' || e.target.isContentEditable);

        // Ctrl+Z / Ctrl+Y — always active even in inputs for the panel
        if ((e.ctrlKey || e.metaKey) && !e.shiftKey && e.code === 'KeyZ' && panel.contains(e.target)) {
            e.preventDefault();
            undo();
            return;
        }
        if ((e.ctrlKey || e.metaKey) && (e.code === 'KeyY' || (e.shiftKey && e.code === 'KeyZ')) && panel.contains(e.target)) {
            e.preventDefault();
            redo();
            return;
        }

        if (inInput) return;

        // Escape → close motion modal if open
        if (e.code === 'Escape') {
            if (motionModal.classList.contains('is-open')) { closeMotionModal(); e.preventDefault(); }
            return;
        }

        // Space → play / pause
        if (e.code === 'Space') {
            e.preventDefault();
            var st = typeof getPresentationState === 'function' ? getPresentationState() : null;
            if (st && st.active && !st.paused) {
                if (typeof pausePresentation === 'function') pausePresentation();
            } else {
                if (typeof playPresentation === 'function') playPresentation(false);
            }
            return;
        }

        // Delete / Backspace → delete selected keyframes
        if (e.code === 'Delete' || e.code === 'Backspace') {
            e.preventDefault();
            deleteSelectedKeys();
            return;
        }

        // Ctrl+A → select all keyframes (on current track if one is selected, else all)
        if ((e.ctrlKey || e.metaKey) && e.code === 'KeyA') {
            e.preventDefault();
            if (selectedTrack && !e.shiftKey) selectAllKeysOnTrack();
            else selectAllKeys();
            return;
        }

        // Ctrl+C → copy selected keyframes to clipboard
        if ((e.ctrlKey || e.metaKey) && !e.shiftKey && e.code === 'KeyC') {
            e.preventDefault();
            copySelectedKeys();
            return;
        }

        // Ctrl+V → paste keyframes at current time
        if ((e.ctrlKey || e.metaKey) && !e.shiftKey && e.code === 'KeyV') {
            e.preventDefault();
            pasteKeys();
            return;
        }

        // Home → seek to start
        if (e.code === 'Home') {
            e.preventDefault();
            if (typeof seekPresentation === 'function') seekPresentation(0);
            updateTimeInputs(); updateScrubber(); updatePlayheads();
            return;
        }

        // End → seek to end
        if (e.code === 'End') {
            e.preventDefault();
            if (typeof seekPresentation === 'function') seekPresentation(getDuration());
            updateTimeInputs(); updateScrubber(); updatePlayheads();
            return;
        }

        // Left/Right arrows → nudge time
        if (e.code === 'ArrowLeft' || e.code === 'ArrowRight') {
            var step = e.shiftKey ? 1.0 : 0.1;
            var dir = (e.code === 'ArrowLeft') ? -1 : 1;
            var newT = clamp(currentTime() + dir * step, 0, getDuration());
            if (typeof seekPresentation === 'function') seekPresentation(newT);
            updateTimeInputs(); updateScrubber(); updatePlayheads();
            e.preventDefault();
            return;
        }
    });

    // ── Public API ──────────────────────────────────────────────────────────
    timelinePanelBinding = {
        open: function() { setPanelOpen(true); },
        close: function() { setPanelOpen(false); },
        toggle: toggle,
        isOpen: function() { return panelOpen; },
        sync: syncFromRuntime,
        loadPreset: loadPresetByName,
        syncRecState: syncRecModal
    };
    return timelinePanelBinding;
}
