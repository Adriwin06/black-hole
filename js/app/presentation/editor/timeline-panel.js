// â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
// Timeline Panel â€” bottom-docked dopesheet for animation editing
// â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
// Provides a full-width bottom panel similar to After Effects / Blender dopesheet
// with transport controls, track list, keyframe lanes, ruler, and key inspector.
// Public API exposed via the module export and runtime UI registry.

import { registerBlackHoleUiBinding } from '../../core/runtime/runtime-registry.js';
import {
    PRESENTATION_PRESETS,
    PRESENTATION_PRESET_ORDER,
    registerPresentationPreset,
    getPresentationTimeline,
    setPresentationTimeline,
    seekPresentation,
    getPresentationState,
    playPresentation,
    pausePresentation,
    stopPresentation,
    listPresentationPresets,
    loadPresentationPreset,
    getPresentationAnnotationsState,
    getPresentationParamHudState,
    isParamInHud,
    toggleParamInHud,
    getPresentationPathValue,
    setPresentationAnnotationsEnabled,
    setPresentationAnnotationsIncludedInRecording,
    setPresentationParamHudEnabled,
    setPresentationParamHudIncludedInRecording,
    removeParamFromHud,
    addParamToHud,
    clearParamHud,
    setParamHudLayout,
    capturePresentationScreenshot,
    startPresentationRecording,
    stopPresentationRecording
} from '../runtime/presentation-controller.js';
import {
    clamp,
    esc,
    clonePlain,
    parseTime,
    normalizeEase,
    formatValue,
    parseValue,
    areValuesEqual,
    normalizeTimelineDraft
} from './timeline-utils.js';
import {
    lookAtOriginQuat,
    quatMul,
    quatSlerp,
    sampleTrackDraft
} from './motion-utils.js';

var PRESENTATION_EDITOR_COMMON_PATHS = [
    'dive.currentR',
    'hover.currentR',
    'observerState.time',
    'observer.distance',
    'observer.orbital_inclination',
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

export var timelinePanelBinding = null;

export function buildTimelinePanel() {
    'use strict';

    // â”€â”€ Utility helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    var panel = document.createElement('div');
    panel.id = 'tl-panel';
    panel.className = 'tl-panel tl-panel--collapsed';
    panel.innerHTML =
        '<div id="tl-resize-handle" class="tl-resize-handle"></div>' +

        // â”€â”€ Transport bar â”€â”€
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
                '<button id="tl-btn-del-preset" class="tl-btn tl-btn--del-preset" type="button" title="Remove this imported preset" style="display:none">&#128465;</button>' +
                '<span class="tl-transport-sep"></span>' +
                '<button id="tl-btn-save" class="tl-btn tl-btn--save" type="button" title="Download the current draft with the linked filename">&#128190;&nbsp;SAVE</button>' +
                '<button id="tl-btn-motion" class="tl-btn tl-btn--motion" type="button" title="Insert a predefined motion function">âŠ•&nbsp;FX</button>' +
                '<button id="tl-btn-rec" class="tl-btn tl-btn--rec" type="button" title="Recording settings">&#9679;&nbsp;REC</button>' +
                '<span class="tl-transport-sep"></span>' +
                '<button id="tl-btn-auto-key" class="tl-btn tl-btn--warn" type="button" title="Auto Keyframe: capture changes. Shift+click to skip camera.">AUTO KEY</button>' +
                '<button id="tl-btn-add-track" class="tl-btn" type="button" title="Add a new track">+ TRACK</button>' +
                '<button id="tl-btn-add-text" class="tl-btn tl-btn--text" type="button" title="Add annotation text at current time">&#9998;&nbsp;TEXT</button>' +
                '<span class="tl-transport-sep"></span>' +
                '<button id="tl-btn-import" class="tl-btn" type="button" title="Import timeline from JSON file">&#8593; IMPORT</button>' +
                '<button id="tl-btn-export" class="tl-btn" type="button" title="Export timeline as JSON file">&#8595; EXPORT</button>' +
                '<span class="tl-transport-sep"></span>' +
                '<button id="tl-btn-close" class="tl-btn" type="button" title="Close panel">&times;</button>' +
            '</div>' +
        '</div>' +

        // â”€â”€ 3-column body â”€â”€
        '<div class="tl-body">' +

            // Left: track list
            '<div class="tl-tracks">' +
                '<div class="tl-tracks-head">' +
                    '<span class="tl-tracks-title">TRACKS</span>' +
                '</div>' +
                '<div id="tl-annot-track-col" class="tl-annot-track-col"></div>' +
                '<div id="tl-track-list" class="tl-track-list"></div>' +
            '</div>' +

            // Center: dopesheet lanes
            '<div class="tl-dopesheet">' +
                '<div id="tl-zoom-ctrls" class="tl-zoom-ctrls">' +
                    '<button id="tl-zoom-out" class="tl-zoom-btn" type="button" title="Zoom out (Ctrl+scroll)">&#x2212;</button>' +
                    '<button id="tl-zoom-reset" class="tl-zoom-btn tl-zoom-level" type="button" title="Reset zoom to fit">1&#x00d7;</button>' +
                    '<button id="tl-zoom-in" class="tl-zoom-btn" type="button" title="Zoom in (Ctrl+scroll)">+</button>' +
                '</div>' +
                '<div id="tl-scroll-wrap" class="tl-scroll-wrap">' +
                    '<div id="tl-time-content" class="tl-time-content">' +
                        '<div id="tl-ruler" class="tl-ruler"></div>' +
                        '<div id="tl-annot-lane-col" class="tl-annot-lane-col"></div>' +
                        '<div id="tl-lanes" class="tl-lanes"></div>' +
                    '</div>' +
                '</div>' +
            '</div>' +

            // Right: key inspector
            '<div class="tl-inspector">' +
                '<div id="tl-key-section">' +
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
                        '<button id="tl-insp-capture" class="tl-mini-btn" type="button" title="Capture the current live value into the inspector">LIVE VALUE</button>' +
                    '</div>' +
                    '<div class="tl-insp-actions">' +
                        '<button id="tl-insp-set" class="tl-mini-btn tl-mini-btn--accent" type="button">SET KEY</button>' +
                        '<button id="tl-insp-del" class="tl-mini-btn tl-mini-btn--danger" type="button">DELETE KEY</button>' +
                    '</div>' +
                    '<datalist id="tl-path-datalist"></datalist>' +
                '</div>' +
                '<div id="tl-ev-section" class="tl-ev-section" style="display:none">' +
                    '<div class="tl-inspector-title">&#9998;&nbsp;TEXT EVENT</div>' +
                    '<div id="tl-ev-end-banner" class="tl-ev-end-banner" style="display:none">' +
                        '<span id="tl-ev-end-label">&#x23f9; End marker</span>' +
                        '<button id="tl-ev-end-remove" class="tl-mini-btn tl-mini-btn--danger" type="button" title="Remove end marker (annotation becomes permanent)">REMOVE END</button>' +
                    '</div>' +
                    '<div id="tl-ev-summary" class="tl-insp-summary">No event selected.</div>' +
                    '<div class="tl-insp-row">' +
                        '<label>Time</label>' +
                        '<input id="tl-ev-time" type="number" min="0" step="0.01" value="0">' +
                    '</div>' +
                    '<div class="tl-insp-row">' +
                        '<label>Dur.</label>' +
                        '<input id="tl-ev-dur" type="number" min="0" step="0.5" value="0" placeholder="0 = forever" title="Duration in seconds (0 = until next event)">' +
                        '<span class="tl-ev-dur-hint">s</span>' +
                    '</div>' +
                    '<div class="tl-insp-row">' +
                        '<label>Title</label>' +
                        '<input id="tl-ev-title" type="text" placeholder="Optional title">' +
                    '</div>' +
                    '<div class="tl-insp-row tl-insp-row--tall">' +
                        '<label>Text</label>' +
                        '<textarea id="tl-ev-body" rows="3" placeholder="Annotation body text..."></textarea>' +
                    '</div>' +
                    '<div class="tl-insp-row">' +
                        '<label>Color</label>' +
                        '<input id="tl-ev-color" type="color" value="#7cc5ff" class="tl-ev-color-input">' +
                        '<label style="margin-left:4px">Width</label>' +
                        '<input id="tl-ev-width" type="number" min="170" max="800" step="10" value="320" style="width:48px">' +
                    '</div>' +
                    '<div class="tl-insp-row">' +
                        '<label>Fade in</label>' +
                        '<input id="tl-ev-fadein" type="number" min="0" max="10" step="0.1" value="0" style="width:52px" title="Fade-in duration in seconds (0 = instant)">' +
                        '<span class="tl-ev-dur-hint">s</span>' +
                    '</div>' +
                    '<div class="tl-insp-row">' +
                        '<label>Place</label>' +
                        '<select id="tl-ev-placement">' +
                            '<option value="auto">Auto</option>' +
                            '<option value="right">Right</option>' +
                            '<option value="left">Left</option>' +
                            '<option value="top">Top</option>' +
                            '<option value="bottom">Bottom</option>' +
                            '<option value="manual">Manual (drag)</option>' +
                        '</select>' +
                    '</div>' +
                    '<div class="tl-insp-actions">' +
                        '<button id="tl-ev-position" class="tl-mini-btn tl-mini-btn--accent" type="button" title="Drag bubble &amp; pointer on screen">&#9995; POSITION</button>' +
                        '<button id="tl-ev-preview" class="tl-mini-btn" type="button" title="Preview this annotation now">&#128065; PREVIEW</button>' +
                    '</div>' +
                    '<div class="tl-insp-actions">' +
                        '<button id="tl-ev-use-time" class="tl-mini-btn" type="button">USE TIME</button>' +
                        '<button id="tl-ev-new" class="tl-mini-btn" type="button">+ NEW</button>' +
                    '</div>' +
                    '<div class="tl-insp-actions">' +
                        '<button id="tl-ev-set" class="tl-mini-btn tl-mini-btn--accent" type="button">SET EVENT</button>' +
                        '<button id="tl-ev-del" class="tl-mini-btn tl-mini-btn--danger" type="button">DELETE</button>' +
                    '</div>' +
                '</div>' +
            '</div>' +

        '</div>' +

        // â”€â”€ Motion Functions modal (floats above panel) â”€â”€
        '<div id="tl-motion-modal" class="tl-motion-modal">' +
            '<div class="tl-motion-hdr">' +
                '<span class="tl-motion-title">MOTION FUNCTION</span>' +
                '<button id="tl-motion-close" class="tl-btn" type="button" title="Close">&times;</button>' +
            '</div>' +
            '<div class="tl-motion-row tl-motion-type-row">' +
                '<label>Type</label>' +
                '<select id="tl-motion-type" class="tl-motion-input">' +
                    '<option value="sky_reveal">Intro: sky reveal</option>' +
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

        // â”€â”€ REC modal (floats above panel) â”€â”€
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
                    '<label class="tl-rec-check"><input type="checkbox" id="tl-rec-param-hud-show" checked> Show param values</label>' +
                    '<label class="tl-rec-check"><input type="checkbox" id="tl-rec-param-hud-record"> Include param values in recording</label>' +
                    '<label class="tl-rec-check"><input type="checkbox" id="tl-rec-reset-sim" checked> Reset simulation on rec start</label>' +
                '</div>' +
                '<div class="tl-rec-row" id="tl-hud-layout-row">' +
                    '<label>HUD pos</label>' +
                    '<input id="tl-hud-x" class="tl-rec-num" type="number" min="0" max="1" step="0.01" value="0" title="Horizontal position (0=left, 1=right)">' +
                    '<input id="tl-hud-y" class="tl-rec-num" type="number" min="0" max="1" step="0.01" value="1" title="Vertical position (0=top, 1=bottom)">' +
                    '<label style="margin-left:6px">px</label>' +
                    '<input id="tl-hud-fontsize" class="tl-rec-num" type="number" min="8" max="48" step="1" value="11" title="Font size in px">' +
                '</div>' +
                '<div class="tl-rec-section">' +
                    '<div class="tl-rec-section-head">' +
                        '<span class="tl-rec-section-title">Visible params</span>' +
                        '<div class="tl-rec-section-actions">' +
                            '<button id="tl-rec-add-selected" class="tl-mini-btn" type="button">ADD SELECTED</button>' +
                            '<button id="tl-rec-clear-hud" class="tl-mini-btn" type="button">CLEAR</button>' +
                        '</div>' +
                    '</div>' +
                    '<div id="tl-rec-hud-empty" class="tl-rec-empty">Use the HUD buttons in the track list, or add the selected track here.</div>' +
                    '<div id="tl-rec-hud-list" class="tl-rec-hud-list"></div>' +
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
                '<button id="tl-rec-shot" class="tl-mini-btn" type="button">PNG SNAPSHOT</button>' +
                '<button id="tl-rec-start" class="tl-mini-btn tl-mini-btn--accent" type="button">&#9679; START REC</button>' +
                '<button id="tl-rec-stop" class="tl-mini-btn tl-mini-btn--danger" type="button">&#9632; STOP REC</button>' +
            '</div>' +
            '<div id="tl-rec-status" class="tl-rec-status">Idle</div>' +
        '</div>' +

        // â”€â”€ Status bar â”€â”€
        '<div id="tl-status" class="tl-status"></div>';

    document.body.appendChild(panel);

    // â”€â”€ Element refs â”€â”€
    var playBtn      = panel.querySelector('#tl-btn-play');
    var pauseBtn     = panel.querySelector('#tl-btn-pause');
    var stopBtn      = panel.querySelector('#tl-btn-stop');
    var autoKeyBtn   = panel.querySelector('#tl-btn-auto-key');
    var addTrackBtn  = panel.querySelector('#tl-btn-add-track');
    var importBtn    = panel.querySelector('#tl-btn-import');
    var exportBtn    = panel.querySelector('#tl-btn-export');
    var saveBtn      = panel.querySelector('#tl-btn-save');
    var delPresetBtn = panel.querySelector('#tl-btn-del-preset');
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
    var recParamHudShowCb   = panel.querySelector('#tl-rec-param-hud-show');
    var recParamHudRecordCb = panel.querySelector('#tl-rec-param-hud-record');
    var hudXInput    = panel.querySelector('#tl-hud-x');
    var hudYInput    = panel.querySelector('#tl-hud-y');
    var hudFsInput   = panel.querySelector('#tl-hud-fontsize');
    var hudLayoutRow = panel.querySelector('#tl-hud-layout-row');
    var recHudEmptyEl    = panel.querySelector('#tl-rec-hud-empty');
    var recHudListEl     = panel.querySelector('#tl-rec-hud-list');
    var recAddSelectedBtn= panel.querySelector('#tl-rec-add-selected');
    var recClearHudBtn   = panel.querySelector('#tl-rec-clear-hud');
    var recResetSimCb    = panel.querySelector('#tl-rec-reset-sim');
    var recQualitySelect = panel.querySelector('#tl-rec-quality');
    var recModeSelect    = panel.querySelector('#tl-rec-mode');
    var recResSelect     = panel.querySelector('#tl-rec-resolution');
    var recFpsInput      = panel.querySelector('#tl-rec-fps');
    var recBitrateInput  = panel.querySelector('#tl-rec-bitrate');
    var recShotBtn       = panel.querySelector('#tl-rec-shot');
    var recStartBtn      = panel.querySelector('#tl-rec-start');
    var recStopBtn       = panel.querySelector('#tl-rec-stop');
    var recStatusEl      = panel.querySelector('#tl-rec-status');

    // â”€â”€ Event (annotation) inspector refs â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    var annotTrackColEl = panel.querySelector('#tl-annot-track-col');
    var annotLaneColEl  = panel.querySelector('#tl-annot-lane-col');
    var scrollWrapEl    = panel.querySelector('#tl-scroll-wrap');
    var timeContentEl   = panel.querySelector('#tl-time-content');
    var zoomOutBtn      = panel.querySelector('#tl-zoom-out');
    var zoomResetBtn    = panel.querySelector('#tl-zoom-reset');
    var zoomInBtn       = panel.querySelector('#tl-zoom-in');
    var keySection    = panel.querySelector('#tl-key-section');
    var evSection     = panel.querySelector('#tl-ev-section');
    var evSummary     = panel.querySelector('#tl-ev-summary');
    var evTimeInput   = panel.querySelector('#tl-ev-time');
    var evTitleInput  = panel.querySelector('#tl-ev-title');
    var evBodyInput   = panel.querySelector('#tl-ev-body');
    var evColorInput  = panel.querySelector('#tl-ev-color');
    var evWidthInput  = panel.querySelector('#tl-ev-width');
    var evFadeInput   = panel.querySelector('#tl-ev-fadein');
    var evDurInput    = panel.querySelector('#tl-ev-dur');
    var evPlacement   = panel.querySelector('#tl-ev-placement');
    var evUseTimeBtn  = panel.querySelector('#tl-ev-use-time');
    var evNewBtn      = panel.querySelector('#tl-ev-new');
    var evSetBtn      = panel.querySelector('#tl-ev-set');
    var evDelBtn      = panel.querySelector('#tl-ev-del');
    var evPositionBtn = panel.querySelector('#tl-ev-position');
    var evPreviewBtn  = panel.querySelector('#tl-ev-preview');
    var addTextBtn    = panel.querySelector('#tl-btn-add-text');
    var evEndBanner   = panel.querySelector('#tl-ev-end-banner');
    var evEndLabel    = panel.querySelector('#tl-ev-end-label');
    var evEndRemove   = panel.querySelector('#tl-ev-end-remove');

    // â”€â”€ State â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    var draft          = null;
    var selectedTrack  = '';
    var selectedKeyT   = NaN;
    // Multi-select: array of { path: string, t: number }
    var selectedKeys   = [];
    var selectedEventIdx = -1;   // index into draft.events; -1 = none
    var selectedEventChannel = 0; // which annotation channel is active in inspector
    var tlZoom = 1.0;             // timeline zoom factor (1 = fit all, >1 = zoomed in)
    var ZOOM_MIN = 1.0, ZOOM_MAX = 20.0;
    var panelOpen      = false;
    var autoKeySnapshot= null;
    var syncTimer      = null;
    var importedPresets = {};     // { presetName: true } for user-imported presets (removable)
    var linkedFileName  = null;   // download filename linked to this draft
    var restoredRecordingUiState = null;
    var PANEL_HEIGHT_KEY  = 'black-hole.tl-panel.height';
    var PANEL_STATE_KEY   = 'black-hole.tl-panel.state';
    var PANEL_DEFAULT_H   = 320;

    // â”€â”€ Clipboard (copy/paste keyframes) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    // Array of { path, t, v, ease } with anchorT = min(t) in the set
    var clipboard = null;

    // â”€â”€ Key drag state â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    // Set when the user starts dragging a diamond; reset on pointerup/cancel
    var keyDragState = null;

    // â”€â”€ Undo/Redo history â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

    // â”€â”€ Persistence â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function saveState() {
        try {
            // Also save the imported preset data so they survive reload
            var importedData = {};
            for (var ipn in importedPresets) {
                if (importedPresets.hasOwnProperty(ipn) && typeof PRESENTATION_PRESETS !== 'undefined' && PRESENTATION_PRESETS[ipn]) {
                    importedData[ipn] = clonePlain(PRESENTATION_PRESETS[ipn]);
                }
            }
            var timelineState = null;
            if (typeof getPresentationTimeline === 'function') {
                timelineState = getPresentationTimeline();
            }
            var state = {
                preset: presetSelect.value || '',
                draft: timelineState ? clonePlain(timelineState) : (draft ? clonePlain(draft) : null),
                selectedTrack: selectedTrack,
                selectedKeys: selectedKeys.slice(),
                selectedEventIdx: selectedEventIdx,
                selectedEventChannel: selectedEventChannel,
                tlZoom: tlZoom,
                tlScrollLeft: scrollWrapEl ? scrollWrapEl.scrollLeft : 0,
                currentTime: currentTime(),
                wasOpen: panelOpen,
                importedPresets: importedPresets,
                importedData: importedData,
                linkedFileName: linkedFileName,
                autoKeySnapshot: autoKeySnapshot ? clonePlain(autoKeySnapshot) : null,
                recordingUi: {
                    quality: recQualitySelect ? recQualitySelect.value : '',
                    mode: recModeSelect ? recModeSelect.value : '',
                    resolution: recResSelect ? recResSelect.value : '',
                    fps: recFpsInput ? recFpsInput.value : '',
                    bitrate: recBitrateInput ? recBitrateInput.value : '',
                    resetSimulation: recResetSimCb ? !!recResetSimCb.checked : true
                }
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

    function rememberState() {
        if (!panelOpen) return;
        saveState();
    }

    function normalizeRecordingUiState(raw) {
        if (!raw || typeof raw !== 'object') return null;
        var out = {};
        if (typeof raw.quality === 'string') out.quality = raw.quality;
        if (typeof raw.mode === 'string') out.mode = raw.mode;
        if (typeof raw.resolution === 'string') out.resolution = raw.resolution;
        if (raw.fps !== undefined) out.fps = String(raw.fps);
        if (raw.bitrate !== undefined) out.bitrate = String(raw.bitrate);
        if (typeof raw.resetSimulation === 'boolean') {
            out.resetSimulation = raw.resetSimulation;
        }
        return out;
    }

    function applySavedRecordingUiState() {
        if (!restoredRecordingUiState) return;
        if (recQualitySelect && restoredRecordingUiState.quality &&
            recQualitySelect.querySelector('option[value="' + CSS.escape(restoredRecordingUiState.quality) + '"]')) {
            recQualitySelect.value = restoredRecordingUiState.quality;
        }
        if (recModeSelect && restoredRecordingUiState.mode &&
            recModeSelect.querySelector('option[value="' + CSS.escape(restoredRecordingUiState.mode) + '"]')) {
            recModeSelect.value = restoredRecordingUiState.mode;
        }
        if (recResSelect && restoredRecordingUiState.resolution &&
            recResSelect.querySelector('option[value="' + CSS.escape(restoredRecordingUiState.resolution) + '"]')) {
            recResSelect.value = restoredRecordingUiState.resolution;
        }
        if (recFpsInput && restoredRecordingUiState.fps) recFpsInput.value = restoredRecordingUiState.fps;
        if (recBitrateInput && restoredRecordingUiState.bitrate) recBitrateInput.value = restoredRecordingUiState.bitrate;
        if (recResetSimCb && typeof restoredRecordingUiState.resetSimulation === 'boolean') {
            recResetSimCb.checked = restoredRecordingUiState.resetSimulation;
        }
    }

    // â”€â”€ Panel open/close â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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
            if (saved) {
                // Restore imported presets first so they are available in the dropdown
                if (saved.importedPresets && typeof saved.importedPresets === 'object') {
                    importedPresets = saved.importedPresets;
                    var iData = saved.importedData || {};
                    for (var ip in importedPresets) {
                        if (importedPresets.hasOwnProperty(ip) && typeof registerPresentationPreset === 'function') {
                            if (iData[ip]) {
                                registerPresentationPreset(clonePlain(iData[ip]), ip);
                            }
                        }
                    }
                }
                if (saved.linkedFileName) linkedFileName = saved.linkedFileName;
                restoredRecordingUiState = normalizeRecordingUiState(saved.recordingUi);
            }
            populatePresets();
            if (saved) {
                // Restore persisted state (saved.preset may be '' for the "new empty" option)
                if (typeof saved.preset === 'string' && presetSelect.querySelector('option[value="' + CSS.escape(saved.preset) + '"]')) {
                    presetSelect.value = saved.preset;
                }
                if (saved.draft) {
                    draft = normalizeTL(saved.draft);
                    if (typeof setPresentationTimeline === 'function') {
                        applyingDraft = true;
                        setPresentationTimeline(clonePlain(draft));
                        applyingDraft = false;
                    }
                    if (saved.autoKeySnapshot && typeof saved.autoKeySnapshot === 'object' &&
                        typeof saved.autoKeySnapshot.values === 'object') {
                        autoKeySnapshot = clonePlain(saved.autoKeySnapshot);
                    }
                    if (typeof saved.currentTime === 'number' && isFinite(saved.currentTime) &&
                        typeof seekPresentation === 'function') {
                        seekPresentation(clamp(saved.currentTime, 0, Math.max(draft.duration || 0, 0.001)));
                    }
                    rebuildAll();
                } else {
                    syncFromRuntime();
                }
                selectedTrack = '';
                selectedKeys = [];
                selectedKeyT = NaN;
                selectedEventIdx = -1;
                if (saved.selectedTrack) selectedTrack = saved.selectedTrack;
                if (saved.selectedKeys && saved.selectedKeys.length) selectedKeys = saved.selectedKeys;
                if (selectedKeys.length === 1) selectedKeyT = selectedKeys[0].t;
                if (typeof saved.selectedEventIdx === 'number' && saved.selectedEventIdx >= 0) selectedEventIdx = saved.selectedEventIdx;
                if (typeof saved.selectedEventChannel === 'number') selectedEventChannel = saved.selectedEventChannel;
                if (typeof saved.tlZoom === 'number' && saved.tlZoom >= 1.0) {
                    tlZoom = saved.tlZoom;
                    applyZoomWidth();
                    updateZoomDisplay();
                    if (typeof saved.tlScrollLeft === 'number') {
                        scrollWrapEl.scrollLeft = saved.tlScrollLeft;
                    }
                }
                rebuildAll();
                updateDelPresetBtn();
            } else {
                // No saved state â€” default to a fresh empty timeline
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

    // â”€â”€ Panel resize â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    var storedH = null;
    try { storedH = parseFloat(localStorage.getItem(PANEL_HEIGHT_KEY)); } catch(e) {}
    if (isFinite(storedH)) panel.style.height = clamp(storedH, PANEL_MIN_H, PANEL_MAX_H) + 'px';
    else panel.style.height = PANEL_DEFAULT_H + 'px';

    window.addEventListener('pagehide', saveState);
    window.addEventListener('beforeunload', saveState);

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

    // â”€â”€ Scrubber â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

    // â”€â”€ Transport â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

    // â”€â”€ Editable time / duration inputs â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

    // â”€â”€ Preset selector â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function populatePresets() {
        var names = (typeof listPresentationPresets === 'function')
            ? listPresentationPresets() : [];
        var currentVal = presetSelect.value;
        presetSelect.innerHTML = '';

        // "New (empty)" option
        var emptyOpt = document.createElement('option');
        emptyOpt.value = '';
        emptyOpt.textContent = '\u2014 new empty \u2014';
        presetSelect.appendChild(emptyOpt);

        for (var i = 0; i < names.length; i++) {
            var opt = document.createElement('option');
            opt.value = names[i];
            opt.textContent = (importedPresets[names[i]] ? '\u2606 ' : '') + names[i];
            presetSelect.appendChild(opt);
        }

        // Restore previous selection if possible (currentVal may be '' for the "new empty" option)
        if (presetSelect.querySelector('option[value="' + CSS.escape(currentVal) + '"]')) {
            presetSelect.value = currentVal;
        } else if (names.length) {
            var st = typeof getPresentationState === 'function' ? getPresentationState() : null;
            var loadedName = st && st.name ? st.name : '';
            if (loadedName && names.indexOf(loadedName) !== -1) {
                presetSelect.value = loadedName;
            } else {
                presetSelect.value = '';
            }
        }
        updateDelPresetBtn();
    }

    function updateDelPresetBtn() {
        var name = presetSelect.value;
        delPresetBtn.style.display = (name && importedPresets[name]) ? '' : 'none';
    }

    function loadPresetByName(name) {
        if (!name) {
            // Empty = start fresh
            draft = normalizeTL({ name: 'Untitled', duration: 12, tracks: [], events: [] });
            if (typeof setPresentationTimeline === 'function') {
                applyingDraft = true;
                setPresentationTimeline(clonePlain(draft));
                applyingDraft = false;
            }
            linkedFileName = null;
            autoKeySnapshot = null;
            rebuildAll();
            updateDelPresetBtn();
            setStatus('New empty timeline created.', '');
            return;
        }
        if (typeof loadPresentationPreset === 'function') {
            if (loadPresentationPreset(name)) {
                syncFromRuntime();
                resetZoom();
                linkedFileName = null;
                autoKeySnapshot = null;
                updateDelPresetBtn();
                setStatus('Loaded preset: ' + name, '');
            } else {
                setStatus('Failed to load preset: ' + name, 'tl-status--warn');
            }
        }
    }

    presetSelect.addEventListener('change', function() {
        loadPresetByName(presetSelect.value);
    });

    // â”€â”€ Sync from runtime â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

    // â”€â”€ Inspector mode toggles â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function showKeyInspector() {
        keySection.style.display = '';
        evSection.style.display  = 'none';
    }
    function showEventInspector(isClearEvent) {
        keySection.style.display = 'none';
        evSection.style.display  = '';
        // Show/hide the end-marker banner and form fields based on event type
        var isEnd = !!isClearEvent;
        evEndBanner.style.display  = isEnd ? '' : 'none';
        // Disable editing fields for end markers so the user can't accidentally re-save
        var fields = [evTimeInput, evTitleInput, evBodyInput, evColorInput, evWidthInput,
                      evFadeInput, evDurInput, evPlacement, evUseTimeBtn, evSetBtn,
                      evPositionBtn, evPreviewBtn];
        for (var _fi = 0; _fi < fields.length; _fi++) {
            if (fields[_fi]) fields[_fi].disabled = isEnd;
        }
    }

    // â”€â”€ Events lane (center) + row label (left) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function rebuildEventsLane() {
        if (!annotTrackColEl || !annotLaneColEl) return;
        var d = getDuration();
        var phPct = clamp(currentTime() / Math.max(d, 0.001) * 100, 0, 100);
        var annTracks = (draft && Array.isArray(draft.annotationTracks)) ? draft.annotationTracks : [{ label: 'Annotation 1' }];
        var numCh = annTracks.length;
        var canDelete = numCh > 1;

        var leftHtml = '<div class="tl-ann-hdr">' +
            '<span class="tl-ann-hdr-title">&#9998;&nbsp;ANNOTATIONS</span>' +
            '<button class="tl-ann-add-btn" id="tl-annot-add-ch" title="Add annotation track" type="button">+</button>' +
            '</div>';
        var rightHtml = '<div class="tl-ann-hdr-spacer"></div>';

        for (var ch = 0; ch < numCh; ch++) {
            var label = annTracks[ch].label || ('Annotation ' + (ch + 1));
            var isSelCh = (ch === selectedEventChannel && selectedEventIdx >= 0);
            leftHtml += '<div class="tl-events-row' + (isSelCh ? ' is-sel-ch' : '') + '" data-ch="' + ch + '">' +
                '<span class="tl-events-row-title">' + esc(label) + '</span>' +
                (canDelete ? '<button class="tl-ann-del-ch" data-ch="' + ch + '" title="Remove this annotation track" type="button">&times;</button>' : '') +
                '</div>';

            var laneHtml = '<span class="tl-lane-playhead" style="left:' + phPct.toFixed(2) + '%"></span>';
            var hasMarkers = false;
            if (draft) {
                for (var i = 0; i < draft.events.length; i++) {
                    var ev = draft.events[i];
                    if (ev.action !== 'annotation' && ev.action !== 'clearAnnotation') continue;
                    if ((ev.channel || 0) !== ch) continue;
                    hasMarkers = true;
                    var kPct = clamp(ev.t / Math.max(d, 0.001) * 100, 0, 100);
                    var isSel = (i === selectedEventIdx);
                    var isClear = (ev.action === 'clearAnnotation');
                    var lbl = (ev.note && ev.note.title) ? ev.note.title : (isClear ? '\u2715 clear' : '(no title)');
                    var color = (ev.note && ev.note.color) ? ev.note.color : '#7cc5ff';
                    var tip = ev.t.toFixed(2) + 's: ' + lbl;
                    var borderStyle = (!isSel && !isClear) ? 'border-color:' + esc(color) + ';' : '';
                    // Compute annotation bar width: end = earliest of paired clearAnnotation,
                    // any clearAnnotation on same channel, or next annotation on same channel
                    // (implicit replace â€” each new annotation fires on the same channel).
                    var durPct = 0;
                    if (ev.action === 'annotation') {
                        var endT = Infinity;
                        for (var ci = 0; ci < draft.events.length; ci++) {
                            var ce = draft.events[ci];
                            if (ce.action === 'clearAnnotation' && ce._pairOf === i && ce.t > ev.t) {
                                endT = Math.min(endT, ce.t); break;
                            }
                        }
                        if (!isFinite(endT)) {
                            for (var ci = 0; ci < draft.events.length; ci++) {
                                var ce = draft.events[ci];
                                if (ce.action === 'clearAnnotation' && (ce.channel || 0) === ch && ce.t > ev.t) {
                                    endT = Math.min(endT, ce.t); break;
                                }
                            }
                        }
                        // Next annotation on same channel is an implicit clear
                        for (var ci = 0; ci < draft.events.length; ci++) {
                            var ce = draft.events[ci];
                            if (ce.action === 'annotation' && (ce.channel || 0) === ch && ce.t > ev.t) {
                                endT = Math.min(endT, ce.t); break;
                            }
                        }
                        if (isFinite(endT)) {
                            durPct = (endT - ev.t) / Math.max(d, 0.001) * 100;
                        }
                    }
                    var widthStyle = durPct > 0 ? 'width:' + durPct.toFixed(3) + '%;' : '';
                    var hasDur = durPct > 0 ? ' has-dur' : '';
                    laneHtml += '<button type="button" class="tl-ev-marker' +
                        (isSel ? ' is-sel' : '') + (isClear ? ' is-clear' : '') + hasDur + '"' +
                        ' style="left:' + kPct.toFixed(3) + '%;' + widthStyle + borderStyle + '"' +
                        ' data-ei="' + i + '" data-ch="' + ch + '" title="' + esc(tip) + '"></button>';
                }
            }
            if (!hasMarkers) {
                laneHtml += '<span class="tl-ev-hint">Empty</span>';
            }
            rightHtml += '<div class="tl-events-lane" data-ch="' + ch + '">' + laneHtml + '</div>';
        }

        annotTrackColEl.innerHTML = leftHtml;
        annotLaneColEl.innerHTML  = rightHtml;

        // Wire the dynamically-injected add/del buttons
        var addChBtn = annotTrackColEl.querySelector('#tl-annot-add-ch');
        if (addChBtn) addChBtn.addEventListener('click', addAnnotationChannel);
        var delBtns = annotTrackColEl.querySelectorAll('.tl-ann-del-ch');
        for (var b = 0; b < delBtns.length; b++) {
            (function(btn) {
                btn.addEventListener('click', function() {
                    removeAnnotationChannel(parseInt(btn.getAttribute('data-ch'), 10));
                });
            })(delBtns[b]);
        }
    }

    // â”€â”€ Add / remove annotation channels â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function addAnnotationChannel() {
        if (!draft) { draft = normalizeTL({ name: 'Untitled', duration: 12, tracks: [], events: [] }); }
        pushUndo();
        if (!Array.isArray(draft.annotationTracks)) draft.annotationTracks = [{ label: 'Annotation 1' }];
        // Find a unique label that doesn't collide with any existing track
        var existingLabels = {};
        for (var _li = 0; _li < draft.annotationTracks.length; _li++) {
            existingLabels[draft.annotationTracks[_li].label] = true;
        }
        var n = draft.annotationTracks.length + 1;
        while (existingLabels['Annotation ' + n]) n++;
        draft.annotationTracks.push({ label: 'Annotation ' + n });
        selectedEventChannel = draft.annotationTracks.length - 1;
        applyDraft();
        rebuildAll();
        setStatus('Annotation track ' + n + ' added.', '');
    }

    function removeAnnotationChannel(ch) {
        if (!draft || !Array.isArray(draft.annotationTracks) || draft.annotationTracks.length <= 1) {
            setStatus('Cannot remove the last annotation track.', 'tl-status--warn');
            return;
        }
        pushUndo();
        // Remove all events on this channel
        for (var i = draft.events.length - 1; i >= 0; i--) {
            if ((draft.events[i].channel || 0) === ch) draft.events.splice(i, 1);
        }
        // Shift channels above the removed one
        for (var j = 0; j < draft.events.length; j++) {
            if (typeof draft.events[j].channel === 'number' && draft.events[j].channel > ch) {
                draft.events[j].channel--;
            }
        }
        draft.annotationTracks.splice(ch, 1);
        if (selectedEventChannel >= draft.annotationTracks.length) {
            selectedEventChannel = Math.max(0, draft.annotationTracks.length - 1);
        }
        if (selectedEventIdx >= 0 && !draft.events[selectedEventIdx]) {
            selectedEventIdx = -1;
            showKeyInspector();
        }
        // Reindex _pairOf after events were spliced out, before applying
        reindexAllPairOf(draft.events);
        applyDraft();
        rebuildAll();
        setStatus('Annotation track removed.', '');
    }

    // â”€â”€ _pairOf re-indexing after sort â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    // Rebuilds all _pairOf indices by looking up annotation objects that each
    // clearAnnotation's _pairOf previously pointed to.  Must be called AFTER
    // the events array has been sorted (and before any further _pairOf reads).
    function reindexAllPairOf(events) {
        for (var i = 0; i < events.length; i++) {
            var e = events[i];
            if (e.action !== 'clearAnnotation') continue;
            // Process all clearAnnotation events, even those without a _pairOf yet
            // _pairOf might be stale; just leave it if we can't find anything better.
            // Heuristic: among annotation events on the same channel that start
            // before this clear, pick the nearest one.
            var bestIdx = -1, bestT = -Infinity;
            var ch = (typeof e.channel === 'number') ? e.channel : 0;
            for (var j = 0; j < events.length; j++) {
                if (events[j].action === 'annotation' &&
                    ((events[j].channel || 0) === ch) &&
                    events[j].t < e.t && events[j].t > bestT) {
                    bestT = events[j].t;
                    bestIdx = j;
                }
            }
            if (bestIdx >= 0) e._pairOf = bestIdx;
        }
    }

    // â”€â”€ Event inspector fill â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function fillEventInspector(ev, idx) {
        selectedEventIdx = idx;
        selectedEventChannel = ev.channel || 0;
        var isClear = (ev.action === 'clearAnnotation');
        showEventInspector(isClear);
        evTimeInput.value  = ev.t.toFixed(2);
        var note = ev.note || {};
        evTitleInput.value = note.title || '';
        evBodyInput.value  = note.text || note.body || '';
        evColorInput.value = note.color || '#7cc5ff';
        evWidthInput.value = note.width || 320;
        evFadeInput.value  = (note.fadeIn && note.fadeIn > 0) ? note.fadeIn : 0;
        evPlacement.value  = (typeof note.boxX === 'number' && typeof note.boxY === 'number')
            ? 'manual' : (note.placement || 'auto');

        // For end markers, find the paired annotation title for the banner
        if (isClear) {
            var pairedTitle = '';
            if (draft && typeof ev._pairOf === 'number' && draft.events[ev._pairOf]
                    && draft.events[ev._pairOf].note) {
                pairedTitle = draft.events[ev._pairOf].note.title || '';
            }
            evEndLabel.textContent = pairedTitle
                ? '\u23f9 End marker for: \u201c' + pairedTitle + '\u201d'
                : '\u23f9 End marker (clears annotation)';
        }

        // Compute duration from the explicitly paired clearAnnotation event only.
        // Never fall back to an arbitrary clear event â€” that causes phantom durations.
        var dur = 0;
        if (draft && ev.action === 'annotation') {
            for (var i = 0; i < draft.events.length; i++) {
                var ce = draft.events[i];
                if (ce.action === 'clearAnnotation' && ce.t > ev.t && ce._pairOf === idx) {
                    dur = ce.t - ev.t;
                    break;
                }
            }
        }
        evDurInput.value = dur > 0 ? dur.toFixed(2) : '';

        var isClear = (ev.action === 'clearAnnotation');
        var trackLabel = (draft && Array.isArray(draft.annotationTracks) && draft.annotationTracks[selectedEventChannel])
            ? draft.annotationTracks[selectedEventChannel].label : ('Track ' + (selectedEventChannel + 1));
        evSummary.textContent = isClear
            ? 'Clear [' + trackLabel + '] @ ' + ev.t.toFixed(2) + 's'
            : '[' + trackLabel + '] @ ' + ev.t.toFixed(2) + 's' + (note.title ? ': ' + note.title : '');
    }

    // â”€â”€ Set / save event from inspector fields â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function doSetEvent() {
        if (!draft) return;
        var t = parseTime(evTimeInput.value, currentTime());
        t = Math.max(0, t);
        var title     = evTitleInput.value.trim();
        var body      = evBodyInput.value;
        var color     = evColorInput.value || '#7cc5ff';
        var placement = evPlacement.value  || 'auto';
        var widthVal  = parseInt(evWidthInput.value, 10);
        if (!isFinite(widthVal) || widthVal < 100) widthVal = 320;
        var durVal    = parseFloat(evDurInput.value);
        if (!isFinite(durVal) || durVal < 0) durVal = 0;

        pushUndo();
        var ev;
        var note = null;
        if (title || body.trim()) {
            // Preserve existing anchor / boxX / boxY if editing
            var prevNote = (selectedEventIdx >= 0 && draft.events[selectedEventIdx])
                ? draft.events[selectedEventIdx].note || {} : {};
            note = {
                anchor: prevNote.anchor || { mode: 'world', target: 'black_hole' },
                placement: placement === 'manual' ? (prevNote.placement || 'auto') : placement,
                color: color,
                width: widthVal
            };
            if (placement === 'manual') {
                if (typeof prevNote.boxX === 'number') note.boxX = prevNote.boxX;
                if (typeof prevNote.boxY === 'number') note.boxY = prevNote.boxY;
            }
            if (title) note.title = title;
            if (body.trim()) note.text = body;
            var fadeInVal = parseFloat(evFadeInput.value);
            if (isFinite(fadeInVal) && fadeInVal > 0) note.fadeIn = fadeInVal;
        } else if (selectedEventIdx >= 0 && draft.events[selectedEventIdx]) {
            // Fields are blank â€” preserve the existing note so that SET EVENT with an empty
            // title/body (e.g. the user only changed the time) never silently converts an
            // annotation into a clearAnnotation.
            note = draft.events[selectedEventIdx].note || null;
        }
        var action = note ? 'annotation' : 'clearAnnotation';

        // Save the object reference of the currently paired clear event (if any)
        // so we can find it reliably after sorting.
        var oldPairedClear = null;
        if (selectedEventIdx >= 0 && draft.events[selectedEventIdx]) {
            var origEv = draft.events[selectedEventIdx];
            for (var i = 0; i < draft.events.length; i++) {
                if (draft.events[i].action === 'clearAnnotation' &&
                    draft.events[i]._pairOf === selectedEventIdx) {
                    oldPairedClear = draft.events[i];
                    break;
                }
            }
        }

        if (selectedEventIdx >= 0 && draft.events[selectedEventIdx]) {
            ev = draft.events[selectedEventIdx];
            ev.t      = t;
            ev.action = action;
            if (note) { ev.note = note; } else { delete ev.note; }
            draft.events.sort(function(a, b) { return a.t - b.t; });
            for (var i = 0; i < draft.events.length; i++) {
                if (draft.events[i] === ev) { selectedEventIdx = i; break; }
            }
        } else {
            ev = { t: t, action: action, channel: selectedEventChannel };
            if (note) ev.note = note;
            draft.events.push(ev);
            draft.events.sort(function(a, b) { return a.t - b.t; });
            for (var i = 0; i < draft.events.length; i++) {
                if (draft.events[i] === ev) { selectedEventIdx = i; break; }
            }
        }

        // â”€â”€ Duration â†’ auto-manage clearAnnotation event â”€â”€
        if (action === 'annotation') {
            // Remove any existing paired clearAnnotation (found by object reference)
            if (oldPairedClear) {
                for (var i = draft.events.length - 1; i >= 0; i--) {
                    if (draft.events[i] === oldPairedClear) {
                        draft.events.splice(i, 1);
                        if (i < selectedEventIdx) selectedEventIdx--;
                        break;
                    }
                }
            }
            var newClearEv = null;
            if (durVal > 0) {
                newClearEv = { t: t + durVal, action: 'clearAnnotation', channel: ev.channel || 0 };
                draft.events.push(newClearEv);
            }
            draft.events.sort(function(a, b) { return a.t - b.t; });
            // Re-find selected index after sort
            for (var i = 0; i < draft.events.length; i++) {
                if (draft.events[i] === ev) { selectedEventIdx = i; break; }
            }
            // Set _pairOf on the newly created clear event only
            if (newClearEv) {
                newClearEv._pairOf = selectedEventIdx;
            }
            // Re-index all other _pairOf references using object identity
            reindexAllPairOf(draft.events);
        }

        draft.duration = Math.max(draft.duration, t + (durVal || 0));
        applyDraft();
        rebuildAll();
        setStatus((action === 'annotation' ? 'Annotation' : 'Clear-annotation') + ' set @ ' + t.toFixed(2) + 's.', '');
    }

    // â”€â”€ Add new annotation event â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function addNewAnnotationEvent() {
        if (!draft) {
            draft = normalizeTL({ name: 'Untitled', duration: 12, tracks: [], events: [] });
        }
        var t = currentTime();
        pushUndo();
        var ev = { t: t, action: 'annotation', channel: selectedEventChannel,
                   note: { title: '', text: '', anchor: { mode: 'world', target: 'black_hole' },
                           placement: 'auto', color: '#7cc5ff' } };
        draft.events.push(ev);
        draft.events.sort(function(a, b) { return a.t - b.t; });
        for (var i = 0; i < draft.events.length; i++) {
            if (draft.events[i] === ev) { selectedEventIdx = i; break; }
        }
        clearMultiSelect();
        applyDraft();
        rebuildAll();
        fillEventInspector(draft.events[selectedEventIdx], selectedEventIdx);
        setStatus('Annotation added at ' + t.toFixed(2) + 's \u2013 edit title/text and click SET EVENT.', 'tl-status--info');
    }

    // â”€â”€ Auto-sync inspector fields into draft event (no undo push) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function syncInspectorToDraft() {
        if (!draft || selectedEventIdx < 0 || !draft.events[selectedEventIdx]) return;
        var ev = draft.events[selectedEventIdx];
        if (ev.action !== 'annotation') return;
        if (!ev.note) ev.note = {};
        var title = evTitleInput.value.trim();
        var body  = evBodyInput.value;
        if (title) ev.note.title = title; else delete ev.note.title;
        if (body.trim()) ev.note.text = body; else delete ev.note.text;
        ev.note.color = evColorInput.value || '#7cc5ff';
        var w = parseInt(evWidthInput.value, 10);
        ev.note.width = (isFinite(w) && w >= 100) ? w : 320;
        var fi = parseFloat(evFadeInput.value);
        if (isFinite(fi) && fi > 0) ev.note.fadeIn = fi; else delete ev.note.fadeIn;
    }

    // â”€â”€ Preview annotation live â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function previewAnnotation() {
        if (typeof setPresentationAnnotation !== 'function') return;
        var note = buildNoteFromInspector();
        if (!note) { if (typeof clearPresentationAnnotation === 'function') clearPresentationAnnotation(); return; }
        setPresentationAnnotation(note);
    }

    function buildNoteFromInspector() {
        var title = evTitleInput.value.trim();
        var body  = evBodyInput.value;
        if (!title && !body.trim()) return null;
        var widthVal = parseInt(evWidthInput.value, 10);
        if (!isFinite(widthVal) || widthVal < 100) widthVal = 320;
        var placement = evPlacement.value || 'auto';
        var prevNote = (selectedEventIdx >= 0 && draft && draft.events[selectedEventIdx])
            ? draft.events[selectedEventIdx].note || {} : {};
        var note = {
            anchor: prevNote.anchor || { mode: 'world', target: 'black_hole' },
            placement: placement === 'manual' ? (prevNote.placement || 'auto') : placement,
            color: evColorInput.value || '#7cc5ff',
            width: widthVal
        };
        if (placement === 'manual') {
            if (typeof prevNote.boxX === 'number') note.boxX = prevNote.boxX;
            if (typeof prevNote.boxY === 'number') note.boxY = prevNote.boxY;
        }
        if (title) note.title = title;
        if (body.trim()) note.text = body;
        var fadeIn = parseFloat(evFadeInput.value);
        if (isFinite(fadeIn) && fadeIn > 0) note.fadeIn = fadeIn;
        return note;
    }

    // â”€â”€ Drag overlay for positioning box and pointer â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function startPositionDrag() {
        if (typeof setPresentationAnnotation !== 'function') return;
        syncInspectorToDraft();
        var note = buildNoteFromInspector();
        if (!note) { setStatus('Add title or text first.', 'tl-status--warn'); return; }

        // Create fullscreen overlay
        var overlay = document.createElement('div');
        overlay.className = 'tl-drag-overlay';
        overlay.innerHTML =
            '<div class="tl-drag-help">Drag <b>box</b> to position bubble. Drag <b>circle</b> to position pointer tip. Press <b>Esc</b> or <b>right-click</b> to finish.</div>' +
            '<div class="tl-drag-handle tl-drag-handle--box" title="Drag to move bubble"></div>' +
            '<div class="tl-drag-handle tl-drag-handle--ptr" title="Drag to move pointer tip"></div>';
        document.body.appendChild(overlay);

        // Hide the timeline panel so it doesn't obstruct the full-screen drag area
        panel.classList.add('tl-panel--collapsed');
        document.body.classList.remove('has-timeline-panel');

        var boxHandle = overlay.querySelector('.tl-drag-handle--box');
        var ptrHandle = overlay.querySelector('.tl-drag-handle--ptr');

        var viewW = window.innerWidth;
        var viewH = window.innerHeight;

        // Initial positions
        var boxX = (typeof note.boxX === 'number') ? note.boxX : 0.7;
        var boxY = (typeof note.boxY === 'number') ? note.boxY : 0.3;
        var anchorX = (note.anchor && typeof note.anchor.x === 'number') ? note.anchor.x : 0.5;
        var anchorY = (note.anchor && typeof note.anchor.y === 'number') ? note.anchor.y : 0.5;

        function updateHandles() {
            boxHandle.style.left = (boxX * viewW) + 'px';
            boxHandle.style.top  = (boxY * viewH) + 'px';
            ptrHandle.style.left = (anchorX * viewW) + 'px';
            ptrHandle.style.top  = (anchorY * viewH) + 'px';
        }

        function updateLivePreview() {
            note.boxX = boxX;
            note.boxY = boxY;
            note.anchor = { mode: 'screen', x: anchorX, y: anchorY };
            note.placement = note.placement || 'auto';
            setPresentationAnnotation(note);
        }

        updateHandles();
        updateLivePreview();

        var dragging = null; // 'box' or 'ptr'
        var dragOffX = 0, dragOffY = 0;

        function onDown(e) {
            e.preventDefault();
            var target = e.target;
            if (target === boxHandle) {
                dragging = 'box';
            } else if (target === ptrHandle) {
                dragging = 'ptr';
            } else {
                return;
            }
            var rect = target.getBoundingClientRect();
            dragOffX = e.clientX - rect.left - rect.width / 2;
            dragOffY = e.clientY - rect.top - rect.height / 2;
        }

        function onMove(e) {
            if (!dragging) return;
            e.preventDefault();
            var x = (e.clientX - dragOffX) / viewW;
            var y = (e.clientY - dragOffY) / viewH;
            x = Math.max(0, Math.min(1, x));
            y = Math.max(0, Math.min(1, y));
            if (dragging === 'box') {
                boxX = x; boxY = y;
            } else {
                anchorX = x; anchorY = y;
            }
            updateHandles();
            updateLivePreview();
        }

        function onUp() {
            dragging = null;
        }

        function finish() {
            overlay.removeEventListener('pointerdown', onDown);
            overlay.removeEventListener('pointermove', onMove);
            overlay.removeEventListener('pointerup', onUp);
            document.removeEventListener('keydown', onKey);
            overlay.removeEventListener('contextmenu', onContext);
            document.body.removeChild(overlay);
            // Restore the timeline panel
            if (panelOpen) {
                panel.classList.remove('tl-panel--collapsed');
                document.body.classList.add('has-timeline-panel');
                updatePushedOffset();
            }
            if (typeof clearPresentationAnnotation === 'function') clearPresentationAnnotation();

            // Apply to draft â€” save the full note (title/text/color from form + position from drag)
            if (draft && selectedEventIdx >= 0 && draft.events[selectedEventIdx]) {
                pushUndo();
                var ev = draft.events[selectedEventIdx];
                // Build complete note from inspector fields + drag position
                var fullNote = buildNoteFromInspector() || {};
                fullNote.boxX = boxX;
                fullNote.boxY = boxY;
                fullNote.anchor = { mode: 'screen', x: anchorX, y: anchorY };
                ev.note = fullNote;
                evPlacement.value = 'manual';
                applyDraft();
                rebuildAll();
                fillEventInspector(draft.events[selectedEventIdx], selectedEventIdx);
                setStatus('Position saved.', 'tl-status--info');
            }
        }

        function onKey(e) {
            if (e.key === 'Escape') { e.preventDefault(); e.stopPropagation(); finish(); }
        }

        function onContext(e) {
            e.preventDefault();
            finish();
        }

        overlay.addEventListener('pointerdown', onDown);
        overlay.addEventListener('pointermove', onMove);
        overlay.addEventListener('pointerup', onUp);
        document.addEventListener('keydown', onKey, true);
        overlay.addEventListener('contextmenu', onContext);
    }

    function rebuildAll() {
        rebuildTrackList();
        rebuildLanes();
        rebuildEventsLane();
        updateRuler();
        updateInspector();
        updatePathDatalist();
    }

    // â”€â”€ Timeline zoom & scroll â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function updateZoomDisplay() {
        if (!zoomResetBtn) return;
        if (tlZoom <= 1.005) {
            zoomResetBtn.textContent = '1\u00d7';
        } else {
            zoomResetBtn.textContent = (tlZoom < 10 ? tlZoom.toFixed(1) : Math.round(tlZoom)) + '\u00d7';
        }
    }

    function applyZoomWidth() {
        if (!timeContentEl) return;
        timeContentEl.style.width = (tlZoom * 100).toFixed(3) + '%';
    }

    function zoomAtClientX(newZoom, clientX) {
        newZoom = clamp(newZoom, ZOOM_MIN, ZOOM_MAX);
        var oldZoom = tlZoom;
        if (Math.abs(newZoom - oldZoom) < 1e-4) return;
        var rect = scrollWrapEl.getBoundingClientRect();
        var cx = (clientX !== undefined) ? clamp(clientX - rect.left, 0, rect.width) : rect.width * 0.5;
        var oldScroll = scrollWrapEl.scrollLeft;
        tlZoom = newZoom;
        applyZoomWidth();
        // Keep the time-point under the cursor at the same screen position
        var ratio = oldZoom / newZoom;
        var newScroll = (oldScroll + cx) * (newZoom / oldZoom) - cx;
        var maxScroll = rect.width * (newZoom - 1);
        scrollWrapEl.scrollLeft = clamp(newScroll, 0, maxScroll);
        updateZoomDisplay();
        updateRuler();
        updatePlayheads();
    }

    function resetZoom() {
        tlZoom = 1.0;
        applyZoomWidth();
        scrollWrapEl.scrollLeft = 0;
        updateZoomDisplay();
        updateRuler();
        updatePlayheads();
    }

    // Scroll the view so the current playhead is visible
    function scrollToPlayhead() {
        if (!scrollWrapEl || tlZoom <= 1.0) return;
        var t = currentTime(), d = getDuration();
        var viewW = scrollWrapEl.clientWidth;
        var contentW = viewW * tlZoom;
        var phPx = (t / Math.max(d, 0.001)) * contentW;
        var sl = scrollWrapEl.scrollLeft;
        var margin = viewW * 0.12;
        if (phPx > sl + viewW - margin) {
            scrollWrapEl.scrollLeft = clamp(phPx - viewW + margin, 0, contentW - viewW);
        } else if (phPx < sl + margin) {
            scrollWrapEl.scrollLeft = clamp(phPx - margin, 0, contentW - viewW);
        }
    }

    // Ctrl+wheel â†’ zoom; plain wheel when zoomed + deltaX=0 â†’ horizontal scroll
    scrollWrapEl.addEventListener('wheel', function(e) {
        if (e.ctrlKey || e.metaKey) {
            e.preventDefault();
            var factor = e.deltaY > 0 ? (1.0 / 1.25) : 1.25;
            zoomAtClientX(tlZoom * factor, e.clientX);
        } else if (tlZoom > 1.0 && e.deltaX === 0 && e.deltaY !== 0) {
            // Regular mouse wheel: redirect vertical scroll to horizontal when zoomed
            e.preventDefault();
            scrollWrapEl.scrollLeft += e.deltaY;
        }
    }, { passive: false });

    zoomOutBtn.addEventListener('click', function() { zoomAtClientX(tlZoom / 1.5); });
    zoomInBtn.addEventListener('click', function() { zoomAtClientX(tlZoom * 1.5); });
    zoomResetBtn.addEventListener('click', resetZoom);

    function startSync() {
        if (syncTimer) return;
        syncTimer = setInterval(function() {
            updateTimeInputs();
            updateScrubber();
            updatePlayheads();
            syncRecModal();
            // Auto-scroll to keep playhead visible during playback
            if (tlZoom > 1.0) {
                var pst = typeof getPresentationState === 'function' ? getPresentationState() : null;
                if (pst && pst.active && !pst.paused) scrollToPlayhead();
            }
        }, 80);
    }
    function stopSync() { clearInterval(syncTimer); syncTimer = null; }

    // â”€â”€ Timeline data helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function normalizeTL(raw) {
        return normalizeTimelineDraft(raw, {
            annotationsState: (typeof getPresentationAnnotationsState === 'function')
                ? getPresentationAnnotationsState()
                : null,
            paramHudState: (typeof getPresentationParamHudState === 'function')
                ? getPresentationParamHudState()
                : null
        });
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

    var applyingDraft = false;
    function applyDraft() {
        if (!draft) return;
        if (typeof setPresentationTimeline === 'function') {
            applyingDraft = true;
            setPresentationTimeline(clonePlain(draft));
            applyingDraft = false;
            setStatus('Timeline applied.', '');
            rememberState();
        }
    }

    function ensureDraftForExternalInsert() {
        if (draft) return draft;
        var runtimeTimeline = (typeof getPresentationTimeline === 'function')
            ? getPresentationTimeline()
            : null;
        if (runtimeTimeline) {
            draft = normalizeTL(runtimeTimeline);
        } else {
            draft = normalizeTL({
                name: 'Captured Motion',
                duration: Math.max(12, currentTime() + 6),
                tracks: [],
                events: []
            });
        }
        return draft;
    }

    function insertAnimationCapture(capture) {
        if (!capture || typeof capture !== 'object') {
            return { ok: false, error: 'Invalid animation capture.' };
        }

        var mode = (capture.mode === 'hover') ? 'hover' : 'dive';
        var samples = Array.isArray(capture.samples) ? capture.samples : [];
        if (!samples.length) {
            return { ok: false, error: 'Nothing was recorded.' };
        }

        ensureDraftForExternalInsert();
        if (!draft) {
            return { ok: false, error: 'Timeline editor is unavailable.' };
        }

        var baseTime = currentTime();
        pushUndo();

        if (presetSelect && presetSelect.querySelector('option[value=""]')) {
            presetSelect.value = '';
        }

        var startEvent = {
            t: baseTime,
            action: (mode === 'dive') ? 'startDive' : 'startHover',
            position: clonePlain(capture.startPosition || { x: 1, y: 0, z: 0 }),
            velocity: clonePlain(capture.startVelocity || { x: 0, y: 0, z: 0 }),
            prevMotionState: !!capture.prevMotionState,
            prevDistance: isFinite(capture.prevDistance) ? capture.prevDistance : undefined,
            observerTime: isFinite(capture.startObserverTime) ? capture.startObserverTime : undefined
        };
        draft.events.push(startEvent);

        var radiusPath = mode + '.currentR';
        var maxT = baseTime;
        var insertedSamples = 0;
        for (var i = 0; i < samples.length; i++) {
            var sample = samples[i];
            if (!sample || typeof sample !== 'object') continue;
            var sampleT = parseFloat(sample.t);
            if (!isFinite(sampleT)) continue;
            var ktime = baseTime + Math.max(0, sampleT);
            var camPos = sample.cameraPosition || {};
            var camQuat = sample.cameraQuaternion || {};
            var sampleRadius = parseFloat(sample.radius);
            if (!isFinite(sampleRadius)) continue;
            var sampleObserverTime = parseFloat(sample.observerTime);
            var panX = isFinite(sample.cameraPanX) ? sample.cameraPanX : 0;
            var panY = isFinite(sample.cameraPanY) ? sample.cameraPanY : 0;
            var posX = isFinite(camPos.x) ? camPos.x : 0;
            var posY = isFinite(camPos.y) ? camPos.y : 0;
            var posZ = isFinite(camPos.z) ? camPos.z : 0;
            var quatX = isFinite(camQuat.x) ? camQuat.x : 0;
            var quatY = isFinite(camQuat.y) ? camQuat.y : 0;
            var quatZ = isFinite(camQuat.z) ? camQuat.z : 0;
            var quatW = isFinite(camQuat.w) ? camQuat.w : 1;

            upsertKey(radiusPath, ktime, sampleRadius, 'linear');
            if (isFinite(sampleObserverTime)) {
                upsertKey('observerState.time', ktime, sampleObserverTime, 'linear');
            }
            upsertKey('cameraPan.x', ktime, panX, 'linear');
            upsertKey('cameraPan.y', ktime, panY, 'linear');
            upsertKey('camera.position.x', ktime, posX, 'linear');
            upsertKey('camera.position.y', ktime, posY, 'linear');
            upsertKey('camera.position.z', ktime, posZ, 'linear');
            upsertKey('camera.quaternion.x', ktime, quatX, 'linear');
            upsertKey('camera.quaternion.y', ktime, quatY, 'linear');
            upsertKey('camera.quaternion.z', ktime, quatZ, 'linear');
            upsertKey('camera.quaternion.w', ktime, quatW, 'linear');
            insertedSamples++;
            if (ktime > maxT) maxT = ktime;
        }

        if (!insertedSamples) {
            redoStack = [];
            draft = undoStack.pop();
            return { ok: false, error: 'Nothing usable was captured.' };
        }

        draft.events.sort(function(a, b) { return a.t - b.t; });
        draft.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });
        draft.duration = Math.max(
            draft.duration,
            maxT,
            baseTime + Math.max(parseFloat(capture.duration) || 0, 0)
        );

        selectedTrack = radiusPath;
        selectedKeyT = baseTime;
        selectedEventIdx = -1;
        clearMultiSelect();
        normalizeQuatSigns();
        applyDraft();
        if (typeof seekPresentation === 'function') seekPresentation(baseTime);
        rebuildAll();
        rememberState();
        setStatus(
            (mode === 'dive' ? 'Dive' : 'Hover') +
            ' capture inserted at t=' + baseTime.toFixed(2) + 's.',
            'tl-status--info'
        );

        return {
            ok: true,
            mode: mode,
            startTime: baseTime,
            endTime: maxT,
            duration: Math.max(0, maxT - baseTime),
            sampleCount: insertedSamples
        };
    }

    // â”€â”€ Path datalist â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

    // â”€â”€ Track list (left column) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function rebuildTrackList() {
        if (!draft || !draft.tracks.length) {
            trackListEl.innerHTML = '<div class="tl-empty">No tracks. Load a preset or add a track.</div>';
            return;
        }
        var html = '';
        for (var i = 0; i < draft.tracks.length; i++) {
            var tr = draft.tracks[i];
            var sel = (tr.path === selectedTrack);
            // Show red border when the row is selected but no individual key is picked (delete-row mode)
            var rowDel = sel && selectedKeys.length === 0;
            var hudActive = (typeof isParamInHud === 'function') ? isParamInHud(tr.path) : false;
            html += '<button type="button" class="tl-track-item' + (sel ? ' is-sel' : '') + (rowDel ? ' is-sel--row' : '') +
                '" data-path="' + esc(tr.path) + '">' +
                '<span class="tl-track-main">' +
                '<span class="tl-track-path">' + esc(tr.path) + '</span>' +
                '<span class="tl-track-hud-btn' + (hudActive ? ' is-active' : '') + '" data-hud-path="' + esc(tr.path) + '" title="' + (hudActive ? 'Hide value from' : 'Show value in') + ' the recording HUD">HUD</span>' +
                '</span>' +
                '<span class="tl-track-keycount">' + tr.keys.length + ' key' + (tr.keys.length === 1 ? '' : 's') + '</span>' +
                '</button>';
        }
        trackListEl.innerHTML = html;
    }

    trackListEl.addEventListener('click', function(e) {
        // HUD toggle â€” must be tested before the track-select logic
        var hudBtn = e.target.closest('[data-hud-path]');
        if (hudBtn) {
            e.stopPropagation();
            var path = hudBtn.getAttribute('data-hud-path');
            if (typeof toggleParamInHud === 'function') toggleParamInHud(path, path);
            rebuildTrackList();
            syncRecModal();
            rememberState();
            return;
        }
        var btn = e.target.closest('[data-path]');
        if (!btn) return;
        selectedTrack = btn.getAttribute('data-path');
        selectedKeyT = NaN;
        clearMultiSelect();
        selectedEventIdx = -1;
        rebuildTrackList();
        rebuildLanes();
        rebuildEventsLane();
        updateInspector();
    });

    // â”€â”€ Ruler (time header in dopesheet) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function updateRuler() {
        var d = getDuration();
        // Compute tick spacing based on the visible time range for nice round numbers
        var visibleD = d / Math.max(1, tlZoom);
        var rawSpacing = visibleD / 10;
        var niceSpacings = [0.1, 0.2, 0.5, 1, 2, 5, 10, 15, 20, 30, 60, 120];
        var spacing = niceSpacings[niceSpacings.length - 1];
        for (var ni = 0; ni < niceSpacings.length; ni++) {
            if (niceSpacings[ni] >= rawSpacing) { spacing = niceSpacings[ni]; break; }
        }
        var html = '';
        var t = 0, safeD = Math.max(d, 0.001);
        while (t <= d + spacing * 0.01) {
            var pct = (t / safeD) * 100;
            var label = t < 60 ? t.toFixed(t % 1 === 0 ? 0 : 1) + 's' : (Math.floor(t / 60) + ':' + ('0' + Math.round(t % 60)).slice(-2));
            html += '<span class="tl-ruler-tick" style="left:' + pct.toFixed(3) + '%">' +
                '<span class="tl-ruler-label">' + label + '</span></span>';
            t = Math.round((t + spacing) * 10000) / 10000;
        }
        var phPct = clamp(currentTime() / safeD * 100, 0, 100);
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

    // â”€â”€ Multi-select helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

    // â”€â”€ Lanes (center dopesheet) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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
            selectedEventIdx = -1;
            showKeyInspector();
            var path = diamond.getAttribute('data-path');
            var ki = parseInt(diamond.getAttribute('data-ki'), 10);
            var track = getTrackByPath(path);
            if (track && track.keys[ki]) {
                var t = track.keys[ki].t;
                var isCtrl = e.ctrlKey || e.metaKey;
                var isShift = e.shiftKey;
                if (isShift) {
                    // Shift+click = select all keys in this column (same time on all tracks)
                    selectedTrack = path;
                    selectedKeyT = t;
                    selectAllKeysAtTime(t);
                } else if (isCtrl) {
                    // Toggle this key in multi-selection
                    if (isKeyMultiSelected(path, t)) {
                        removeFromMultiSelect(path, t);
                    } else {
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
                } else {
                    // Replace selection with just this key
                    clearMultiSelect();
                    addToMultiSelect(path, t);
                    selectedTrack = path;
                    selectedKeyT = t;
                    rebuildTrackList();
                    rebuildLanes();
                    fillInspector(track, ki);
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
            selectedEventIdx = -1;
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

    // â”€â”€ Keyframe drag (move diamonds by dragging) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    lanesEl.addEventListener('pointerdown', function(e) {
        if (e.button !== 0) return;
        var diamond = e.target.closest('.tl-diamond');
        if (!diamond) return;
        var path = diamond.getAttribute('data-path');
        var ki   = parseInt(diamond.getAttribute('data-ki'), 10);
        var track = getTrackByPath(path);
        if (!track || !track.keys[ki]) return;

        var clickedT = track.keys[ki].t;
        // Handle modifier clicks here directly (click event is suppressed
        // by e.preventDefault below, so the click handler won't fire).
        if (e.shiftKey) {
            // Shift+click = select all keys in this column (same time on all tracks)
            selectedEventIdx = -1;
            showKeyInspector();
            selectedTrack = path;
            selectedKeyT = clickedT;
            selectAllKeysAtTime(clickedT);
            return;
        }
        if (e.ctrlKey || e.metaKey) {
            // Toggle this key in multi-selection
            selectedEventIdx = -1;
            showKeyInspector();
            if (isKeyMultiSelected(path, clickedT)) {
                removeFromMultiSelect(path, clickedT);
            } else {
                addToMultiSelect(path, clickedT);
            }
            selectedTrack = path;
            selectedKeyT = clickedT;
            rebuildTrackList();
            rebuildLanes();
            if (selectedKeys.length === 1) {
                fillInspector(track, ki);
            } else {
                inspSummary.textContent = selectedKeys.length + ' keyframes selected.';
            }
            return;
        }
        // Plain click: replace selection with just this key
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

    // â”€â”€ Playhead live-update â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function updatePlayheads() {
        var d = getDuration();
        var phPct = clamp(currentTime() / Math.max(d, 0.001) * 100, 0, 100).toFixed(2) + '%';
        var rulerPH = panel.querySelector('#tl-ruler-playhead');
        if (rulerPH) rulerPH.style.left = phPct;
        var lanePHs = panel.querySelectorAll('.tl-lane-playhead');
        for (var i = 0; i < lanePHs.length; i++) lanePHs[i].style.left = phPct;
    }

    // â”€â”€ Inspector (right column) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function updateInspector() {
        // If an event is selected, keep the event inspector showing
        if (selectedEventIdx >= 0) {
            if (!draft || selectedEventIdx >= draft.events.length) {
                selectedEventIdx = -1;
                showKeyInspector();
            } else {
                fillEventInspector(draft.events[selectedEventIdx], selectedEventIdx);
                return;
            }
        }
        showKeyInspector();
        // Multi-selection summary
        if (selectedKeys.length > 1) {
            var pathSet = {};
            for (var s = 0; s < selectedKeys.length; s++) pathSet[selectedKeys[s].path] = 1;
            var pathCount = Object.keys(pathSet).length;
            inspSummary.textContent = selectedKeys.length + ' keyframes selected (' + pathCount + ' track' + (pathCount > 1 ? 's' : '') + ')';
            inspDel.textContent = 'DELETE KEY';
            inspDel.title = '';
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
        inspSummary.textContent = 'Track: ' + track.path + ' (' + track.keys.length + ' key' + (track.keys.length === 1 ? '' : 's') + ')  â€” Del to remove';
        inspPath.value = track.path;
        inspTime.value = currentTime().toFixed(2);
        if (document.activeElement !== inspValue || !inspValue.value) {
            var liveVal = captureLivePathValue(track.path);
            inspValue.value = (typeof liveVal !== 'undefined') ? formatValue(liveVal) : '';
        }
        inspDel.textContent = 'DELETE TRACK';
        inspDel.title = 'Delete the entire track (all keys)';
    }

    function fillInspector(track, ki) {
        inspDel.textContent = 'DELETE KEY';
        inspDel.title = '';
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

    function captureLivePathValue(path) {
        if (!path || typeof getPresentationPathValue !== 'function') return undefined;
        return getPresentationPathValue(path);
    }

    function captureInspectorValue(path) {
        var v = captureLivePathValue(path);
        if (typeof v === 'undefined') return false;
        inspValue.value = formatValue(v);
        return true;
    }

    function resolveInspectorKeyValue(path) {
        var raw = inspValue.value;
        if (typeof raw === 'string' && raw.trim() === '') {
            var live = captureLivePathValue(path);
            if (typeof live !== 'undefined') {
                inspValue.value = formatValue(live);
                return { value: live, fromLive: true, valid: true };
            }
            return { value: undefined, fromLive: false, valid: false };
        }
        return { value: parseValue(raw), fromLive: false, valid: true };
    }

    // Inspector actions
    inspUseTime.addEventListener('click', function() {
        inspTime.value = currentTime().toFixed(2);
    });
    inspCapture.addEventListener('click', function() {
        var path = (inspPath.value || '').trim();
        if (!path) return;
        captureInspectorValue(path);
    });
    function doSetKey() {
        if (!draft) return;
        var path = (inspPath.value || '').trim();
        if (!path) { setStatus('Path is required.', 'tl-status--warn'); return; }
        var t = parseTime(inspTime.value, currentTime());
        var ease = normalizeEase(inspEase.value);
        var resolvedValue = resolveInspectorKeyValue(path);
        if (!resolvedValue.valid) {
            setStatus('No live value found for this path. Type a value or click LIVE VALUE first.', 'tl-status--warn');
            return;
        }
        var value = resolvedValue.value;

        pushUndo();
        var track = getTrackByPath(path);
        if (!track) {
            track = { path: path, compile: false, keys: [] };
            draft.tracks.push(track);
            draft.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });
        }
        // If the time changed, remove the old keyframe so it moves rather than duplicates
        if (isFinite(selectedKeyT) && path === selectedTrack && Math.abs(t - selectedKeyT) > 1e-4) {
            var oldKi = getKeyAt(track, selectedKeyT);
            if (oldKi >= 0) track.keys.splice(oldKi, 1);
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
        setStatus((resolvedValue.fromLive ? 'Live value keyed: ' : 'Key set: ') + path + ' @ ' + t.toFixed(2) + 's', '');
    }
    inspSet.addEventListener('click', doSetKey);
    inspTime.addEventListener('keydown', function(e) {
        if (e.key === 'Enter') { e.preventDefault(); doSetKey(); }
    });
    inspDel.addEventListener('click', function() {
        if (!draft) return;
        // Track selected but no individual key â€” delete the whole track
        if (selectedKeys.length === 0 && selectedTrack) {
            deleteTrack(selectedTrack);
            return;
        }
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

    // â”€â”€ Annotation lane pointer handler (select + drag markers, multi-channel) â”€â”€
    (function() {
        var dragEi = -1;
        var dragOrigCh = -1;
        var dragTargetCh = -1;
        var dragLaneEl = null;
        var dragStartX = 0;
        var dragged = false;
        var DRAG_THRESHOLD = 4;

        function clearLaneHighlight() {
            var marked = annotLaneColEl.querySelectorAll('.tl-lane--drag-target');
            for (var i = 0; i < marked.length; i++) marked[i].classList.remove('tl-lane--drag-target');
        }

        annotLaneColEl.addEventListener('pointerdown', function(e) {
            if (e.button !== 0) return;
            var marker = e.target.closest('.tl-ev-marker');
            if (marker) {
                dragEi = parseInt(marker.getAttribute('data-ei'), 10);
                dragOrigCh = parseInt(marker.getAttribute('data-ch'), 10);
                dragTargetCh = dragOrigCh;
                dragLaneEl = marker.closest('.tl-events-lane');
                dragStartX = e.clientX;
                dragged = false;
                annotLaneColEl.setPointerCapture(e.pointerId);
                e.preventDefault();
            }
        });

        annotLaneColEl.addEventListener('pointermove', function(e) {
            if (dragEi < 0) return;
            var dx = e.clientX - dragStartX;
            if (!dragged && Math.abs(dx) < DRAG_THRESHOLD) return;
            dragged = true;
            if (!dragLaneEl || !draft) return;
            var rect = dragLaneEl.getBoundingClientRect();
            if (rect.width <= 0) return;
            var dur = getDuration();
            var ratio = clamp((e.clientX - rect.left) / rect.width, 0, 1);
            var newT = ratio * dur;
            draft.events[dragEi].t = Math.max(0, newT);
            var markerEl = dragLaneEl.querySelector('.tl-ev-marker[data-ei="' + dragEi + '"]');
            if (markerEl && dur > 0) markerEl.style.left = ((newT / dur) * 100) + '%';

            // Detect which annotation lane the pointer is hovering over for cross-track drop.
            var overEl = document.elementFromPoint(e.clientX, e.clientY);
            var overLane = overEl && overEl.closest('.tl-events-lane');
            var newTargetCh = dragOrigCh;
            if (overLane) {
                var lch = parseInt(overLane.getAttribute('data-ch'), 10);
                if (isFinite(lch)) newTargetCh = lch;
            }
            if (newTargetCh !== dragTargetCh) {
                clearLaneHighlight();
                dragTargetCh = newTargetCh;
                if (newTargetCh !== dragOrigCh) {
                    var targetLaneEl = annotLaneColEl.querySelector('.tl-events-lane[data-ch="' + newTargetCh + '"]');
                    if (targetLaneEl) targetLaneEl.classList.add('tl-lane--drag-target');
                }
            }
        });

        annotLaneColEl.addEventListener('pointerup', function(e) {
            if (dragEi < 0) return;
            var ei = dragEi;
            var origCh = dragOrigCh;
            var targetCh = dragTargetCh;
            dragEi = -1; dragOrigCh = -1; dragTargetCh = -1;
            clearLaneHighlight();
            annotLaneColEl.releasePointerCapture(e.pointerId);

            if (dragged && draft) {
                pushUndo();
                var ev = draft.events[ei];
                // Reassign channel when dropped onto a different track.
                if (targetCh >= 0 && targetCh !== origCh) {
                    ev.channel = targetCh;
                    // Also move any paired clearAnnotation to the new channel.
                    for (var i = 0; i < draft.events.length; i++) {
                        if (draft.events[i].action === 'clearAnnotation' && draft.events[i]._pairOf === ei) {
                            draft.events[i].channel = targetCh;
                        }
                    }
                }
                draft.events.sort(function(a, b) { return a.t - b.t; });
                for (var i = 0; i < draft.events.length; i++) {
                    if (draft.events[i] === ev) { selectedEventIdx = i; break; }
                }
                // Re-index all _pairOf references invalidated by the sort.
                reindexAllPairOf(draft.events);
                dragLaneEl = null;
                applyDraft();
                rebuildAll();
                fillEventInspector(draft.events[selectedEventIdx], selectedEventIdx);
                var msg = 'Event moved to ' + ev.t.toFixed(2) + 's.';
                if (targetCh >= 0 && targetCh !== origCh) msg += ' Track: ' + (targetCh + 1) + '.';
                setStatus(msg, '');
            } else {
                dragLaneEl = null;
                if (draft && draft.events[ei]) {
                    clearMultiSelect();
                    selectedTrack = '';
                    var clickedEv = draft.events[ei];
                    // If the user clicked a clearAnnotation end-marker, redirect
                    // focus to the paired annotation event instead.
                    if (clickedEv.action === 'clearAnnotation') {
                        var pairIdx = -1;
                        if (typeof clickedEv._pairOf === 'number' && draft.events[clickedEv._pairOf]
                                && draft.events[clickedEv._pairOf].action === 'annotation') {
                            pairIdx = clickedEv._pairOf;
                        } else {
                            // Fallback: nearest preceding annotation on same channel
                            var ch2 = clickedEv.channel || 0;
                            for (var pi = ei - 1; pi >= 0; pi--) {
                                if (draft.events[pi].action === 'annotation'
                                        && (draft.events[pi].channel || 0) === ch2) {
                                    pairIdx = pi; break;
                                }
                            }
                        }
                        if (pairIdx >= 0) {
                            fillEventInspector(draft.events[pairIdx], pairIdx);
                        } else {
                            fillEventInspector(clickedEv, ei);
                        }
                    } else {
                        fillEventInspector(clickedEv, ei);
                    }
                    rebuildEventsLane();
                    rebuildTrackList();
                    rebuildLanes();
                }
            }
        });

        annotLaneColEl.addEventListener('pointercancel', function() {
            if (dragEi < 0) return;
            dragEi = -1; dragOrigCh = -1; dragTargetCh = -1;
            clearLaneHighlight();
            dragLaneEl = null;
            rebuildEventsLane();
        });

        // Background click: deselect + seek, and track active channel
        annotLaneColEl.addEventListener('click', function(e) {
            if (e.target.closest('.tl-ev-marker')) return;
            var lane = e.target.closest('.tl-events-lane');
            if (lane) {
                var lch = parseInt(lane.getAttribute('data-ch'), 10);
                if (isFinite(lch)) selectedEventChannel = lch;
            }
            if (selectedEventIdx >= 0) {
                selectedEventIdx = -1;
                showKeyInspector();
                rebuildEventsLane();
            }
            if (lane) {
                var rect = lane.getBoundingClientRect();
                if (rect.width > 0) {
                    var ratio = clamp((e.clientX - rect.left) / rect.width, 0, 1);
                    var t2 = ratio * getDuration();
                    if (typeof seekPresentation === 'function') seekPresentation(t2);
                    updateTimeInputs(); updateScrubber(); updatePlayheads();
                }
            }
        });
    })();

    // â”€â”€ Text event inspector handlers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    evUseTimeBtn.addEventListener('click', function() {
        evTimeInput.value = currentTime().toFixed(2);
    });
    evSetBtn.addEventListener('click', doSetEvent);
    evTimeInput.addEventListener('keydown', function(e) {
        if (e.key === 'Enter') { e.preventDefault(); doSetEvent(); }
    });
    // Auto-sync inspector fields to draft while typing
    evTitleInput.addEventListener('input', syncInspectorToDraft);
    evBodyInput.addEventListener('input', syncInspectorToDraft);
    evColorInput.addEventListener('input', syncInspectorToDraft);
    evWidthInput.addEventListener('change', syncInspectorToDraft);
    evFadeInput.addEventListener('change', syncInspectorToDraft);
    // Duration field: immediately update (or remove) the paired clearAnnotation event
    evDurInput.addEventListener('change', function() {
        if (!draft || selectedEventIdx < 0 || !draft.events[selectedEventIdx]) return;
        var ev = draft.events[selectedEventIdx];
        if (ev.action !== 'annotation') return;
        var durVal = parseFloat(evDurInput.value);
        if (!isFinite(durVal) || durVal < 0) durVal = 0;
        pushUndo();
        // Remove any existing paired clearAnnotation for this event
        for (var i = draft.events.length - 1; i >= 0; i--) {
            if (draft.events[i].action === 'clearAnnotation' &&
                    draft.events[i]._pairOf === selectedEventIdx) {
                draft.events.splice(i, 1);
                if (i < selectedEventIdx) selectedEventIdx--;
            }
        }
        // Create a new one if duration > 0
        if (durVal > 0) {
            var newClear = { t: ev.t + durVal, action: 'clearAnnotation', channel: ev.channel || 0 };
            draft.events.push(newClear);
            draft.events.sort(function(a, b) { return a.t - b.t; });
            for (var j = 0; j < draft.events.length; j++) {
                if (draft.events[j] === ev) { selectedEventIdx = j; break; }
            }
            newClear._pairOf = selectedEventIdx;
        }
        reindexAllPairOf(draft.events);
        applyDraft();
        rebuildAll();
    });
    evDelBtn.addEventListener('click', function() {
        if (!draft || selectedEventIdx < 0) return;
        pushUndo();
        draft.events.splice(selectedEventIdx, 1);
        reindexAllPairOf(draft.events);
        selectedEventIdx = -1;
        showKeyInspector();
        applyDraft();
        rebuildAll();
        setStatus('Annotation event deleted.', '');
    });
    // Remove End button: deletes only the clearAnnotation end marker so the
    // annotation that owns it becomes permanent (no automatic clear).
    evEndRemove.addEventListener('click', function() {
        if (!draft || selectedEventIdx < 0) return;
        var ev = draft.events[selectedEventIdx];
        if (!ev || ev.action !== 'clearAnnotation') return;
        // Find the paired annotation to keep selected after deletion
        var pairIdx = -1;
        if (typeof ev._pairOf === 'number' && draft.events[ev._pairOf]
                && draft.events[ev._pairOf].action === 'annotation') {
            pairIdx = ev._pairOf;
        }
        pushUndo();
        draft.events.splice(selectedEventIdx, 1);
        reindexAllPairOf(draft.events);
        // Select the paired annotation (index may have shifted after splice)
        if (pairIdx >= 0 && pairIdx < selectedEventIdx) {
            selectedEventIdx = pairIdx;
        } else if (pairIdx > selectedEventIdx) {
            selectedEventIdx = pairIdx - 1;
        } else {
            selectedEventIdx = -1;
        }
        if (selectedEventIdx >= 0 && draft.events[selectedEventIdx]) {
            fillEventInspector(draft.events[selectedEventIdx], selectedEventIdx);
        } else {
            selectedEventIdx = -1;
            showKeyInspector();
        }
        applyDraft();
        rebuildAll();
        setStatus('End marker removed. Annotation will now persist.', '');
    });
    evNewBtn.addEventListener('click', addNewAnnotationEvent);
    addTextBtn.addEventListener('click', addNewAnnotationEvent);
    evPositionBtn.addEventListener('click', startPositionDrag);
    evPreviewBtn.addEventListener('click', previewAnnotation);

    // â”€â”€ Add track â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

    // â”€â”€ Auto Keyframe â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    // Parameters that are UI/performance meta-settings and should never be
    // captured as timeline keyframes (they would override the user's choice on play).
    // These are all the parameters owned by the quality preset system.
    var SNAPSHOT_EXCLUDED_PATHS = {
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

    function captureSnapshot() {
        var out = {};
        function add(p, v) {
            if (SNAPSHOT_EXCLUDED_PATHS[p]) return;
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
    // â”€â”€ Auto Key initial-state modal â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function showAutoKeyModal(t, snap, paths, onFull, onDiff) {
        var overlay = document.createElement('div');
        overlay.className = 'tl-modal-overlay';
        var tLabel = t < 0.01 ? 't = 0' : 't = ' + t.toFixed(2) + 's';
        overlay.innerHTML =
            '<div class="tl-modal">' +
            '<div class="tl-modal-title">Auto Key &#8212; Initial State</div>' +
            '<p class="tl-modal-body">First Auto Key capture at <b>' + esc(tLabel) + '</b>.<br>' +
            'How should the initial state be recorded?</p>' +
            '<div class="tl-modal-choices">' +
            '<button class="tl-modal-btn tl-modal-btn--full" type="button">' +
            '<span class="tl-modal-btn-label">Full initial state</span>' +
            '<span class="tl-modal-btn-sub">Key every parameter now (' + paths.length + ' values) &mdash; guarantees a clean reset on play</span>' +
            '</button>' +
            '<button class="tl-modal-btn tl-modal-btn--diff" type="button">' +
            '<span class="tl-modal-btn-label">Changes only (diff mode)</span>' +
            '<span class="tl-modal-btn-sub">Only record what changes on the next press &mdash; lighter, manual approach</span>' +
            '</button>' +
            '</div>' +
            '<button class="tl-modal-btn tl-modal-btn--cancel" type="button">Cancel</button>' +
            '</div>';
        document.body.appendChild(overlay);

        function close() { if (overlay.parentNode) overlay.parentNode.removeChild(overlay); }
        overlay.querySelector('.tl-modal-btn--full').addEventListener('click', function() { close(); onFull(); });
        overlay.querySelector('.tl-modal-btn--diff').addEventListener('click', function() { close(); onDiff(); });
        overlay.querySelector('.tl-modal-btn--cancel').addEventListener('click', close);
        overlay.addEventListener('mousedown', function(e) { if (e.target === overlay) close(); });
        overlay.addEventListener('keydown', function(e) { if (e.key === 'Escape') close(); });
        // Focus the overlay so Escape works immediately
        overlay.setAttribute('tabindex', '-1');
        overlay.focus();
    }

    autoKeyBtn.addEventListener('click', function(e) {
        if (!draft) { setStatus('Load a timeline first.', 'tl-status--warn'); return; }
        var skipCamera = e.shiftKey;
        var t = currentTime();
        var snap = captureSnapshot();
        var paths = Object.keys(snap);
        if (skipCamera) {
            paths = paths.filter(function(p) {
                return p.indexOf('camera.position.') !== 0 &&
                       p.indexOf('camera.quaternion.') !== 0;
            });
        }
        if (!paths.length) { setStatus('Nothing to capture.', 'tl-status--warn'); return; }

        if (!autoKeySnapshot) {
            // First press â€” check if this looks like an initial state (tâ‰ˆ0 or no tracks yet)
            var isInitial = (t < 0.01 || !draft.tracks.length);
            if (isInitial) {
                showAutoKeyModal(t, snap, paths,
                    // Full state
                    function() {
                        pushUndo();
                        var count = 0;
                        for (var i = 0; i < paths.length; i++) {
                            upsertKey(paths[i], 0, snap[paths[i]], 'linear');
                            count++;
                        }
                        draft.tracks.sort(function(a, b) { return a.path.localeCompare(b.path); });
                        normalizeQuatSigns();
                        autoKeySnapshot = { time: t, values: clonePlain(snap) };
                        applyDraft();
                        rebuildAll();
                        setStatus('Full initial state saved: ' + count + ' parameters keyed at t=0. Now change controls & press AUTO KEY.', 'tl-status--info');
                    },
                    // Diff mode
                    function() {
                        autoKeySnapshot = { time: t, values: clonePlain(snap) };
                        setStatus('Baseline captured at ' + t.toFixed(2) + 's. Change controls, press AUTO KEY again.', 'tl-status--info');
                    }
                );
                return;
            }
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
            if (typeof seekPresentation === 'function') seekPresentation(t);
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

    // â”€â”€ Quaternion sign normalisation â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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
                    // Derive a name from the JSON name field or the filename
                    var presetName = (obj && obj.name && obj.name.trim()) ? obj.name.trim()
                        : file.name.replace(/\.json$/i, '').replace(/[_-]+/g, ' ').trim();
                    if (!presetName) presetName = 'Imported';
                    obj.name = presetName;

                    if (typeof registerPresentationPreset === 'function') {
                        registerPresentationPreset(obj, presetName);
                    }
                    importedPresets[presetName] = true;

                    if (typeof setPresentationTimeline === 'function') {
                        applyingDraft = true;
                        if (setPresentationTimeline(obj)) {
                            applyingDraft = false;
                            syncFromRuntime();
                            linkedFileName = file.name;
                            populatePresets();
                            presetSelect.value = presetName;
                            updateDelPresetBtn();
                            setStatus('Imported: ' + file.name, '');
                        } else {
                            applyingDraft = false;
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
        var defaultName = (tl.name && tl.name !== 'Untitled' && tl.name !== 'Custom')
            ? tl.name : '';
        var exportName = prompt('Name for the exported timeline:', defaultName);
        if (exportName === null) return; // cancelled
        exportName = exportName.trim();
        if (!exportName) exportName = 'timeline';
        tl.name = exportName;
        if (draft) draft.name = exportName;
        var json = JSON.stringify(tl, null, 2);
        var blob = new Blob([json], { type: 'application/json' });
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        var safeFilename = exportName.replace(/[^a-zA-Z0-9_-]/g, '_');
        a.download = safeFilename + '.json';
        document.body.appendChild(a);
        a.click();
        a.remove();
        URL.revokeObjectURL(url);
        linkedFileName = safeFilename + '.json';

        // Register in the preset list (same as import) so it can be reloaded
        if (typeof registerPresentationPreset === 'function') {
            registerPresentationPreset(clonePlain(tl), exportName);
        }
        importedPresets[exportName] = true;
        populatePresets();
        presetSelect.value = exportName;
        updateDelPresetBtn();

        setStatus('Exported: ' + a.download, '');
    });

    // â”€â”€ Delete imported preset â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    delPresetBtn.addEventListener('click', function() {
        var name = presetSelect.value;
        if (!name || !importedPresets[name]) return;
        if (!confirm('Remove imported preset "' + name + '" from the list?')) return;
        delete importedPresets[name];
        if (typeof PRESENTATION_PRESETS !== 'undefined') delete PRESENTATION_PRESETS[name];
        if (typeof PRESENTATION_PRESET_ORDER !== 'undefined') {
            var idx = PRESENTATION_PRESET_ORDER.indexOf(name);
            if (idx !== -1) PRESENTATION_PRESET_ORDER.splice(idx, 1);
        }
        presetSelect.value = '';
        loadPresetByName('');
        populatePresets();
        setStatus('Preset "' + name + '" removed.', '');
    });

    // â”€â”€ Save button â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    saveBtn.addEventListener('click', function() {
        var tl = typeof getPresentationTimeline === 'function' ? getPresentationTimeline() : null;
        if (!tl && draft) tl = clonePlain(draft);
        if (!tl) { setStatus('No timeline to save.', 'tl-status--warn'); return; }

        if (!linkedFileName) {
            // Not linked to any fileâ€”ask for a name and export as a new file
            var saveName = prompt('No file linked. Enter a name to save as:', (tl.name && tl.name !== 'Untitled' && tl.name !== 'Custom') ? tl.name : '');
            if (saveName === null) return; // cancelled
            saveName = saveName.trim();
            if (!saveName) { setStatus('Save cancelled.', 'tl-status--warn'); return; }
            tl.name = saveName;
            if (draft) draft.name = saveName;
            linkedFileName = saveName.replace(/[^a-zA-Z0-9_-]/g, '_') + '.json';
        }

        var json = JSON.stringify(tl, null, 2);
        var blob = new Blob([json], { type: 'application/json' });
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = linkedFileName;
        document.body.appendChild(a);
        a.click();
        a.remove();
        URL.revokeObjectURL(url);

        // Also update the registered preset if it exists
        var name = tl.name || (draft && draft.name) || '';
        if (name && typeof registerPresentationPreset === 'function') {
            registerPresentationPreset(tl, name);
            importedPresets[name] = true;
            populatePresets();
            presetSelect.value = name;
            updateDelPresetBtn();
        }
        setStatus('Saved: ' + a.download, '');
    });

    // â”€â”€ Listen for external timeline loads â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    window.addEventListener('presentation:timeline-panel-sync', function() {
        if (panelOpen && !applyingDraft) syncFromRuntime();
    });

    // â•â• Motion Functions â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
    // Each entry: params[] with { id, label, type, min, step, def, defaultFn, options }
    var MOTION_TYPES = {
        sky_reveal: {
            params: [
                { id: 'start',    label: 'Start time (s)',  type: 'number', min: 0,   step: 0.1, defaultFn: function() { return currentTime().toFixed(2); } },
                { id: 'duration', label: 'Duration (s)',    type: 'number', min: 0.1, step: 1,   def: 12 },
                { id: 'from',     label: 'BH enters from', type: 'select', options: [['bottom','Bottom \u2014 see sky above'],['top','Top \u2014 see sky below'],['left','Left \u2014 see sky to the right'],['right','Right \u2014 see sky to the left']] },
                { id: 'ease',     label: 'Ease',            type: 'select', options: [['smoother','smoother'],['smooth','smooth'],['linear','linear']] }
            ]
        },
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

            // Seed camera position â€” prefer draft track value at start time, then runtime
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
            // Use 32 steps per orbit â€” all linear so velocity is perfectly constant.
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
                // Always linear â€” smooth/smoother ease would decelerate at every
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

        } else if (type === 'sky_reveal') {
            var fromDir = getMotionParam('from') || 'bottom';

            // cam_pan shifts the rendered image so the BH centre appears at
            // screen coords  screen_pos = -cam_pan.
            // Screen x \u2208 [-1, 1], screen y \u2248 [-0.56, 0.56] (1920\u00d71080 aspect).
            // Starting offsets large enough to push the BH fully outside any frame.
            var startPanX = 0, startPanY = 0;
            if      (fromDir === 'bottom') { startPanY = +1.0; }  // BH below \u2192 starfield above
            else if (fromDir === 'top')    { startPanY = -1.0; }  // BH above \u2192 equatorial sky below
            else if (fromDir === 'left')   { startPanX = +2.0; }  // BH to the left \u2192 sky to the right
            else if (fromDir === 'right')  { startPanX = -2.0; }  // BH to the right \u2192 sky to the left

            // End value: current cameraPan (or 0)
            var endPanX = 0, endPanY = 0;
            if (typeof getPresentationPathValue === 'function') {
                var epx = getPresentationPathValue('cameraPan.x');
                var epy = getPresentationPathValue('cameraPan.y');
                if (typeof epx === 'number') endPanX = epx;
                if (typeof epy === 'number') endPanY = epy;
            }

            upsertKey('cameraPan.x', start,            endPanX + startPanX, 'linear');
            upsertKey('cameraPan.x', start + duration, endPanX,             ease);
            upsertKey('cameraPan.y', start,            endPanY + startPanY, 'linear');
            upsertKey('cameraPan.y', start + duration, endPanY,             ease);
            selectedTrack = 'cameraPan.x';
            setStatus('Sky reveal: BH enters from ' + fromDir + ' over ' + duration + 's at t=' + start.toFixed(2) + 's.', '');
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

    // â”€â”€ Recording modal â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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
            { v: '3840x2160', l: '3840\xd72160 (4K)' },
            { v: '7680x4320', l: '7680\xd74320 (8K)' }
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

    function renderRecHudList() {
        if (!recHudListEl || !recHudEmptyEl) return;
        var hs = (typeof getPresentationParamHudState === 'function')
            ? getPresentationParamHudState()
            : { items: [] };
        var items = Array.isArray(hs.items) ? hs.items : [];
        if (recClearHudBtn) recClearHudBtn.disabled = items.length === 0;
        if (recAddSelectedBtn) {
            recAddSelectedBtn.disabled = !selectedTrack ||
                (typeof isParamInHud === 'function' && isParamInHud(selectedTrack));
        }
        if (!items.length) {
            recHudEmptyEl.style.display = '';
            recHudListEl.innerHTML = '';
            return;
        }
        recHudEmptyEl.style.display = 'none';
        var html = '';
        for (var i = 0; i < items.length; i++) {
            var item = items[i];
            var label = item.label || item.path || '';
            html += '<button type="button" class="tl-rec-chip" data-remove-hud-path="' + esc(item.path) + '"' +
                ' title="Remove ' + esc(item.path) + ' from the HUD">' +
                '<span class="tl-rec-chip-label">' + esc(label) + '</span>' +
                '<span class="tl-rec-chip-x">&times;</span>' +
                '</button>';
        }
        recHudListEl.innerHTML = html;
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
        if (recParamHudShowCb)   recParamHudShowCb.checked   = !!s.param_hud_enabled;
        if (recParamHudRecordCb) recParamHudRecordCb.checked = !!s.param_hud_in_recording;
        // Sync HUD position/size inputs and show/hide layout row
        if (hudLayoutRow) hudLayoutRow.style.display = (s.param_hud_count > 0) ? '' : 'none';
        if (typeof getPresentationParamHudState === 'function') {
            var hs = getPresentationParamHudState();
            if (hudXInput  && document.activeElement !== hudXInput)  hudXInput.value  = hs.anchorX.toFixed(2);
            if (hudYInput  && document.activeElement !== hudYInput)  hudYInput.value  = hs.anchorY.toFixed(2);
            if (hudFsInput && document.activeElement !== hudFsInput) hudFsInput.value = hs.fontSize;
        }
        renderRecHudList();

        // Refresh resolution label with live viewport size
        refreshRecResolutionLabel();

        // Keep REC button indicator in sync
        if (recBtn) recBtn.classList.toggle('is-recording', !!s.recording);
        if (recShotBtn) recShotBtn.disabled = !!s.recording;

        if (!s.recording) {
            recStartBtn.disabled = false;
            recStopBtn.disabled  = true;
            if (s.recording_offline_unavailable_reason) {
                setRecStatus(s.recording_offline_unavailable_reason, 'is-warning');
            } else {
                setRecStatus('Idle', '');
            }
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
        } else if (phase === 'finalizing-encode' || phase === 'finalizing-mux' || phase === 'finalizing-download') {
            var fp = s.recording_offline_finalizing_progress;
            var label = phase === 'finalizing-encode' ? 'Encoding' : phase === 'finalizing-mux' ? 'Muxing' : 'Downloading';
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
        applySavedRecordingUiState();
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
        rememberState();
    });
    recAnnotCb.addEventListener('change', function() {
        if (typeof setPresentationAnnotationsEnabled === 'function') setPresentationAnnotationsEnabled(recAnnotCb.checked);
        rememberState();
    });
    recAnnotRecordCb.addEventListener('change', function() {
        if (typeof setPresentationAnnotationsIncludedInRecording === 'function') setPresentationAnnotationsIncludedInRecording(recAnnotRecordCb.checked);
        rememberState();
    });
    recParamHudShowCb.addEventListener('change', function() {
        if (typeof setPresentationParamHudEnabled === 'function') setPresentationParamHudEnabled(recParamHudShowCb.checked);
        rememberState();
    });
    recParamHudRecordCb.addEventListener('change', function() {
        if (typeof setPresentationParamHudIncludedInRecording === 'function') setPresentationParamHudIncludedInRecording(recParamHudRecordCb.checked);
        rememberState();
    });

    if (recHudListEl) {
        recHudListEl.addEventListener('click', function(e) {
            var chip = e.target.closest('[data-remove-hud-path]');
            if (!chip) return;
            var path = chip.getAttribute('data-remove-hud-path');
            if (path && typeof removeParamFromHud === 'function') removeParamFromHud(path);
            syncRecModal();
            rebuildTrackList();
            rememberState();
        });
    }
    if (recAddSelectedBtn) {
        recAddSelectedBtn.addEventListener('click', function() {
            if (!selectedTrack || typeof addParamToHud !== 'function') return;
            addParamToHud(selectedTrack, selectedTrack);
            syncRecModal();
            rebuildTrackList();
            rememberState();
        });
    }
    if (recClearHudBtn) {
        recClearHudBtn.addEventListener('click', function() {
            if (typeof clearParamHud !== 'function') return;
            clearParamHud();
            syncRecModal();
            rebuildTrackList();
            rememberState();
        });
    }

    function applyHudLayout() {
        if (typeof setParamHudLayout !== 'function') return;
        var x  = parseFloat(hudXInput.value);
        var y  = parseFloat(hudYInput.value);
        var fs = parseFloat(hudFsInput.value);
        setParamHudLayout({
            anchorX:  isFinite(x)  ? x  : undefined,
            anchorY:  isFinite(y)  ? y  : undefined,
            fontSize: isFinite(fs) ? fs : undefined
        });
        syncRecModal();
        rememberState();
    }
    if (hudXInput)    hudXInput.addEventListener('input', applyHudLayout);
    if (hudYInput)    hudYInput.addEventListener('input', applyHudLayout);
    if (hudFsInput)   hudFsInput.addEventListener('input', applyHudLayout);
    if (recQualitySelect) recQualitySelect.addEventListener('change', rememberState);
    if (recModeSelect)    recModeSelect.addEventListener('change', rememberState);
    if (recResSelect)     recResSelect.addEventListener('change', rememberState);
    if (recFpsInput)      recFpsInput.addEventListener('input', rememberState);
    if (recBitrateInput)  recBitrateInput.addEventListener('input', rememberState);
    if (recResetSimCb)    recResetSimCb.addEventListener('change', rememberState);

    recShotBtn.addEventListener('click', function() {
        if (typeof capturePresentationScreenshot !== 'function') return;
        if (recQualitySelect && recQualitySelect.querySelector('option[value="cinematic"]')) {
            recQualitySelect.value = 'cinematic';
        }
        var captured = capturePresentationScreenshot({
            qualityPreset: 'cinematic',
            recordingResolution: recResSelect ? recResSelect.value : 'current',
            includeAnnotationsInScreenshot: !!recAnnotRecordCb.checked
        });
        if (!captured) {
            var st = (typeof getPresentationState === 'function') ? getPresentationState() : {};
            setRecStatus(st.recording_offline_unavailable_reason || 'Failed to capture screenshot.', 'is-warning');
        } else {
            setRecStatus('Offline screenshot downloaded (.png).', '');
        }
        syncRecModal();
    });

    recStartBtn.addEventListener('click', async function() {
        if (typeof startPresentationRecording !== 'function') return;
        if (recResetSimCb && recResetSimCb.checked) {
            if (typeof observer !== 'undefined' && observer) observer.time = 0.0;
            if (typeof shader !== 'undefined' && shader) shader.needsUpdate = true;
        }
        var fps     = clampNum(parseFloat(recFpsInput.value) || 60, 24, 120);
        var bitrate = clampNum(parseFloat(recBitrateInput.value) || 20, 4, 80);
        recFpsInput.value     = Math.round(fps);
        recBitrateInput.value = Math.round(bitrate);
        var recMode = recModeSelect ? recModeSelect.value : 'offline';
        var recOptions = {
            fps: fps,
            bitrateMbps: bitrate,
            autoStopOnPresentationEnd: true,
            recordingMode:       recMode,
            recordingResolution: recResSelect     ? recResSelect.value     : 'current',
            qualityPreset:       recQualitySelect ? recQualitySelect.value : 'optimal',
            includeAnnotationsInRecording: !!recAnnotRecordCb.checked
        };

        // For offline mode, try to use the File System Access API to stream encoded
        // chunks directly to disk. This avoids accumulating the entire video as a
        // growing ArrayBuffer in RAM, which causes GPU/browser crashes on long renders.
        if (recMode === 'offline' &&
            typeof showSaveFilePicker === 'function' &&
            typeof window.WebMMuxer !== 'undefined' &&
            typeof window.WebMMuxer.FileSystemWritableFileStreamTarget === 'function') {
            try {
                var suggestedName = (typeof presentationCaptureFilename === 'function')
                    ? presentationCaptureFilename('black-hole-presentation', 'video/webm')
                    : 'black-hole-presentation.webm';
                var fileHandle = await showSaveFilePicker({
                    suggestedName: suggestedName,
                    types: [{ description: 'WebM video', accept: { 'video/webm': ['.webm'] } }]
                });
                recOptions.writableFileStream = await fileHandle.createWritable();
            } catch (pickerErr) {
                if (pickerErr.name === 'AbortError') return; // user cancelled the dialog
                // File picker failed for another reason â€” fall back to ArrayBuffer mode silently
                console.warn('File picker failed, falling back to in-memory recording:', pickerErr);
            }
        }

        var started = startPresentationRecording(recOptions);
        if (!started) {
            if (recOptions.writableFileStream) {
                try { recOptions.writableFileStream.abort(); } catch (e) {}
            }
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

    // â”€â”€ Keyboard shortcuts â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    function deleteTrack(path) {
        if (!draft || !path) return;
        pushUndo();
        draft.tracks = draft.tracks.filter(function(tr) { return tr.path !== path; });
        // Remove any selectedKeys that were on this track
        selectedKeys = selectedKeys.filter(function(sk) { return sk.path !== path; });
        if (selectedTrack === path) {
            selectedTrack = draft.tracks.length ? draft.tracks[0].path : '';
            selectedKeyT = NaN;
        }
        applyDraft();
        rebuildAll();
        setStatus('Track deleted: ' + path, '');
    }

    function deleteSelectedKeys() {
        if (!draft) return;
        // If no keys are selected but a track is active, delete the whole track
        if (!selectedKeys.length && selectedTrack) {
            deleteTrack(selectedTrack);
            return;
        }
        if (!selectedKeys.length) return;
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

    // â”€â”€ Copy / Paste keyframes â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

        // Ctrl+Z / Ctrl+Y â€” fire whenever the panel is open, unless focus is
        // in an input that lives *outside* the timeline panel itself.
        var inExternalInput = inInput && !panel.contains(e.target);
        if ((e.ctrlKey || e.metaKey) && !e.shiftKey && e.code === 'KeyZ' && !inExternalInput) {
            e.preventDefault();
            e.stopPropagation();
            undo();
            return;
        }
        if ((e.ctrlKey || e.metaKey) && (e.code === 'KeyY' || (e.shiftKey && e.code === 'KeyZ')) && !inExternalInput) {
            e.preventDefault();
            e.stopPropagation();
            redo();
            return;
        }

        if (inInput) return;

        // Escape â†’ close motion modal if open, or deselect event
        if (e.code === 'Escape') {
            if (motionModal.classList.contains('is-open')) { closeMotionModal(); e.preventDefault(); return; }
            if (selectedEventIdx >= 0) {
                selectedEventIdx = -1;
                showKeyInspector();
                rebuildEventsLane();
                e.preventDefault();
            }
            return;
        }

        // Space â†’ play / pause
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

        // K â†’ key the selected track at the playhead using the current live value
        if (!e.ctrlKey && !e.metaKey && !e.altKey && e.code === 'KeyK') {
            var keyPath = selectedTrack || (inspPath.value || '').trim();
            if (!keyPath) {
                setStatus('Select a track first, then press K to key its live value.', 'tl-status--warn');
                e.preventDefault();
                return;
            }
            inspPath.value = keyPath;
            inspTime.value = currentTime().toFixed(2);
            captureInspectorValue(keyPath);
            doSetKey();
            e.preventDefault();
            return;
        }

        // Delete / Backspace â†’ delete selected keyframes, or the selected track if none selected
        if (e.code === 'Delete' || e.code === 'Backspace') {
            e.preventDefault();
            deleteSelectedKeys();
            return;
        }

        // Shift+A â†’ select all keys in the same column as the current playhead time
        if (e.shiftKey && e.code === 'KeyA') {
            e.preventDefault();
            selectAllKeysAtTime(currentTime());
            return;
        }

        // Ctrl+A â†’ select all keyframes (on current track if one is selected, else all)
        if ((e.ctrlKey || e.metaKey) && e.code === 'KeyA') {
            e.preventDefault();
            if (selectedTrack && !e.shiftKey) selectAllKeysOnTrack();
            else selectAllKeys();
            return;
        }

        // Ctrl+C â†’ copy selected keyframes to clipboard
        if ((e.ctrlKey || e.metaKey) && !e.shiftKey && e.code === 'KeyC') {
            e.preventDefault();
            copySelectedKeys();
            return;
        }

        // Ctrl+V â†’ paste keyframes at current time
        if ((e.ctrlKey || e.metaKey) && !e.shiftKey && e.code === 'KeyV') {
            e.preventDefault();
            pasteKeys();
            return;
        }

        // Home â†’ seek to start
        if (e.code === 'Home') {
            e.preventDefault();
            if (typeof seekPresentation === 'function') seekPresentation(0);
            updateTimeInputs(); updateScrubber(); updatePlayheads();
            return;
        }

        // End â†’ seek to end
        if (e.code === 'End') {
            e.preventDefault();
            if (typeof seekPresentation === 'function') seekPresentation(getDuration());
            updateTimeInputs(); updateScrubber(); updatePlayheads();
            return;
        }

        // Left/Right arrows â†’ nudge time
        if (e.code === 'ArrowLeft' || e.code === 'ArrowRight') {
            var step = e.shiftKey ? 1.0 : 0.1;
            var dir = (e.code === 'ArrowLeft') ? -1 : 1;
            var newT = clamp(currentTime() + dir * step, 0, getDuration());
            if (typeof seekPresentation === 'function') seekPresentation(newT);
            updateTimeInputs(); updateScrubber(); updatePlayheads();
            e.preventDefault();
            return;
        }
    }, true);

    // â”€â”€ Public API â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    timelinePanelBinding = {
        open: function() { setPanelOpen(true); },
        close: function() { setPanelOpen(false); },
        toggle: toggle,
        isOpen: function() { return panelOpen; },
        sync: syncFromRuntime,
        loadPreset: loadPresetByName,
        syncRecState: syncRecModal,
        insertAnimationCapture: insertAnimationCapture
    };
    if (typeof registerBlackHoleUiBinding === 'function') {
        registerBlackHoleUiBinding('timelinePanel', timelinePanelBinding);
    }
    return timelinePanelBinding;
}



