// Role: Presentation controls embedded in the ANIMATIONS panel.
//       Handles timeline preset selection, playback, recording, and
//       annotation toggles without using dat.GUI rows.

function createPresentationAnimationSectionHtml() {
    return '' +
        '<div class="anim-section" id="presentation-section">' +
            '<button class="anim-section-toggle" id="presentation-section-toggle" ' +
                'type="button" aria-expanded="false" aria-controls="presentation-section-body">' +
                '<span class="anim-section-arrow">&#9654;</span> PRESENTATION TIMELINE' +
            '</button>' +
            '<div class="anim-section-body" id="presentation-section-body">' +
                '<div class="presentation-desc">Scripted walkthrough and recorder for slides/videos.</div>' +

                '<div class="presentation-row">' +
                    '<label for="presentation-preset-select">Preset</label>' +
                    '<select id="presentation-preset-select" class="presentation-select"></select>' +
                '</div>' +

                '<div class="presentation-row presentation-toggles">' +
                    '<label><input type="checkbox" id="presentation-loop"> Loop timeline</label>' +
                '</div>' +
                '<div class="presentation-row presentation-toggles">' +
                    '<label><input type="checkbox" id="presentation-annotations-enabled" checked> Show explanations</label>' +
                '</div>' +
                '<div class="presentation-row presentation-toggles">' +
                    '<label><input type="checkbox" id="presentation-annotations-recording"> Include text in recording</label>' +
                '</div>' +
                '<div class="presentation-row presentation-toggles">' +
                    '<label><input type="checkbox" id="presentation-record-reset-sim" checked> Reset simulation on rec start</label>' +
                '</div>' +

                '<div class="presentation-row presentation-timeline-row">' +
                    '<label for="presentation-scrub">Time</label>' +
                    '<div class="presentation-timeline-stack">' +
                        '<div id="presentation-timeline" class="presentation-timeline" aria-hidden="true">' +
                            '<div id="presentation-timeline-fill" class="presentation-timeline-fill"></div>' +
                            '<div id="presentation-timeline-markers" class="presentation-timeline-markers"></div>' +
                            '<div id="presentation-timeline-playhead" class="presentation-timeline-playhead"></div>' +
                        '</div>' +
                        '<input id="presentation-scrub" type="range" min="0" max="1" step="0.01" value="0">' +
                    '</div>' +
                    '<span id="presentation-scrub-label">0.0s</span>' +
                '</div>' +

                '<div class="presentation-btn-row">' +
                    '<button id="presentation-play-btn" class="presentation-btn presentation-btn-primary" type="button">&#9654; PLAY</button>' +
                    '<button id="presentation-pause-btn" class="presentation-btn" type="button">&#10074;&#10074; PAUSE</button>' +
                    '<button id="presentation-stop-btn" class="presentation-btn" type="button">&#9632; STOP</button>' +
                '</div>' +

                (typeof createPresentationEditorHtml === 'function'
                    ? createPresentationEditorHtml()
                    : '') +

                '<div class="presentation-row">' +
                    '<label for="presentation-record-quality">Rec quality</label>' +
                    '<select id="presentation-record-quality" class="presentation-select"></select>' +
                '</div>' +

                '<div class="presentation-row">' +
                    '<label for="presentation-record-mode">Rec mode</label>' +
                    '<select id="presentation-record-mode" class="presentation-select"></select>' +
                '</div>' +

                '<div class="presentation-row">' +
                    '<label for="presentation-record-resolution">Rec resolution</label>' +
                    '<select id="presentation-record-resolution" class="presentation-select"></select>' +
                '</div>' +

                '<div class="presentation-row">' +
                    '<label for="presentation-fps">Rec FPS</label>' +
                    '<input id="presentation-fps" type="number" min="24" max="120" step="1" value="60">' +
                    '<label for="presentation-bitrate" class="presentation-inline-label">Mbps</label>' +
                    '<input id="presentation-bitrate" type="number" min="4" max="80" step="1" value="20">' +
                '</div>' +

                '<div class="presentation-btn-row">' +
                    '<button id="presentation-record-start-btn" class="presentation-btn presentation-btn-record" type="button">&#9679; START REC</button>' +
                    '<button id="presentation-record-stop-btn" class="presentation-btn" type="button">&#9632; STOP REC</button>' +
                '</div>' +

                '<div id="presentation-status" class="presentation-status">Idle</div>' +
            '</div>' +
        '</div>';
}

function bindPresentationAnimationSection(panelRoot) {
    if (!panelRoot) return;

    var section = panelRoot.querySelector('#presentation-section');
    if (!section) return;

    var presetSelect = section.querySelector('#presentation-preset-select');
    var loopCheckbox = section.querySelector('#presentation-loop');
    var annotationsCheckbox = section.querySelector('#presentation-annotations-enabled');
    var recordAnnotationsCheckbox = section.querySelector('#presentation-annotations-recording');
    var scrubInput = section.querySelector('#presentation-scrub');
    var scrubLabel = section.querySelector('#presentation-scrub-label');
    var timelineEl = section.querySelector('#presentation-timeline');
    var timelineFillEl = section.querySelector('#presentation-timeline-fill');
    var timelineMarkersEl = section.querySelector('#presentation-timeline-markers');
    var timelinePlayheadEl = section.querySelector('#presentation-timeline-playhead');
    var playBtn = section.querySelector('#presentation-play-btn');
    var pauseBtn = section.querySelector('#presentation-pause-btn');
    var stopBtn = section.querySelector('#presentation-stop-btn');
    var recordQualitySelect = section.querySelector('#presentation-record-quality');
    var recordModeSelect = section.querySelector('#presentation-record-mode');
    var recordResolutionSelect = section.querySelector('#presentation-record-resolution');
    var fpsInput = section.querySelector('#presentation-fps');
    var bitrateInput = section.querySelector('#presentation-bitrate');
    var recordStartBtn = section.querySelector('#presentation-record-start-btn');
    var recordStopBtn = section.querySelector('#presentation-record-stop-btn');
    var recordResetSimCheckbox = section.querySelector('#presentation-record-reset-sim');
    var statusEl = section.querySelector('#presentation-status');
    var presentationEditorBinding = null;
    var NEW_PRESET_OPTION_VALUE = '__new_preset__';
    var timelineMarkerSignature = '';
    var timelineMarkersDirty = true;
    var timelinePointerDragging = false;
    var timelinePointerId = null;
    var timelineLoadedState = false;
    var timelineDurationState = 0;

    function clamp(v, lo, hi) {
        return Math.max(lo, Math.min(hi, v));
    }

    function setStatus(text, mode) {
        if (!statusEl) return;
        statusEl.textContent = text;
        statusEl.classList.remove('is-playing', 'is-recording', 'is-warning');
        if (mode) statusEl.classList.add(mode);
    }

    function updateScrubLabel(seconds) {
        if (!scrubLabel) return;
        scrubLabel.textContent = (seconds || 0).toFixed(1) + 's';
    }

    function getScrubDurationSeconds(fallbackSeconds) {
        var duration = parseFloat(scrubInput ? scrubInput.max : 0);
        if (!isFinite(duration) || duration <= 0) {
            duration = parseFloat(fallbackSeconds);
        }
        if (!isFinite(duration) || duration <= 0) {
            duration = 1.0;
        }
        return duration;
    }

    function updateTimelineProgressVisual(timeSeconds, durationSeconds) {
        if (!timelineEl) return;
        var safeDuration = Math.max(0.001, durationSeconds || 0.0);
        var normalized = clamp((timeSeconds || 0.0) / safeDuration, 0.0, 1.0);
        if (timelineFillEl) {
            timelineFillEl.style.transform = 'translateY(-50%) scaleX(' + normalized.toFixed(5) + ')';
        }
        if (timelinePlayheadEl) {
            timelinePlayheadEl.style.left = (normalized * 100.0).toFixed(3) + '%';
        }
    }

    function collectTimelineKeyframeMarkers(durationSeconds) {
        if (typeof getPresentationTimeline !== 'function') return [];
        var timeline = getPresentationTimeline();
        if (!timeline || !Array.isArray(timeline.tracks)) return [];

        var bucketMap = {};
        var markerIds = [];
        var safeDuration = Math.max(0.001, durationSeconds || 0.0);
        for (var i = 0; i < timeline.tracks.length; i++) {
            var track = timeline.tracks[i];
            if (!track || !Array.isArray(track.keys)) continue;
            for (var k = 0; k < track.keys.length; k++) {
                var key = track.keys[k];
                var t = parseFloat(key && key.t);
                if (!isFinite(t)) continue;
                t = clamp(t, 0.0, safeDuration);
                var bucketId = String(Math.round(t * 100.0));
                if (!bucketMap[bucketId]) {
                    bucketMap[bucketId] = { t: t, count: 1 };
                    markerIds.push(bucketId);
                } else {
                    bucketMap[bucketId].count += 1;
                    if (t < bucketMap[bucketId].t) bucketMap[bucketId].t = t;
                }
            }
        }

        var markers = [];
        for (var m = 0; m < markerIds.length; m++) {
            markers.push(bucketMap[markerIds[m]]);
        }
        markers.sort(function(a, b) { return a.t - b.t; });
        return markers;
    }

    function updateTimelineMarkers(durationSeconds) {
        if (!timelineMarkersEl) return;
        var markers = collectTimelineKeyframeMarkers(durationSeconds);
        if (!markers.length) {
            timelineMarkerSignature = '';
            timelineMarkersEl.innerHTML = '';
            return;
        }

        var signatureParts = [durationSeconds.toFixed(3)];
        for (var s = 0; s < markers.length; s++) {
            signatureParts.push(markers[s].t.toFixed(2) + ':' + markers[s].count);
        }
        var signature = signatureParts.join('|');
        if (signature === timelineMarkerSignature) return;
        timelineMarkerSignature = signature;

        var html = '';
        var safeDuration = Math.max(0.001, durationSeconds || 0.0);
        for (var i = 0; i < markers.length; i++) {
            var marker = markers[i];
            var leftPct = clamp((marker.t / safeDuration) * 100.0, 0.0, 100.0);
            var markerClass = 'presentation-timeline-marker' + (marker.count > 1 ? ' is-stack' : '');
            var markerTitle = marker.t.toFixed(2) + 's';
            if (marker.count > 1) {
                markerTitle += ' (' + marker.count + ' keys)';
            }
            html += '<span class="' + markerClass + '" style="left:' + leftPct.toFixed(3) +
                '%" title="' + markerTitle + '"></span>';
        }
        timelineMarkersEl.innerHTML = html;
    }

    function syncScrubVisual(timeSeconds, durationSeconds, forceMarkers) {
        var safeDuration = getScrubDurationSeconds(durationSeconds);
        var t = clamp(parseFloat(timeSeconds) || 0.0, 0.0, safeDuration);
        updateScrubLabel(t);
        if (forceMarkers) {
            updateTimelineMarkers(safeDuration);
            timelineMarkersDirty = false;
        }
        updateTimelineProgressVisual(t, safeDuration);
        return t;
    }

    function seekTimelineFromClientX(clientX, commitSeek) {
        if (!timelineEl || !scrubInput) return;
        var rect = timelineEl.getBoundingClientRect();
        if (!rect.width) return;
        var x = clamp(clientX - rect.left, 0.0, rect.width);
        var duration = getScrubDurationSeconds();
        var t = (x / rect.width) * duration;
        scrubInput.value = t.toFixed(3);
        syncScrubVisual(t, duration, false);
        if (commitSeek && typeof seekPresentation === 'function') {
            seekPresentation(t);
        }
    }

    function formatShortDuration(seconds) {
        var s = Math.max(0, Math.floor(seconds || 0));
        var h = Math.floor(s / 3600);
        var m = Math.floor((s % 3600) / 60);
        var sec = s % 60;
        function pad2(v) { return (v < 10 ? '0' : '') + v; }
        if (h > 0) return h + ':' + pad2(m) + ':' + pad2(sec);
        return m + ':' + pad2(sec);
    }

    function populatePresets() {
        var names = (typeof listPresentationPresets === 'function')
            ? listPresentationPresets() : [];

        presetSelect.innerHTML = '';
        var newPresetOpt = document.createElement('option');
        newPresetOpt.value = NEW_PRESET_OPTION_VALUE;
        newPresetOpt.textContent = 'New Preset';
        presetSelect.appendChild(newPresetOpt);

        for (var i = 0; i < names.length; i++) {
            var opt = document.createElement('option');
            opt.value = names[i];
            opt.textContent = names[i];
            presetSelect.appendChild(opt);
        }

        if (!names.length) {
            presetSelect.value = NEW_PRESET_OPTION_VALUE;
            return;
        }

        var defaultName = (names.indexOf('Full Feature Tour') !== -1)
            ? 'Full Feature Tour' : names[0];
        presetSelect.value = defaultName;
    }

    function updateEditorVisibilityForPresetSelection() {
        if (!presentationEditorBinding) return;
        var isNewPreset = !!presetSelect && presetSelect.value === NEW_PRESET_OPTION_VALUE;
        if (typeof presentationEditorBinding.setVisible === 'function') {
            presentationEditorBinding.setVisible(isNewPreset);
        }
        if (!isNewPreset && typeof presentationEditorBinding.setOpen === 'function') {
            presentationEditorBinding.setOpen(false);
        }
    }

    function hasRealPresetOptions() {
        if (typeof listPresentationPresets !== 'function') return false;
        var names = listPresentationPresets();
        return Array.isArray(names) && names.length > 0;
    }

    function populateRecordingQualityOptions() {
        if (!recordQualitySelect) return;

        var options = [
            { label: 'Current (no override)', value: 'current' }
        ];
        var order = ['optimal', 'high', 'ultra', 'cinematic', 'medium', 'mobile'];
        var labels = {
            mobile: 'Mobile (fastest)',
            medium: 'Medium',
            high: 'High',
            optimal: 'Optimal (recommended)',
            ultra: 'Ultra',
            cinematic: 'Cinematic (offline)'
        };
        for (var i = 0; i < order.length; i++) {
            var id = order[i];
            if (typeof QUALITY_PRESETS !== 'undefined' && QUALITY_PRESETS && QUALITY_PRESETS[id]) {
                options.push({ label: labels[id] || id, value: id });
            }
        }

        recordQualitySelect.innerHTML = '';
        for (var j = 0; j < options.length; j++) {
            var opt = document.createElement('option');
            opt.value = options[j].value;
            opt.textContent = options[j].label;
            recordQualitySelect.appendChild(opt);
        }
        recordQualitySelect.value = 'current';
    }

    function populateRecordingModeOptions() {
        if (!recordModeSelect) return;

        var options = [
            { label: 'Offline (fixed FPS)', value: 'offline' },
            { label: 'Realtime (screen capture)', value: 'realtime' }
        ];

        recordModeSelect.innerHTML = '';
        for (var i = 0; i < options.length; i++) {
            var opt = document.createElement('option');
            opt.value = options[i].value;
            opt.textContent = options[i].label;
            recordModeSelect.appendChild(opt);
        }
        recordModeSelect.value = 'offline';
    }

    function populateRecordingResolutionOptions() {
        if (!recordResolutionSelect) return;

        var options = [
            { value: 'current', label: 'Current viewport' },
            { value: '1280x720', label: '1280x720 (720p)' },
            { value: '1920x1080', label: '1920x1080 (1080p)' },
            { value: '2560x1440', label: '2560x1440 (1440p)' },
            { value: '3840x2160', label: '3840x2160 (4K UHD)' }
        ];

        recordResolutionSelect.innerHTML = '';
        for (var i = 0; i < options.length; i++) {
            var opt = document.createElement('option');
            opt.value = options[i].value;
            opt.textContent = options[i].label;
            recordResolutionSelect.appendChild(opt);
        }
        recordResolutionSelect.value = 'current';
        refreshCurrentResolutionOptionLabel();
    }

    function refreshCurrentResolutionOptionLabel() {
        if (!recordResolutionSelect) return;
        var currentOption = null;
        for (var i = 0; i < recordResolutionSelect.options.length; i++) {
            if (recordResolutionSelect.options[i].value === 'current') {
                currentOption = recordResolutionSelect.options[i];
                break;
            }
        }
        if (!currentOption) return;

        var currentWidth = 0;
        var currentHeight = 0;
        if (typeof renderer !== 'undefined' && renderer && renderer.domElement) {
            currentWidth = Math.max(1, renderer.domElement.width || 0);
            currentHeight = Math.max(1, renderer.domElement.height || 0);
        }
        if (currentWidth <= 0 || currentHeight <= 0) {
            currentWidth = Math.max(1, Math.floor((window.innerWidth || 1280) * (window.devicePixelRatio || 1)));
            currentHeight = Math.max(1, Math.floor((window.innerHeight || 720) * (window.devicePixelRatio || 1)));
        }
        currentOption.textContent = 'Current viewport (' + currentWidth + 'x' + currentHeight + ')';
    }

    function syncFromState() {
        if (typeof getPresentationState !== 'function') return;
        var state = getPresentationState();
        function setRecordingStatus(text, mode) {
            if (state.recording && state.recording_background_throttle_detected) {
                var reason = state.recording_background_throttle_reason ||
                    'Tab is hidden; browser may throttle rendering.';
                setStatus(text + ' | Warning: ' + reason, 'is-warning');
                return;
            }
            setStatus(text, mode || 'is-recording');
        }

        var durationSeconds = (typeof state.duration === 'number')
            ? Math.max(1.0, state.duration)
            : getScrubDurationSeconds();
        if (typeof state.duration === 'number') {
            scrubInput.max = durationSeconds.toFixed(2);
        }
        if (!scrubInput.matches(':active') && !timelinePointerDragging) {
            scrubInput.value = String(state.time || 0.0);
        }
        var shouldRefreshMarkers = timelineMarkersDirty ||
            timelineLoadedState !== !!state.loaded ||
            Math.abs(timelineDurationState - durationSeconds) > 1e-4;
        var scrubUiTime = parseFloat(scrubInput.value);
        if (!isFinite(scrubUiTime)) {
            scrubUiTime = state.time || 0.0;
        }
        syncScrubVisual(scrubUiTime, durationSeconds, shouldRefreshMarkers);
        timelineLoadedState = !!state.loaded;
        timelineDurationState = durationSeconds;

        if (typeof state.loop === 'boolean') {
            loopCheckbox.checked = state.loop;
        }

        if (typeof state.annotations_enabled === 'boolean') {
            annotationsCheckbox.checked = state.annotations_enabled;
        }
        if (typeof state.annotations_in_recording === 'boolean') {
            recordAnnotationsCheckbox.checked = state.annotations_in_recording;
        }
        if (recordQualitySelect) {
            if (state.recording) {
                var desiredQuality = state.recording_quality_preset || 'current';
                var hasDesired = false;
                for (var i = 0; i < recordQualitySelect.options.length; i++) {
                    if (recordQualitySelect.options[i].value === desiredQuality) {
                        hasDesired = true;
                        break;
                    }
                }
                if (hasDesired) recordQualitySelect.value = desiredQuality;
            }
            recordQualitySelect.disabled = !!state.recording;
        }
        if (recordModeSelect) {
            var offlineSupported = state.recording_offline_supported !== false;
            for (var j = 0; j < recordModeSelect.options.length; j++) {
                if (recordModeSelect.options[j].value === 'offline') {
                    recordModeSelect.options[j].disabled = !offlineSupported;
                    break;
                }
            }

            if (state.recording) {
                var desiredMode = state.recording_mode || 'realtime';
                if (!offlineSupported && desiredMode === 'offline') {
                    desiredMode = 'realtime';
                }
                var hasDesiredMode = false;
                for (var k = 0; k < recordModeSelect.options.length; k++) {
                    if (recordModeSelect.options[k].value === desiredMode) {
                        hasDesiredMode = true;
                        break;
                    }
                }
                if (hasDesiredMode) recordModeSelect.value = desiredMode;
            }
            recordModeSelect.disabled = !!state.recording;
        }
        if (recordResolutionSelect) {
            refreshCurrentResolutionOptionLabel();
            if (state.recording) {
                var desiredRes = state.recording_resolution_preset || 'current';
                var hasDesiredRes = false;
                for (var r = 0; r < recordResolutionSelect.options.length; r++) {
                    if (recordResolutionSelect.options[r].value === desiredRes) {
                        hasDesiredRes = true;
                        break;
                    }
                }
                if (hasDesiredRes) recordResolutionSelect.value = desiredRes;
            }
            recordResolutionSelect.disabled = !!state.recording;
        }

        if (state.recording) {
            if (state.recording_mode === 'offline') {
                var framesDone = state.recording_offline_frames_done || 0;
                var framesTotal = state.recording_offline_frames_total || 0;
                var renderFps = state.recording_offline_render_fps || 0;
                var queueSize = state.recording_offline_encode_queue || 0;
                var phase = state.recording_offline_phase || 'rendering';
                var sinceLastFrame = state.recording_offline_since_last_frame_s || 0;
                var outputWidth = state.recording_output_width || 0;
                var outputHeight = state.recording_output_height || 0;
                var outputResText = (outputWidth > 0 && outputHeight > 0)
                    ? (' @ ' + outputWidth + 'x' + outputHeight)
                    : '';

                if (phase === 'finalizing-encode' || phase === 'finalizing-mux' ||
                    phase === 'finalizing-download' || phase === 'done') {
                    var finalPct = Math.max(0, Math.min(100, Math.round(
                        (state.recording_offline_finalizing_progress || 0) * 100
                    )));
                    var stepText = '';
                    if (phase === 'finalizing-encode') {
                        stepText = 'flushing encoder';
                    } else if (phase === 'finalizing-mux') {
                        stepText = 'muxing WebM';
                    } else if (phase === 'finalizing-download') {
                        stepText = 'preparing download';
                    } else {
                        stepText = 'done';
                    }
                    setRecordingStatus(
                        'Finalizing video... ' + finalPct + '% (' + stepText +
                        ', ' + framesDone + ' frames rendered' + outputResText + ')',
                        'is-recording'
                    );
                } else if (phase === 'failed') {
                    setStatus('Offline render failed. Check console for details.', 'is-warning');
                } else if (sinceLastFrame >= 5.0) {
                    setRecordingStatus(
                        'Offline rendering paused (' + framesDone + ' frames, waiting ' +
                        sinceLastFrame.toFixed(1) + 's, queue ' + queueSize + outputResText + ')',
                        'is-warning'
                    );
                } else if (framesTotal > 0) {
                    var pct = Math.max(0, Math.min(100, Math.round((state.recording_offline_progress || 0) * 100)));
                    var eta = state.recording_offline_eta_s;
                    var etaText = (typeof eta === 'number' && eta >= 0)
                        ? formatShortDuration(eta)
                        : '--:--';
                    setRecordingStatus(
                        'Offline ' + pct + '% (' + framesDone + '/' + framesTotal +
                        ' frames, ' + renderFps.toFixed(1) + ' fps, ETA ' + etaText +
                        ', queue ' + queueSize + outputResText + ')',
                        'is-recording'
                    );
                } else {
                    var timelineDone = state.recording_offline_timeline_done_s || 0;
                    setRecordingStatus(
                        'Offline rendering... ' + framesDone + ' frames (' +
                        timelineDone.toFixed(1) + 's timeline, ' + renderFps.toFixed(1) +
                        ' fps, queue ' + queueSize + outputResText + ')',
                        'is-recording'
                    );
                }
            } else {
                var realtimeW = state.recording_output_width || 0;
                var realtimeH = state.recording_output_height || 0;
                if (realtimeW > 0 && realtimeH > 0) {
                    setRecordingStatus('Realtime screen capture @ ' + realtimeW + 'x' + realtimeH, 'is-recording');
                } else {
                    setRecordingStatus('Realtime screen capture...', 'is-recording');
                }
            }
        } else if (recordModeSelect && recordModeSelect.value === 'offline' &&
            state.recording_offline_supported === false) {
            var reason = state.recording_offline_unavailable_reason || 'Offline mode not supported in this browser.';
            setStatus(reason, 'is-warning');
        } else if (state.playing) {
            setStatus('Playing ' + (state.name || 'timeline'), 'is-playing');
        } else if (state.loaded) {
            setStatus('Loaded: ' + (state.name || 'timeline'), '');
        } else if (presetSelect && presetSelect.value === NEW_PRESET_OPTION_VALUE) {
            setStatus('New preset mode: edit timeline then press APPLY.', '');
        } else if (state.presets_loading) {
            setStatus('Loading presets...', '');
        } else if (hasRealPresetOptions()) {
            setStatus('Preset ready. Press PLAY to apply.', '');
        } else if (presetSelect && presetSelect.options && presetSelect.options.length > 0) {
            setStatus('No presets found.', 'is-warning');
        } else {
            setStatus('Idle', '');
        }
    }

    function loadSelectedPreset() {
        if (typeof loadPresentationPreset !== 'function') return;
        if (!presetSelect.value) {
            setStatus('No preset selected.', 'is-warning');
            return;
        }

        if (presetSelect.value === NEW_PRESET_OPTION_VALUE) {
            if (presentationEditorBinding) {
                if (typeof presentationEditorBinding.setVisible === 'function') {
                    presentationEditorBinding.setVisible(true);
                }
                if (typeof presentationEditorBinding.startNewPresetDraft === 'function') {
                    presentationEditorBinding.startNewPresetDraft();
                } else if (typeof presentationEditorBinding.setOpen === 'function') {
                    presentationEditorBinding.setOpen(true);
                }
            }
            setStatus('New preset mode: edit timeline then press APPLY.', '');
            timelineMarkersDirty = true;
            syncFromState();
            return;
        }

        updateEditorVisibilityForPresetSelection();

        if (!loadPresentationPreset(presetSelect.value)) {
            setStatus('Failed to load preset.', 'is-warning');
            return;
        }
        if (typeof setPresentationLoop === 'function') {
            setPresentationLoop(loopCheckbox.checked);
        }
        if (typeof setPresentationAnnotationsEnabled === 'function') {
            setPresentationAnnotationsEnabled(annotationsCheckbox.checked);
        }
        if (typeof setPresentationAnnotationsIncludedInRecording === 'function') {
            setPresentationAnnotationsIncludedInRecording(recordAnnotationsCheckbox.checked);
        }
        if (section && typeof CustomEvent === 'function') {
            section.dispatchEvent(new CustomEvent('presentation:timeline-loaded', {
                detail: {
                    source: 'preset',
                    name: presetSelect.value
                }
            }));
        }
        timelineMarkersDirty = true;
        syncFromState();
    }

    presetSelect.addEventListener('change', loadSelectedPreset);

    loopCheckbox.addEventListener('change', function() {
        if (typeof setPresentationLoop === 'function') {
            setPresentationLoop(loopCheckbox.checked);
        }
    });

    annotationsCheckbox.addEventListener('change', function() {
        if (typeof setPresentationAnnotationsEnabled === 'function') {
            setPresentationAnnotationsEnabled(annotationsCheckbox.checked);
        }
    });

    recordAnnotationsCheckbox.addEventListener('change', function() {
        if (typeof setPresentationAnnotationsIncludedInRecording === 'function') {
            setPresentationAnnotationsIncludedInRecording(recordAnnotationsCheckbox.checked);
        }
    });

    scrubInput.addEventListener('input', function() {
        syncScrubVisual(parseFloat(scrubInput.value) || 0.0, getScrubDurationSeconds(), false);
    });
    scrubInput.addEventListener('change', function() {
        var t = parseFloat(scrubInput.value) || 0.0;
        syncScrubVisual(t, getScrubDurationSeconds(), false);
        if (typeof seekPresentation === 'function') {
            seekPresentation(t);
        }
    });

    if (timelineEl) {
        function endTimelinePointerDrag(event, commitSeek) {
            if (!timelinePointerDragging) return;
            if (timelinePointerId !== null && event.pointerId !== timelinePointerId) return;
            if (commitSeek) {
                seekTimelineFromClientX(event.clientX, true);
            }
            timelinePointerDragging = false;
            timelinePointerId = null;
            timelineEl.classList.remove('is-dragging');
            if (timelineEl.releasePointerCapture) {
                try {
                    timelineEl.releasePointerCapture(event.pointerId);
                } catch (err) {
                    // Ignore when capture is already released.
                }
            }
            event.preventDefault();
        }

        timelineEl.addEventListener('pointerdown', function(event) {
            if (event.button !== 0) return;
            timelinePointerDragging = true;
            timelinePointerId = event.pointerId;
            timelineEl.classList.add('is-dragging');
            if (timelineEl.setPointerCapture) {
                timelineEl.setPointerCapture(event.pointerId);
            }
            seekTimelineFromClientX(event.clientX, true);
            event.preventDefault();
        });

        timelineEl.addEventListener('pointermove', function(event) {
            if (!timelinePointerDragging) return;
            if (timelinePointerId !== null && event.pointerId !== timelinePointerId) return;
            seekTimelineFromClientX(event.clientX, true);
            event.preventDefault();
        });

        timelineEl.addEventListener('pointerup', function(event) {
            endTimelinePointerDrag(event, true);
        });

        timelineEl.addEventListener('pointercancel', function(event) {
            endTimelinePointerDrag(event, false);
        });
    }

    playBtn.addEventListener('click', function() {
        if (typeof getPresentationState !== 'function') return;
        var state = getPresentationState();
        if (!state.loaded) {
            if (typeof ensurePresentationPresetsLoaded === 'function') {
                setStatus('Loading presets...', '');
                var waitForPresets = ensurePresentationPresetsLoaded();
                if (waitForPresets && typeof waitForPresets.then === 'function') {
                    waitForPresets.then(function() {
                        loadSelectedPreset();
                        if (presetSelect.value !== NEW_PRESET_OPTION_VALUE &&
                            typeof playPresentation === 'function') {
                            playPresentation(false);
                        }
                        syncFromState();
                    });
                    return;
                }
            }
            loadSelectedPreset();
            if (presetSelect.value === NEW_PRESET_OPTION_VALUE) {
                syncFromState();
                return;
            }
        }
        if (typeof playPresentation === 'function') playPresentation(false);
        syncFromState();
    });

    pauseBtn.addEventListener('click', function() {
        if (typeof pausePresentation === 'function') pausePresentation();
        syncFromState();
    });

    stopBtn.addEventListener('click', function() {
        if (typeof stopPresentation === 'function') stopPresentation();
        syncFromState();
    });

    recordStartBtn.addEventListener('click', function() {
        if (typeof startPresentationRecording !== 'function') return;
        if (recordResetSimCheckbox && recordResetSimCheckbox.checked) {
            if (typeof observer !== 'undefined' && observer) {
                observer.time = 0.0;
            }
            if (typeof shader !== 'undefined' && shader) {
                shader.needsUpdate = true;
            }
        }
        var fps = clamp(parseFloat(fpsInput.value) || 60, 24, 120);
        var bitrate = clamp(parseFloat(bitrateInput.value) || 20, 4, 80);
        fpsInput.value = Math.round(fps);
        bitrateInput.value = Math.round(bitrate);
        var started = startPresentationRecording({
            fps: fps,
            bitrateMbps: bitrate,
            autoStopOnPresentationEnd: true,
            recordingMode: recordModeSelect ? recordModeSelect.value : 'offline',
            recordingResolution: recordResolutionSelect ? recordResolutionSelect.value : 'current',
            qualityPreset: recordQualitySelect ? recordQualitySelect.value : 'current',
            includeAnnotationsInRecording: !!recordAnnotationsCheckbox.checked
        });
        if (!started && typeof getPresentationState === 'function') {
            var stateAfterStart = getPresentationState();
            if (stateAfterStart.recording_mode_preferred === 'offline' &&
                stateAfterStart.recording_offline_unavailable_reason) {
                setStatus(stateAfterStart.recording_offline_unavailable_reason, 'is-warning');
            } else {
                setStatus('Failed to start recording.', 'is-warning');
            }
        }
        syncFromState();
    });

    recordStopBtn.addEventListener('click', function() {
        if (typeof stopPresentationRecording === 'function') {
            stopPresentationRecording();
        }
        syncFromState();
    });

    function initializePresetsUI() {
        populateRecordingQualityOptions();
        populateRecordingModeOptions();
        populateRecordingResolutionOptions();
        if (typeof ensurePresentationPresetsLoaded === 'function') {
            setStatus('Loading presets...', '');
            var loading = ensurePresentationPresetsLoaded();
            if (loading && typeof loading.then === 'function') {
                loading.then(function() {
                    populatePresets();
                    if (!hasRealPresetOptions()) {
                        setStatus('No presets found.', 'is-warning');
                    } else {
                        setStatus('Preset ready. Press PLAY to apply.', '');
                    }
                    syncFromState();
                }).catch(function() {
                    setStatus('Failed to load preset files.', 'is-warning');
                });
                return;
            }
        }

        populatePresets();
        if (!hasRealPresetOptions()) {
            setStatus('No presets found.', 'is-warning');
        } else {
            setStatus('Preset ready. Press PLAY to apply.', '');
        }
        syncFromState();
    }

    initializePresetsUI();

    if (typeof bindPresentationTimelineEditor === 'function') {
        presentationEditorBinding = bindPresentationTimelineEditor(section, {
            setStatus: setStatus,
            syncFromState: syncFromState,
            getSelectedPresetName: function() {
                return presetSelect ? presetSelect.value : '';
            },
            loadSelectedPreset: loadSelectedPreset
        });
    }

    if (presentationEditorBinding && typeof presentationEditorBinding.syncFromRuntime === 'function') {
        presentationEditorBinding.syncFromRuntime(true);
    }
    updateEditorVisibilityForPresetSelection();
    section.addEventListener('presentation:timeline-loaded', function() {
        timelineMarkersDirty = true;
        syncFromState();
    });

    if (section._presentationSyncTimer) {
        clearInterval(section._presentationSyncTimer);
    }
    section._presentationSyncTimer = window.setInterval(syncFromState, 120);
}
