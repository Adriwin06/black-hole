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

                '<div class="presentation-row">' +
                    '<label for="presentation-scrub">Time</label>' +
                    '<input id="presentation-scrub" type="range" min="0" max="1" step="0.01" value="0">' +
                    '<span id="presentation-scrub-label">0.0s</span>' +
                '</div>' +

                '<div class="presentation-btn-row">' +
                    '<button id="presentation-play-btn" class="presentation-btn presentation-btn-primary" type="button">&#9654; PLAY</button>' +
                    '<button id="presentation-pause-btn" class="presentation-btn" type="button">&#10074;&#10074; PAUSE</button>' +
                    '<button id="presentation-stop-btn" class="presentation-btn" type="button">&#9632; STOP</button>' +
                '</div>' +

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
    var statusEl = section.querySelector('#presentation-status');

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
        for (var i = 0; i < names.length; i++) {
            var opt = document.createElement('option');
            opt.value = names[i];
            opt.textContent = names[i];
            presetSelect.appendChild(opt);
        }

        if (!names.length) return;

        var defaultName = (names.indexOf('Full Feature Tour') !== -1)
            ? 'Full Feature Tour' : names[0];
        presetSelect.value = defaultName;
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
            if (state.recording && state.recording_background_throttle_risk) {
                var reason = state.recording_background_throttle_reason ||
                    'Tab is hidden; browser may throttle rendering.';
                setStatus(text + ' | Warning: ' + reason, 'is-warning');
                return;
            }
            setStatus(text, mode || 'is-recording');
        }

        if (typeof state.duration === 'number') {
            scrubInput.max = Math.max(1.0, state.duration).toFixed(2);
        }
        if (!scrubInput.matches(':active')) {
            scrubInput.value = state.time || 0.0;
        }
        updateScrubLabel(state.time || 0.0);

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
        } else if (state.presets_loading) {
            setStatus('Loading presets...', '');
        } else if (presetSelect && presetSelect.options && presetSelect.options.length > 0) {
            setStatus('Preset ready. Press PLAY to apply.', '');
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
        updateScrubLabel(parseFloat(scrubInput.value) || 0.0);
    });
    scrubInput.addEventListener('change', function() {
        if (typeof seekPresentation === 'function') {
            seekPresentation(parseFloat(scrubInput.value) || 0.0);
        }
    });

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
                        if (typeof playPresentation === 'function') playPresentation(false);
                        syncFromState();
                    });
                    return;
                }
            }
            loadSelectedPreset();
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
                    if (!presetSelect.options.length) {
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
        if (!presetSelect.options.length) {
            setStatus('No presets found.', 'is-warning');
        } else {
            setStatus('Preset ready. Press PLAY to apply.', '');
        }
        syncFromState();
    }

    initializePresetsUI();

    if (section._presentationSyncTimer) {
        clearInterval(section._presentationSyncTimer);
    }
    section._presentationSyncTimer = window.setInterval(syncFromState, 120);
}
