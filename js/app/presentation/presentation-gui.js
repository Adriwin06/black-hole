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

    function syncFromState() {
        if (typeof getPresentationState !== 'function') return;
        var state = getPresentationState();

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

        if (state.recording) {
            setStatus('Recording timeline output...', 'is-recording');
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
        startPresentationRecording({
            fps: fps,
            bitrateMbps: bitrate,
            autoStopOnPresentationEnd: true,
            includeAnnotationsInRecording: !!recordAnnotationsCheckbox.checked
        });
        syncFromState();
    });

    recordStopBtn.addEventListener('click', function() {
        if (typeof stopPresentationRecording === 'function') {
            stopPresentationRecording();
        }
        syncFromState();
    });

    function initializePresetsUI() {
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
