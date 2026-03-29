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
                '<div class="presentation-desc">Playback, toggles &amp; recording are in the <strong>&#9650; TIMELINE</strong> panel at the bottom.</div>' +

                '<div class="presentation-row">' +
                    '<label for="presentation-preset-select">Preset</label>' +
                    '<select id="presentation-preset-select" class="presentation-select"></select>' +
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
    var statusEl     = section.querySelector('#presentation-status');
    var NEW_PRESET_OPTION_VALUE = '__new_preset__';
    var timelineMarkersDirty = true; // kept so external code that sets this still works

    function setStatus(text, mode) {
        if (!statusEl) return;
        statusEl.textContent = text;
        statusEl.classList.remove('is-playing', 'is-recording', 'is-warning');
        if (mode) statusEl.classList.add(mode);
    }

    function hasRealPresetOptions() {
        if (typeof listPresentationPresets !== 'function') return false;
        var names = listPresentationPresets();
        return Array.isArray(names) && names.length > 0;
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
        // Old editor removed — now handled by the bottom timeline panel
    }

    function syncFromState() {
        if (typeof getPresentationState !== 'function') return;
        var s = getPresentationState();
        if (s.recording) setStatus('Recording\u2026', 'is-recording');
        else if (s.playing) setStatus('Playing \u2013 ' + (s.name || ''), 'is-playing');
        else if (s.loaded) setStatus('Loaded: ' + (s.name || ''), '');
        else if (s.presets_loading) setStatus('Loading presets\u2026', '');
        else setStatus('Idle', '');
    }

    function loadSelectedPreset() {
        var v = presetSelect.value;
        if (!v || v === NEW_PRESET_OPTION_VALUE) {
            if (typeof timelinePanelBinding !== 'undefined' && timelinePanelBinding) {
                timelinePanelBinding.open();
            }
            setStatus('Open the \u25b2 TIMELINE panel to edit.', '');
            return;
        }
        if (typeof loadPresentationPreset !== 'function' || !loadPresentationPreset(v)) {
            setStatus('Failed to load: ' + v, 'is-warning');
            return;
        }
        syncFromState();
    }

    presetSelect.addEventListener('change', loadSelectedPreset);

    function initializePresetsUI() {
        if (typeof ensurePresentationPresetsLoaded === 'function') {
            setStatus('Loading presets\u2026', '');
            var loading = ensurePresentationPresetsLoaded();
            if (loading && typeof loading.then === 'function') {
                loading.then(function() {
                    populatePresets();
                    setStatus(hasRealPresetOptions() ? 'Select a preset below.' : 'No presets found.', hasRealPresetOptions() ? '' : 'is-warning');
                    syncFromState();
                }).catch(function() {
                    setStatus('Failed to load preset files.', 'is-warning');
                });
                return;
            }
        }
        populatePresets();
        setStatus(hasRealPresetOptions() ? 'Select a preset below.' : 'No presets found.', hasRealPresetOptions() ? '' : 'is-warning');
        syncFromState();
    }

    initializePresetsUI();

    if (section._presentationSyncTimer) clearInterval(section._presentationSyncTimer);
    section._presentationSyncTimer = window.setInterval(syncFromState, 500);
}
