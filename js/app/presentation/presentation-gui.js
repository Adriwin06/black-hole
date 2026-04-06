// Role: Presentation controls embedded in the ANIMATIONS panel.
//       Handles timeline preset selection, playback, recording, and
//       annotation toggles without using dat.GUI rows.

import {
    ensurePresentationPresetsLoaded,
    listPresentationPresets,
    getPresentationState,
    loadPresentationPreset
} from './presentation-controller.js';
import { getTimelinePanelBinding } from '../core/ui-bindings.js';

export function createPresentationAnimationSectionHtml() {
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

export function bindPresentationAnimationSection(panelRoot) {
    if (!panelRoot) return;

    var section = panelRoot.querySelector('#presentation-section');
    if (!section) return;

    var presetSelect = section.querySelector('#presentation-preset-select');
    var statusEl = section.querySelector('#presentation-status');
    var NEW_PRESET_OPTION_VALUE = '__new_preset__';

    function setStatus(text, mode) {
        if (!statusEl) return;
        statusEl.textContent = text;
        statusEl.classList.remove('is-playing', 'is-recording', 'is-warning');
        if (mode) statusEl.classList.add(mode);
    }

    function hasRealPresetOptions() {
        var names = listPresentationPresets();
        return Array.isArray(names) && names.length > 0;
    }

    function populatePresets() {
        var names = listPresentationPresets();

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
            ? 'Full Feature Tour'
            : names[0];
        presetSelect.value = defaultName;
    }

    function syncFromState() {
        var s = getPresentationState();
        if (s.recording) setStatus('Recording...', 'is-recording');
        else if (s.playing) setStatus('Playing - ' + (s.name || ''), 'is-playing');
        else if (s.loaded) setStatus('Loaded: ' + (s.name || ''), '');
        else if (s.presets_loading) setStatus('Loading presets...', '');
        else setStatus('Idle', '');
    }

    function loadSelectedPreset() {
        var value = presetSelect.value;
        if (!value || value === NEW_PRESET_OPTION_VALUE) {
            var timelineBinding = getTimelinePanelBinding();
            if (timelineBinding && typeof timelineBinding.open === 'function') {
                timelineBinding.open();
            }
            setStatus('Open the ▲ TIMELINE panel to edit.', '');
            return;
        }
        if (!loadPresentationPreset(value)) {
            setStatus('Failed to load: ' + value, 'is-warning');
            return;
        }
        syncFromState();
    }

    presetSelect.addEventListener('change', loadSelectedPreset);

    function initializePresetsUI() {
        setStatus('Loading presets...', '');
        var loading = ensurePresentationPresetsLoaded();
        if (loading && typeof loading.then === 'function') {
            loading.then(function() {
                populatePresets();
                setStatus(
                    hasRealPresetOptions() ? 'Select a preset below.' : 'No presets found.',
                    hasRealPresetOptions() ? '' : 'is-warning'
                );
                syncFromState();
            }).catch(function() {
                setStatus('Failed to load preset files.', 'is-warning');
            });
            return;
        }

        populatePresets();
        setStatus(
            hasRealPresetOptions() ? 'Select a preset below.' : 'No presets found.',
            hasRealPresetOptions() ? '' : 'is-warning'
        );
        syncFromState();
    }

    initializePresetsUI();

    if (section._presentationSyncTimer) clearInterval(section._presentationSyncTimer);
    section._presentationSyncTimer = window.setInterval(syncFromState, 500);
}
