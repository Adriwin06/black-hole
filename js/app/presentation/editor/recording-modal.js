import { observer, shader } from '../../core/runtime/runtime-state.js';
import { QUALITY_PRESETS } from '../../ui/quality-presets.js';

export function setupRecordingModal(options) {
    var panel = options.panel;
    var esc = options.esc;
    var rememberState = options.rememberState;
    var applySavedRecordingUiState = options.applySavedRecordingUiState;
    var rebuildTrackList = options.rebuildTrackList;
    var getSelectedTrack = options.getSelectedTrack;
    var getPresentationState = options.getPresentationState;
    var getPresentationParamHudState = options.getPresentationParamHudState;
    var isParamInHud = options.isParamInHud;
    var setPresentationLoop = options.setPresentationLoop;
    var setPresentationAnnotationsEnabled = options.setPresentationAnnotationsEnabled;
    var setPresentationAnnotationsIncludedInRecording =
        options.setPresentationAnnotationsIncludedInRecording;
    var setPresentationParamHudEnabled = options.setPresentationParamHudEnabled;
    var setPresentationParamHudIncludedInRecording =
        options.setPresentationParamHudIncludedInRecording;
    var removeParamFromHud = options.removeParamFromHud;
    var addParamToHud = options.addParamToHud;
    var clearParamHud = options.clearParamHud;
    var setParamHudLayout = options.setParamHudLayout;
    var capturePresentationScreenshot = options.capturePresentationScreenshot;
    var startPresentationRecording = options.startPresentationRecording;
    var stopPresentationRecording = options.stopPresentationRecording;

    var elements = options.elements || {};
    var recBtn = elements.recBtn;
    var recModal = elements.recModal;
    var recCloseBtn = elements.recCloseBtn;
    var recLoopCb = elements.recLoopCb;
    var recAnnotCb = elements.recAnnotCb;
    var recAnnotRecordCb = elements.recAnnotRecordCb;
    var recParamHudShowCb = elements.recParamHudShowCb;
    var recParamHudRecordCb = elements.recParamHudRecordCb;
    var hudXInput = elements.hudXInput;
    var hudYInput = elements.hudYInput;
    var hudFsInput = elements.hudFsInput;
    var hudLayoutRow = elements.hudLayoutRow;
    var recHudEmptyEl = elements.recHudEmptyEl;
    var recHudListEl = elements.recHudListEl;
    var recAddSelectedBtn = elements.recAddSelectedBtn;
    var recClearHudBtn = elements.recClearHudBtn;
    var recResetSimCb = elements.recResetSimCb;
    var recQualitySelect = elements.recQualitySelect;
    var recModeSelect = elements.recModeSelect;
    var recResSelect = elements.recResSelect;
    var recFpsInput = elements.recFpsInput;
    var recBitrateInput = elements.recBitrateInput;
    var recShotBtn = elements.recShotBtn;
    var recStartBtn = elements.recStartBtn;
    var recStopBtn = elements.recStopBtn;
    var recStatusEl = elements.recStatusEl;

    function clampNum(v, lo, hi) {
        return Math.min(hi, Math.max(lo, v));
    }

    function formatRecDuration(s) {
        if (isNaN(s) || s < 0) return '--:--';
        var m = Math.floor(s / 60);
        var sec = Math.floor(s % 60);
        return (m < 10 ? '0' : '') + m + ':' + (sec < 10 ? '0' : '') + sec;
    }

    function setRecStatus(text, cls) {
        if (!recStatusEl) return;
        recStatusEl.textContent = text;
        recStatusEl.className = 'tl-rec-status' + (cls ? ' ' + cls : '');
    }

    function populateRecQuality() {
        if (!recQualitySelect) return;
        var cur = recQualitySelect.value;
        recQualitySelect.innerHTML = '';
        var presets = (QUALITY_PRESETS && typeof QUALITY_PRESETS === 'object')
            ? Object.keys(QUALITY_PRESETS)
            : [];
        var keys = presets.length
            ? presets
            : ['optimal', 'high', 'ultra', 'cinematic', 'medium', 'mobile'];
        for (var i = 0; i < keys.length; i++) {
            var option = document.createElement('option');
            option.value = keys[i];
            option.textContent = keys[i].charAt(0).toUpperCase() + keys[i].slice(1);
            recQualitySelect.appendChild(option);
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
        var modes = [
            { v: 'offline', l: 'Offline (fixed FPS)' },
            { v: 'realtime', l: 'Realtime (screen capture)' }
        ];
        for (var i = 0; i < modes.length; i++) {
            var option = document.createElement('option');
            option.value = modes[i].v;
            option.textContent = modes[i].l;
            recModeSelect.appendChild(option);
        }
        if (cur) recModeSelect.value = cur;
        if (!recModeSelect.value) recModeSelect.value = 'offline';
    }

    function refreshRecResolutionLabel() {
        if (!recResSelect || typeof getPresentationState !== 'function') return;
        var state = getPresentationState() || {};
        for (var i = 0; i < recResSelect.options.length; i++) {
            var option = recResSelect.options[i];
            if (option.value === 'current') {
                var w = state.recording_output_width || window.innerWidth;
                var h = state.recording_output_height || window.innerHeight;
                option.textContent = 'Current viewport (' + w + '\xd7' + h + ')';
                break;
            }
        }
    }

    function populateRecResolution() {
        if (!recResSelect) return;
        var cur = recResSelect.value;
        recResSelect.innerHTML = '';
        var options = [
            { v: 'current', l: 'Current viewport' },
            { v: '1280x720', l: '1280\xd7720 (HD)' },
            { v: '1920x1080', l: '1920\xd71080 (Full HD)' },
            { v: '2560x1440', l: '2560\xd71440 (2K)' },
            { v: '3840x2160', l: '3840\xd72160 (4K)' },
            { v: '7680x4320', l: '7680\xd74320 (8K)' }
        ];
        for (var i = 0; i < options.length; i++) {
            var option = document.createElement('option');
            option.value = options[i].v;
            option.textContent = options[i].l;
            recResSelect.appendChild(option);
        }
        if (cur) recResSelect.value = cur;
        if (!recResSelect.value) recResSelect.value = 'current';
        refreshRecResolutionLabel();
    }

    function renderRecHudList() {
        if (!recHudListEl || !recHudEmptyEl || typeof getPresentationParamHudState !== 'function') {
            return;
        }
        var hudState = getPresentationParamHudState() || { items: [] };
        var items = Array.isArray(hudState.items) ? hudState.items : [];
        if (recClearHudBtn) recClearHudBtn.disabled = items.length === 0;

        var selectedTrack = typeof getSelectedTrack === 'function' ? getSelectedTrack() : '';
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
        if (!recModal || !recBtn || typeof getPresentationState !== 'function') return;

        if (!recModal.classList.contains('is-open')) {
            var closedState = getPresentationState() || {};
            recBtn.classList.toggle('is-recording', !!closedState.recording);
            return;
        }

        var state = getPresentationState() || {};

        if (recLoopCb) recLoopCb.checked = !!state.loop;
        if (recAnnotCb) recAnnotCb.checked = !!state.annotations_enabled;
        if (recAnnotRecordCb) recAnnotRecordCb.checked = !!state.annotations_in_recording;
        if (recParamHudShowCb) recParamHudShowCb.checked = !!state.param_hud_enabled;
        if (recParamHudRecordCb) recParamHudRecordCb.checked = !!state.param_hud_in_recording;

        if (hudLayoutRow) hudLayoutRow.style.display = (state.param_hud_count > 0) ? '' : 'none';
        if (typeof getPresentationParamHudState === 'function') {
            var hudState = getPresentationParamHudState();
            if (hudXInput && document.activeElement !== hudXInput) {
                hudXInput.value = hudState.anchorX.toFixed(2);
            }
            if (hudYInput && document.activeElement !== hudYInput) {
                hudYInput.value = hudState.anchorY.toFixed(2);
            }
            if (hudFsInput && document.activeElement !== hudFsInput) {
                hudFsInput.value = hudState.fontSize;
            }
        }

        renderRecHudList();
        refreshRecResolutionLabel();

        recBtn.classList.toggle('is-recording', !!state.recording);
        if (recShotBtn) recShotBtn.disabled = !!state.recording;

        if (!state.recording) {
            if (recStartBtn) recStartBtn.disabled = false;
            if (recStopBtn) recStopBtn.disabled = true;
            if (state.recording_offline_unavailable_reason) {
                setRecStatus(state.recording_offline_unavailable_reason, 'is-warning');
            } else {
                setRecStatus('Idle', '');
            }
            return;
        }

        if (recStartBtn) recStartBtn.disabled = true;
        if (recStopBtn) recStopBtn.disabled = false;

        if (state.recording_mode === 'realtime') {
            if (state.recording_background_throttle_detected) {
                setRecStatus('\u26a0 Background throttle detected \u2014 keep window focused!', 'is-warning');
            } else {
                setRecStatus('Recording\u2026 (realtime)', 'is-recording');
            }
            return;
        }

        var phase = state.recording_offline_phase || '';
        if (phase === 'rendering') {
            var pct = Math.round((state.recording_offline_progress || 0) * 100);
            var done = state.recording_offline_frames_done || 0;
            var total = state.recording_offline_frames_total || 0;
            var fps = state.recording_offline_render_fps
                ? (' @ ' + state.recording_offline_render_fps.toFixed(1) + ' fps')
                : '';
            var eta = (state.recording_offline_eta_s != null && state.recording_offline_eta_s >= 0)
                ? (' ETA ' + formatRecDuration(state.recording_offline_eta_s))
                : '';
            setRecStatus('Rendering ' + pct + '% (' + done + '/' + total + ')' + fps + eta, 'is-recording');
        } else if (
            phase === 'finalizing-encode' ||
            phase === 'finalizing-mux' ||
            phase === 'finalizing-download'
        ) {
            var progress = state.recording_offline_finalizing_progress;
            var label = phase === 'finalizing-encode'
                ? 'Encoding'
                : phase === 'finalizing-mux'
                    ? 'Muxing'
                    : 'Downloading';
            var pctStr = (progress != null && progress >= 0)
                ? (' ' + Math.round(progress * 100) + '%')
                : '\u2026';
            setRecStatus(label + pctStr, 'is-recording');
        } else {
            setRecStatus('Recording\u2026', 'is-recording');
        }
    }

    function openRecModal() {
        populateRecQuality();
        populateRecMode();
        populateRecResolution();
        if (typeof applySavedRecordingUiState === 'function') {
            applySavedRecordingUiState();
        }
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
        if (typeof setPresentationLoop === 'function') {
            setPresentationLoop(recLoopCb.checked);
        }
        rememberState();
    });
    recAnnotCb.addEventListener('change', function() {
        if (typeof setPresentationAnnotationsEnabled === 'function') {
            setPresentationAnnotationsEnabled(recAnnotCb.checked);
        }
        rememberState();
    });
    recAnnotRecordCb.addEventListener('change', function() {
        if (typeof setPresentationAnnotationsIncludedInRecording === 'function') {
            setPresentationAnnotationsIncludedInRecording(recAnnotRecordCb.checked);
        }
        rememberState();
    });
    recParamHudShowCb.addEventListener('change', function() {
        if (typeof setPresentationParamHudEnabled === 'function') {
            setPresentationParamHudEnabled(recParamHudShowCb.checked);
        }
        rememberState();
    });
    recParamHudRecordCb.addEventListener('change', function() {
        if (typeof setPresentationParamHudIncludedInRecording === 'function') {
            setPresentationParamHudIncludedInRecording(recParamHudRecordCb.checked);
        }
        rememberState();
    });

    if (recHudListEl) {
        recHudListEl.addEventListener('click', function(e) {
            var chip = e.target.closest('[data-remove-hud-path]');
            if (!chip) return;
            var path = chip.getAttribute('data-remove-hud-path');
            if (path && typeof removeParamFromHud === 'function') {
                removeParamFromHud(path);
            }
            syncRecModal();
            if (typeof rebuildTrackList === 'function') rebuildTrackList();
            rememberState();
        });
    }
    if (recAddSelectedBtn) {
        recAddSelectedBtn.addEventListener('click', function() {
            var selectedTrack = typeof getSelectedTrack === 'function' ? getSelectedTrack() : '';
            if (!selectedTrack || typeof addParamToHud !== 'function') return;
            addParamToHud(selectedTrack, selectedTrack);
            syncRecModal();
            if (typeof rebuildTrackList === 'function') rebuildTrackList();
            rememberState();
        });
    }
    if (recClearHudBtn) {
        recClearHudBtn.addEventListener('click', function() {
            if (typeof clearParamHud !== 'function') return;
            clearParamHud();
            syncRecModal();
            if (typeof rebuildTrackList === 'function') rebuildTrackList();
            rememberState();
        });
    }

    function applyHudLayout() {
        if (typeof setParamHudLayout !== 'function') return;
        var x = parseFloat(hudXInput.value);
        var y = parseFloat(hudYInput.value);
        var fs = parseFloat(hudFsInput.value);
        setParamHudLayout({
            anchorX: isFinite(x) ? x : undefined,
            anchorY: isFinite(y) ? y : undefined,
            fontSize: isFinite(fs) ? fs : undefined
        });
        syncRecModal();
        rememberState();
    }

    if (hudXInput) hudXInput.addEventListener('input', applyHudLayout);
    if (hudYInput) hudYInput.addEventListener('input', applyHudLayout);
    if (hudFsInput) hudFsInput.addEventListener('input', applyHudLayout);
    if (recQualitySelect) recQualitySelect.addEventListener('change', rememberState);
    if (recModeSelect) recModeSelect.addEventListener('change', rememberState);
    if (recResSelect) recResSelect.addEventListener('change', rememberState);
    if (recFpsInput) recFpsInput.addEventListener('input', rememberState);
    if (recBitrateInput) recBitrateInput.addEventListener('input', rememberState);
    if (recResetSimCb) recResetSimCb.addEventListener('change', rememberState);

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
            var shotState = (typeof getPresentationState === 'function')
                ? getPresentationState()
                : {};
            setRecStatus(
                shotState.recording_offline_unavailable_reason ||
                    'Failed to capture screenshot.',
                'is-warning'
            );
        } else {
            setRecStatus('Offline screenshot downloaded (.png).', '');
        }
        syncRecModal();
    });

    recStartBtn.addEventListener('click', async function() {
        if (typeof startPresentationRecording !== 'function') return;

        if (recResetSimCb && recResetSimCb.checked) {
            observer.time = 0.0;
            shader.needsUpdate = true;
        }

        var fps = clampNum(parseFloat(recFpsInput.value) || 60, 24, 120);
        var bitrate = clampNum(parseFloat(recBitrateInput.value) || 20, 4, 80);
        recFpsInput.value = Math.round(fps);
        recBitrateInput.value = Math.round(bitrate);

        var recMode = recModeSelect ? recModeSelect.value : 'offline';
        var recOptions = {
            fps: fps,
            bitrateMbps: bitrate,
            autoStopOnPresentationEnd: true,
            recordingMode: recMode,
            recordingResolution: recResSelect ? recResSelect.value : 'current',
            qualityPreset: recQualitySelect ? recQualitySelect.value : 'optimal',
            includeAnnotationsInRecording: !!recAnnotRecordCb.checked
        };

        if (
            recMode === 'offline' &&
            typeof showSaveFilePicker === 'function' &&
            typeof window.WebMMuxer !== 'undefined' &&
            typeof window.WebMMuxer.FileSystemWritableFileStreamTarget === 'function'
        ) {
            try {
                var fileHandle = await showSaveFilePicker({
                    suggestedName: 'black-hole-presentation.webm',
                    types: [{
                        description: 'WebM video',
                        accept: { 'video/webm': ['.webm'] }
                    }]
                });
                recOptions.writableFileStream = await fileHandle.createWritable();
            } catch (pickerErr) {
                if (pickerErr.name === 'AbortError') return;
                console.warn(
                    'File picker failed, falling back to in-memory recording:',
                    pickerErr
                );
            }
        }

        var started = startPresentationRecording(recOptions);
        if (!started) {
            if (recOptions.writableFileStream) {
                try { recOptions.writableFileStream.abort(); } catch (e) {}
            }
            var startState = (typeof getPresentationState === 'function')
                ? getPresentationState()
                : {};
            setRecStatus(
                startState.recording_offline_unavailable_reason ||
                    'Failed to start recording.',
                'is-warning'
            );
        }
        syncRecModal();
    });

    recStopBtn.addEventListener('click', function() {
        if (typeof stopPresentationRecording === 'function') {
            stopPresentationRecording();
        }
        syncRecModal();
    });

    panel.addEventListener('click', function(e) {
        if (
            recModal.classList.contains('is-open') &&
            !recModal.contains(e.target) &&
            e.target !== recBtn
        ) {
            closeRecModal();
        }
    }, true);

    return {
        sync: syncRecModal,
        open: openRecModal,
        close: closeRecModal
    };
}
