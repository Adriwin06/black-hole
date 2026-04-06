import { clamp } from './timeline-utils.js';

export function attachTimelineKeyboardShortcuts(options) {
    options = options || {};

    var panel = options.panel || null;
    var motionModal = options.motionModal || null;
    var inspPath = options.inspPath || null;
    var inspTime = options.inspTime || null;

    var isPanelOpen = typeof options.isPanelOpen === 'function'
        ? options.isPanelOpen
        : function() { return false; };
    var closeMotionModal = typeof options.closeMotionModal === 'function'
        ? options.closeMotionModal
        : function() {};
    var getSelectedEventIndex = typeof options.getSelectedEventIndex === 'function'
        ? options.getSelectedEventIndex
        : function() { return -1; };
    var setSelectedEventIndex = typeof options.setSelectedEventIndex === 'function'
        ? options.setSelectedEventIndex
        : function() {};
    var showKeyInspector = typeof options.showKeyInspector === 'function'
        ? options.showKeyInspector
        : function() {};
    var rebuildEventsLane = typeof options.rebuildEventsLane === 'function'
        ? options.rebuildEventsLane
        : function() {};
    var getPresentationState = typeof options.getPresentationState === 'function'
        ? options.getPresentationState
        : function() { return null; };
    var pausePresentation = typeof options.pausePresentation === 'function'
        ? options.pausePresentation
        : function() {};
    var playPresentation = typeof options.playPresentation === 'function'
        ? options.playPresentation
        : function() {};
    var getSelectedTrack = typeof options.getSelectedTrack === 'function'
        ? options.getSelectedTrack
        : function() { return ''; };
    var currentTime = typeof options.currentTime === 'function'
        ? options.currentTime
        : function() { return 0; };
    var getDuration = typeof options.getDuration === 'function'
        ? options.getDuration
        : function() { return 0; };
    var setStatus = typeof options.setStatus === 'function'
        ? options.setStatus
        : function() {};
    var captureInspectorValue = typeof options.captureInspectorValue === 'function'
        ? options.captureInspectorValue
        : function() {};
    var doSetKey = typeof options.doSetKey === 'function'
        ? options.doSetKey
        : function() {};
    var deleteSelectedKeys = typeof options.deleteSelectedKeys === 'function'
        ? options.deleteSelectedKeys
        : function() {};
    var selectAllKeysAtTime = typeof options.selectAllKeysAtTime === 'function'
        ? options.selectAllKeysAtTime
        : function() {};
    var selectAllKeysOnTrack = typeof options.selectAllKeysOnTrack === 'function'
        ? options.selectAllKeysOnTrack
        : function() {};
    var selectAllKeys = typeof options.selectAllKeys === 'function'
        ? options.selectAllKeys
        : function() {};
    var copySelectedKeys = typeof options.copySelectedKeys === 'function'
        ? options.copySelectedKeys
        : function() {};
    var pasteKeys = typeof options.pasteKeys === 'function'
        ? options.pasteKeys
        : function() {};
    var seekPresentation = typeof options.seekPresentation === 'function'
        ? options.seekPresentation
        : function() {};
    var updateTimeInputs = typeof options.updateTimeInputs === 'function'
        ? options.updateTimeInputs
        : function() {};
    var updateScrubber = typeof options.updateScrubber === 'function'
        ? options.updateScrubber
        : function() {};
    var updatePlayheads = typeof options.updatePlayheads === 'function'
        ? options.updatePlayheads
        : function() {};
    var undo = typeof options.undo === 'function'
        ? options.undo
        : function() {};
    var redo = typeof options.redo === 'function'
        ? options.redo
        : function() {};

    function handleKeyDown(e) {
        if (!isPanelOpen()) return;

        var target = e.target || null;
        var tag = ((target && target.tagName) || '').toLowerCase();
        var inInput = (tag === 'input' || tag === 'textarea' || tag === 'select' ||
            !!(target && target.isContentEditable));

        var inExternalInput = inInput && panel && !panel.contains(target);
        if ((e.ctrlKey || e.metaKey) && !e.shiftKey && e.code === 'KeyZ' && !inExternalInput) {
            e.preventDefault();
            e.stopPropagation();
            undo();
            return;
        }
        if ((e.ctrlKey || e.metaKey) &&
            (e.code === 'KeyY' || (e.shiftKey && e.code === 'KeyZ')) &&
            !inExternalInput) {
            e.preventDefault();
            e.stopPropagation();
            redo();
            return;
        }

        if (inInput) return;

        if (e.code === 'Escape') {
            if (motionModal && motionModal.classList.contains('is-open')) {
                closeMotionModal();
                e.preventDefault();
                return;
            }
            if (getSelectedEventIndex() >= 0) {
                setSelectedEventIndex(-1);
                showKeyInspector();
                rebuildEventsLane();
                e.preventDefault();
            }
            return;
        }

        if (e.code === 'Space') {
            e.preventDefault();
            var state = getPresentationState();
            if (state && state.active && !state.paused) {
                pausePresentation();
            } else {
                playPresentation(false);
            }
            return;
        }

        if (!e.ctrlKey && !e.metaKey && !e.altKey && e.code === 'KeyK') {
            var keyPath = getSelectedTrack() || (inspPath && inspPath.value ? inspPath.value.trim() : '');
            if (!keyPath) {
                setStatus('Select a track first, then press K to key its live value.', 'tl-status--warn');
                e.preventDefault();
                return;
            }
            if (inspPath) inspPath.value = keyPath;
            if (inspTime) inspTime.value = currentTime().toFixed(2);
            captureInspectorValue(keyPath);
            doSetKey();
            e.preventDefault();
            return;
        }

        if (e.code === 'Delete' || e.code === 'Backspace') {
            e.preventDefault();
            deleteSelectedKeys();
            return;
        }

        if (e.shiftKey && e.code === 'KeyA') {
            e.preventDefault();
            selectAllKeysAtTime(currentTime());
            return;
        }

        if ((e.ctrlKey || e.metaKey) && e.code === 'KeyA') {
            e.preventDefault();
            if (getSelectedTrack() && !e.shiftKey) selectAllKeysOnTrack();
            else selectAllKeys();
            return;
        }

        if ((e.ctrlKey || e.metaKey) && !e.shiftKey && e.code === 'KeyC') {
            e.preventDefault();
            copySelectedKeys();
            return;
        }

        if ((e.ctrlKey || e.metaKey) && !e.shiftKey && e.code === 'KeyV') {
            e.preventDefault();
            pasteKeys();
            return;
        }

        if (e.code === 'Home') {
            e.preventDefault();
            seekPresentation(0);
            updateTimeInputs();
            updateScrubber();
            updatePlayheads();
            return;
        }

        if (e.code === 'End') {
            e.preventDefault();
            seekPresentation(getDuration());
            updateTimeInputs();
            updateScrubber();
            updatePlayheads();
            return;
        }

        if (e.code === 'ArrowLeft' || e.code === 'ArrowRight') {
            var step = e.shiftKey ? 1.0 : 0.1;
            var dir = (e.code === 'ArrowLeft') ? -1 : 1;
            var newT = clamp(currentTime() + dir * step, 0, getDuration());
            seekPresentation(newT);
            updateTimeInputs();
            updateScrubber();
            updatePlayheads();
            e.preventDefault();
        }
    }

    window.addEventListener('keydown', handleKeyDown, true);
    return function detachTimelineKeyboardShortcuts() {
        window.removeEventListener('keydown', handleKeyDown, true);
    };
}
