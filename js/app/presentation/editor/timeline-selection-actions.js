import { clamp, clonePlain } from './timeline-utils.js';

export function createTimelineSelectionActions(options) {
    options = options || {};

    var getDraft = options.getDraft || function() { return null; };
    var getSelectedTrack = options.getSelectedTrack || function() { return ''; };
    var setSelectedTrack = options.setSelectedTrack || function() {};
    var setSelectedKeyT = options.setSelectedKeyT || function() {};
    var getSelectedKeys = options.getSelectedKeys || function() { return []; };
    var setSelectedKeys = options.setSelectedKeys || function() {};
    var pushUndo = options.pushUndo || function() {};
    var getTrackByPath = options.getTrackByPath || function() { return null; };
    var getKeyAt = options.getKeyAt || function() { return -1; };
    var applyDraft = options.applyDraft || function() {};
    var rebuildAll = options.rebuildAll || function() {};
    var rebuildLanes = options.rebuildLanes || function() {};
    var addToMultiSelect = options.addToMultiSelect || function() {};
    var clearMultiSelect = options.clearMultiSelect || function() {};
    var currentTime = options.currentTime || function() { return 0; };
    var getDuration = options.getDuration || function() { return 0; };
    var setStatus = options.setStatus || function() {};
    var inspSummary = options.inspSummary || null;

    var clipboard = null;

    function deleteTrack(path) {
        var draft = getDraft();
        if (!draft || !path) return;
        pushUndo();
        draft.tracks = draft.tracks.filter(function(track) {
            return track.path !== path;
        });
        setSelectedKeys(getSelectedKeys().filter(function(key) {
            return key.path !== path;
        }));
        if (getSelectedTrack() === path) {
            setSelectedTrack(draft.tracks.length ? draft.tracks[0].path : '');
            setSelectedKeyT(NaN);
        }
        applyDraft();
        rebuildAll();
        setStatus('Track deleted: ' + path, '');
    }

    function deleteSelectedKeys() {
        var draft = getDraft();
        var selectedTrack = getSelectedTrack();
        var selectedKeys = getSelectedKeys();
        if (!draft) return;
        if (!selectedKeys.length && selectedTrack) {
            deleteTrack(selectedTrack);
            return;
        }
        if (!selectedKeys.length) return;
        pushUndo();
        var count = 0;
        for (var i = 0; i < selectedKeys.length; i++) {
            var selectedKey = selectedKeys[i];
            var track = getTrackByPath(selectedKey.path);
            if (!track) continue;
            var keyIndex = getKeyAt(track, selectedKey.t);
            if (keyIndex >= 0) {
                track.keys.splice(keyIndex, 1);
                count++;
            }
        }
        draft.tracks = draft.tracks.filter(function(track) {
            return track.keys.length > 0;
        });
        setSelectedKeyT(NaN);
        clearMultiSelect();
        applyDraft();
        rebuildAll();
        if (count) {
            setStatus(count + ' key' + (count > 1 ? 's' : '') + ' deleted.', '');
        }
    }

    function selectAllKeysOnTrack() {
        var draft = getDraft();
        var selectedTrack = getSelectedTrack();
        if (!draft || !selectedTrack) return;
        var track = getTrackByPath(selectedTrack);
        if (!track) return;
        clearMultiSelect();
        for (var i = 0; i < track.keys.length; i++) {
            addToMultiSelect(track.path, track.keys[i].t);
        }
        rebuildLanes();
        if (inspSummary) {
            inspSummary.textContent = getSelectedKeys().length + ' keyframes selected on ' + selectedTrack;
        }
    }

    function selectAllKeys() {
        var draft = getDraft();
        if (!draft) return;
        clearMultiSelect();
        for (var i = 0; i < draft.tracks.length; i++) {
            var track = draft.tracks[i];
            for (var k = 0; k < track.keys.length; k++) {
                addToMultiSelect(track.path, track.keys[k].t);
            }
        }
        rebuildLanes();
        if (inspSummary) {
            inspSummary.textContent = getSelectedKeys().length + ' keyframes selected (all).';
        }
    }

    function copySelectedKeys() {
        var draft = getDraft();
        var selectedKeys = getSelectedKeys();
        if (!draft || !selectedKeys.length) return;
        var anchorT = Infinity;
        for (var i = 0; i < selectedKeys.length; i++) {
            if (selectedKeys[i].t < anchorT) anchorT = selectedKeys[i].t;
        }
        var entries = [];
        for (var s = 0; s < selectedKeys.length; s++) {
            var selectedKey = selectedKeys[s];
            var track = getTrackByPath(selectedKey.path);
            if (!track) continue;
            var keyIndex = getKeyAt(track, selectedKey.t);
            if (keyIndex < 0) continue;
            var key = track.keys[keyIndex];
            entries.push({
                path: selectedKey.path,
                relT: selectedKey.t - anchorT,
                v: clonePlain(key.v),
                ease: key.ease
            });
        }
        if (!entries.length) return;
        clipboard = { anchorT: anchorT, entries: entries };
        setStatus(entries.length + ' key' + (entries.length > 1 ? 's' : '') + ' copied.', '');
    }

    function pasteKeys() {
        var draft = getDraft();
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

    return {
        deleteTrack: deleteTrack,
        deleteSelectedKeys: deleteSelectedKeys,
        selectAllKeysOnTrack: selectAllKeysOnTrack,
        selectAllKeys: selectAllKeys,
        copySelectedKeys: copySelectedKeys,
        pasteKeys: pasteKeys
    };
}
