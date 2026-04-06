export function setupTimelineFileIo(options) {
    options = options || {};

    var elements = options.elements || {};
    var importBtn = elements.importBtn;
    var exportBtn = elements.exportBtn;
    var saveBtn = elements.saveBtn;
    var delPresetBtn = elements.delPresetBtn;
    var presetSelect = elements.presetSelect;

    var getDraft = options.getDraft || function() { return null; };
    var setDraftName = options.setDraftName || function() {};
    var getImportedPresets = options.getImportedPresets || function() { return {}; };
    var getLinkedFileName = options.getLinkedFileName || function() { return null; };
    var setLinkedFileName = options.setLinkedFileName || function() {};
    var getApplyingDraft = options.getApplyingDraft || function() { return false; };
    var setApplyingDraft = options.setApplyingDraft || function() {};
    var isPanelOpen = options.isPanelOpen || function() { return false; };
    var getPresentationTimeline = options.getPresentationTimeline || function() { return null; };
    var setPresentationTimeline = options.setPresentationTimeline || function() { return false; };
    var registerPresentationPreset = options.registerPresentationPreset || function() {};
    var loadPresetByName = options.loadPresetByName || function() {};
    var populatePresets = options.populatePresets || function() {};
    var updateDelPresetBtn = options.updateDelPresetBtn || function() {};
    var syncFromRuntime = options.syncFromRuntime || function() {};
    var setStatus = options.setStatus || function() {};
    var esc = options.esc || function(value) { return value; };
    var clonePlain = options.clonePlain || function(value) { return value; };
    var presentationPresets = options.presentationPresets || null;
    var presentationPresetOrder = options.presentationPresetOrder || null;

    function getTimelineForExport() {
        var timeline = getPresentationTimeline();
        if (!timeline) {
            var draft = getDraft();
            if (draft) timeline = clonePlain(draft);
        }
        return timeline;
    }

    if (importBtn) {
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
                        var presetName = (obj && obj.name && obj.name.trim()) ? obj.name.trim()
                            : file.name.replace(/\.json$/i, '').replace(/[_-]+/g, ' ').trim();
                        if (!presetName) presetName = 'Imported';
                        obj.name = presetName;

                        registerPresentationPreset(obj, presetName);
                        getImportedPresets()[presetName] = true;

                        setApplyingDraft(true);
                        if (setPresentationTimeline(obj)) {
                            setApplyingDraft(false);
                            syncFromRuntime();
                            setLinkedFileName(file.name);
                            populatePresets();
                            if (presetSelect) presetSelect.value = presetName;
                            updateDelPresetBtn();
                            setStatus('Imported: ' + file.name, '');
                        } else {
                            setApplyingDraft(false);
                            setStatus('Invalid timeline data in ' + file.name, 'tl-status--error');
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
    }

    if (exportBtn) {
        exportBtn.addEventListener('click', function() {
            var timeline = getTimelineForExport();
            if (!timeline) {
                setStatus('No timeline to export.', 'tl-status--warn');
                return;
            }
            var defaultName = (timeline.name && timeline.name !== 'Untitled' && timeline.name !== 'Custom')
                ? timeline.name
                : '';
            var exportName = prompt('Name for the exported timeline:', defaultName);
            if (exportName === null) return;
            exportName = exportName.trim();
            if (!exportName) exportName = 'timeline';
            timeline.name = exportName;
            setDraftName(exportName);

            var json = JSON.stringify(timeline, null, 2);
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

            setLinkedFileName(safeFilename + '.json');
            registerPresentationPreset(clonePlain(timeline), exportName);
            getImportedPresets()[exportName] = true;
            populatePresets();
            if (presetSelect) presetSelect.value = exportName;
            updateDelPresetBtn();
            setStatus('Exported: ' + a.download, '');
        });
    }

    if (delPresetBtn) {
        delPresetBtn.addEventListener('click', function() {
            var name = presetSelect ? presetSelect.value : '';
            if (!name || !getImportedPresets()[name]) return;
            if (!confirm('Remove imported preset "' + esc(name) + '" from the list?')) return;
            delete getImportedPresets()[name];
            if (presentationPresets) delete presentationPresets[name];
            if (presentationPresetOrder) {
                var idx = presentationPresetOrder.indexOf(name);
                if (idx !== -1) presentationPresetOrder.splice(idx, 1);
            }
            if (presetSelect) presetSelect.value = '';
            loadPresetByName('');
            populatePresets();
            setStatus('Preset "' + name + '" removed.', '');
        });
    }

    if (saveBtn) {
        saveBtn.addEventListener('click', function() {
            var timeline = getTimelineForExport();
            if (!timeline) {
                setStatus('No timeline to save.', 'tl-status--warn');
                return;
            }

            var linkedFileName = getLinkedFileName();
            if (!linkedFileName) {
                var suggestedName = (timeline.name && timeline.name !== 'Untitled' && timeline.name !== 'Custom')
                    ? timeline.name
                    : '';
                var saveName = prompt('No file linked. Enter a name to save as:', suggestedName);
                if (saveName === null) return;
                saveName = saveName.trim();
                if (!saveName) {
                    setStatus('Save cancelled.', 'tl-status--warn');
                    return;
                }
                timeline.name = saveName;
                setDraftName(saveName);
                linkedFileName = saveName.replace(/[^a-zA-Z0-9_-]/g, '_') + '.json';
                setLinkedFileName(linkedFileName);
            }

            var json = JSON.stringify(timeline, null, 2);
            var blob = new Blob([json], { type: 'application/json' });
            var url = URL.createObjectURL(blob);
            var a = document.createElement('a');
            a.href = url;
            a.download = linkedFileName;
            document.body.appendChild(a);
            a.click();
            a.remove();
            URL.revokeObjectURL(url);

            var presetName = timeline.name || '';
            if (presetName) {
                registerPresentationPreset(timeline, presetName);
                getImportedPresets()[presetName] = true;
                populatePresets();
                if (presetSelect) presetSelect.value = presetName;
                updateDelPresetBtn();
            }
            setStatus('Saved: ' + a.download, '');
        });
    }

    window.addEventListener('presentation:timeline-panel-sync', function() {
        if (isPanelOpen() && !getApplyingDraft()) syncFromRuntime();
    });
}
