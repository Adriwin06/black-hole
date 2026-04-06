import { $ } from '../../vendor.js';

var PRESENTATION_PRESET_MANIFEST_PATH = 'js/app/presentation/presets/manifest.json';

export var PRESENTATION_PRESETS = {};
export var PRESENTATION_PRESET_ORDER = [];

export var presentationPresetLoadState = {
    loading: false,
    loaded: false,
    error: null,
    promise: null
};

function clonePresentationPresetData(value) {
    return JSON.parse(JSON.stringify(value));
}

function requestPresentationJson(path) {
    return new Promise(function(resolve, reject) {
        if (typeof $ !== 'undefined' && $ && typeof $.getJSON === 'function') {
            $.getJSON(path)
                .done(function(data) { resolve(data); })
                .fail(function(jqXHR, textStatus, err) {
                    reject(new Error('Failed to load ' + path + ': ' + (err || textStatus || 'unknown error')));
                });
            return;
        }

        if (typeof XMLHttpRequest === 'undefined') {
            reject(new Error('No JSON loader available for ' + path));
            return;
        }

        var req = new XMLHttpRequest();
        req.open('GET', path, true);
        req.onreadystatechange = function() {
            if (req.readyState !== 4) return;
            if (req.status >= 200 && req.status < 300) {
                try {
                    resolve(JSON.parse(req.responseText));
                } catch (parseErr) {
                    reject(parseErr);
                }
                return;
            }
            reject(new Error('Failed to load ' + path + ': HTTP ' + req.status));
        };
        req.send();
    });
}

export function registerPresentationPreset(preset, fallbackName) {
    if (!preset || typeof preset !== 'object') return false;

    var name = (typeof preset.name === 'string' && preset.name.trim())
        ? preset.name.trim()
        : (fallbackName || '');
    if (!name) return false;

    var copy = clonePresentationPresetData(preset);
    copy.name = name;
    PRESENTATION_PRESETS[name] = copy;
    if (PRESENTATION_PRESET_ORDER.indexOf(name) === -1) {
        PRESENTATION_PRESET_ORDER.push(name);
    }
    return true;
}

export function ensurePresentationPresetsLoaded() {
    if (presentationPresetLoadState.loaded) {
        return Promise.resolve(PRESENTATION_PRESETS);
    }
    if (presentationPresetLoadState.promise) {
        return presentationPresetLoadState.promise;
    }

    presentationPresetLoadState.loading = true;
    presentationPresetLoadState.error = null;

    presentationPresetLoadState.promise = requestPresentationJson(PRESENTATION_PRESET_MANIFEST_PATH)
        .then(function(manifest) {
            if (!Array.isArray(manifest)) {
                throw new Error('Invalid presentation preset manifest format.');
            }

            PRESENTATION_PRESETS = {};
            PRESENTATION_PRESET_ORDER = [];

            var jobs = manifest.map(function(entry) {
                if (!entry || typeof entry.file !== 'string') {
                    return Promise.resolve(false);
                }
                return requestPresentationJson(entry.file)
                    .then(function(preset) {
                        var fallbackName = (typeof entry.name === 'string') ? entry.name : '';
                        if (!registerPresentationPreset(preset, fallbackName)) {
                            throw new Error('Invalid preset data in ' + entry.file);
                        }
                        return true;
                    });
            });
            return Promise.all(jobs);
        })
        .then(function() {
            presentationPresetLoadState.loading = false;
            presentationPresetLoadState.loaded = true;
            return PRESENTATION_PRESETS;
        })
        .catch(function(err) {
            presentationPresetLoadState.loading = false;
            presentationPresetLoadState.loaded = false;
            presentationPresetLoadState.error = err;
            console.warn('Presentation presets load failed:', err);
            return PRESENTATION_PRESETS;
        });

    return presentationPresetLoadState.promise;
}


