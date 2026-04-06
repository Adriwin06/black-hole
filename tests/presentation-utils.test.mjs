import test from 'node:test';
import assert from 'node:assert/strict';

import {
    presentationEasing,
    samplePresentationTrack
} from '../js/app/presentation/runtime/track-sampling.js';
import {
    normalizeTimelineDraft,
    normalizeTimelineParamHudConfig
} from '../js/app/presentation/editor/timeline-utils.js';

test('presentationEasing preserves endpoints and smooths midpoints', function() {
    assert.equal(presentationEasing(0, 'smooth'), 0);
    assert.equal(presentationEasing(1, 'smooth'), 1);
    assert.equal(presentationEasing(0, 'smoother'), 0);
    assert.equal(presentationEasing(1, 'smoother'), 1);
    assert.equal(presentationEasing(0.5, 'linear'), 0.5);
    assert.ok(presentationEasing(0.5, 'smooth') > 0.49);
    assert.ok(presentationEasing(0.5, 'smooth') < 0.51);
});

test('samplePresentationTrack interpolates numeric keys and clamps to edges', function() {
    var track = {
        keys: [
            { t: 0, v: 10 },
            { t: 2, v: 30, ease: 'smooth' }
        ]
    };

    assert.equal(samplePresentationTrack(track, -1), 10);
    assert.equal(samplePresentationTrack(track, 3), 30);
    assert.equal(samplePresentationTrack(track, 1), 20);
});

test('samplePresentationTrack keeps discrete values until the next key', function() {
    var track = {
        keys: [
            { t: 0, v: 'idle' },
            { t: 1, v: 'active' }
        ]
    };

    assert.equal(samplePresentationTrack(track, 0.5), 'idle');
    assert.equal(samplePresentationTrack(track, 1), 'active');
});

test('normalizeTimelineParamHudConfig clamps values and deduplicates items', function() {
    var config = normalizeTimelineParamHudConfig({
        anchorX: 2,
        anchorY: -1,
        fontSize: 200,
        items: [
            { path: ' params.look.glow ', label: ' Glow ' },
            { path: 'params.look.glow', label: 'Duplicate' },
            { path: 'params.look.disk_gain' },
            null
        ]
    }, null);

    assert.equal(config.anchorX, 1);
    assert.equal(config.anchorY, 0);
    assert.equal(config.fontSize, 48);
    assert.deepEqual(config.items, [
        { path: 'params.look.glow', label: 'Glow' },
        { path: 'params.look.disk_gain', label: 'params.look.disk_gain' }
    ]);
});

test('normalizeTimelineDraft sorts tracks, keys, and events while preserving defaults', function() {
    var draft = normalizeTimelineDraft({
        name: 'Demo',
        duration: 0.1,
        tracks: [
            {
                path: ' params.b ',
                keys: [
                    { t: 2, v: 2, ease: 'bogus' },
                    { t: 1, v: 1, ease: 'smooth' }
                ]
            },
            {
                path: 'params.a',
                compile: true,
                keys: [
                    { t: 'bad', v: 0 },
                    { t: 3, v: 9, ease: 'smoother' }
                ]
            }
        ],
        events: [
            { t: 5, action: 'later' },
            { t: 1, action: 'first' },
            { t: 'bad', action: 'drop' }
        ],
        annotationTracks: [],
        annotations: { includeInRecording: true },
        paramHud: {
            enabled: false,
            items: [{ path: 'params.look.glow' }]
        }
    }, {
        annotationsState: { enabled: true, includeInRecording: false },
        paramHudState: { enabled: true, includeInRecording: false, anchorX: 0, anchorY: 1, fontSize: 11, items: [] }
    });

    assert.equal(draft.name, 'Demo');
    assert.equal(draft.duration, 0.5);
    assert.deepEqual(draft.tracks.map(function(track) { return track.path; }), ['params.a', 'params.b']);
    assert.deepEqual(draft.tracks[1].keys, [
        { t: 1, v: 1, ease: 'smooth' },
        { t: 2, v: 2, ease: 'linear' }
    ]);
    assert.deepEqual(draft.events.map(function(event) { return event.action; }), ['first', 'later']);
    assert.equal(draft.annotationTracks.length, 1);
    assert.equal(draft.annotationTracks[0].label, 'Annotation 1');
    assert.equal(draft.annotations.enabled, true);
    assert.equal(draft.annotations.includeInRecording, true);
    assert.equal(draft.paramHud.enabled, false);
    assert.deepEqual(draft.paramHud.items, [
        { path: 'params.look.glow', label: 'params.look.glow' }
    ]);
});
