---
---

# Presentation JSON Guide

This document explains the JSON format used by the presentation timeline system.

- Controller: `js/app/presentation/presentation-controller.js`
- Presets folder: `js/app/presentation/presets/`
- Preset manifest: `js/app/presentation/presets/manifest.json`
- UI workflow guide: `docs/presentation-editor.md`

## 1. Minimal preset

```json
{
  "name": "My Preset",
  "duration": 12,
  "loop": false,
  "tracks": [
    {
      "path": "observer.distance",
      "keys": [
        { "t": 0, "v": 18 },
        { "t": 12, "v": 10, "ease": "smooth" }
      ]
    }
  ],
  "events": [
    { "t": 0, "action": "set", "path": "observer.motion", "value": false },
    { "t": 0.01, "action": "updateShader" }
  ]
}
```

## 2. Root fields

- `name` (string): preset name shown in the UI
- `duration` (number, seconds): total timeline duration
- `loop` (boolean): default loop flag when loaded
- `tracks` (array): continuously sampled keyframed values
- `events` (array): discrete actions executed at exact times
- `annotationTracks` (array, optional): annotation-lane metadata
- `annotations` (object, optional): annotation-overlay defaults
- `paramHud` (object, optional): parameter-HUD defaults and visible item list

Notes:

- If `duration` is missing or invalid, it is inferred from the largest `t` in tracks and events.
- Negative times are clamped to `0`.
- Keys and events are sorted by `t` on load.
- If `annotationTracks` is omitted, a default lane is created.
- If `annotations` or `paramHud` is omitted, the current runtime/editor overlay settings are kept.

`paramHud` supports:

- `enabled` (boolean)
- `includeInRecording` (boolean)
- `anchorX` / `anchorY` (number, normalized `0..1`)
- `fontSize` (number, pixels)
- `items` (array of `{ "path": "...", "label": "..." }`)

## 3. Track format

Each track object contains:

- `path` (string): property path to animate
- `compile` (optional boolean): force shader recompile when this track changes
- `keys` (array): keyframes

Each keyframe contains:

- `t` (number): time in seconds
- `v` (any): value
- `ease` (optional string): `linear`, `smooth`, or `smoother`

Interpolation rules:

- number to number: interpolated
- non-number values: stepped

Important:

- Quality-owned paths such as `quality`, `n_steps`, `sample_count`, `resolution_scale`, and the TAA tuning fields are ignored on load in both tracks and `set` events, so timelines cannot silently override the active render preset.

## 4. Event format

Each event object contains:

- `t` (number): time in seconds
- `action` (string): one of the supported actions below
- optional action-specific fields such as `path`, `value`, `compile`, `note`, `channel`, `position`, or `velocity`

Supported actions:

- `set`
- `updateShader`
- `startDive`
- `pauseDive`
- `resetDive`
- `startHover`
- `pauseHover`
- `resetHover`
- `annotation`
- `clearAnnotation`

Action fields:

- `set`: requires `path` and `value`, optional `compile`
- `updateShader`: no extra fields
- `pauseDive`, `resetDive`, `pauseHover`, `resetHover`: no extra fields
- `annotation`: requires `note`, optional `channel`
- `clearAnnotation`: optional `channel`

### Captured dive / hover start events

`startDive` and `startHover` can carry extra fields. These are produced automatically by the `RECORD TO TIMELINE` workflow in the `ANIMATIONS` panel and are preserved by JSON export/import.

Optional fields:

- `position`: `{ "x": number, "y": number, "z": number }`
- `velocity`: `{ "x": number, "y": number, "z": number }`
- `prevMotionState`: boolean
- `prevDistance`: number
- `observerTime`: number

These fields let playback restart the live dive/hover mode from the same observer state that was used during capture before the radius track takes over.

## 5. Path resolution

Supported path prefixes:

- `cameraPan.*` -> camera pan vector
- `camera.*` -> live Three.js camera transform
- `observerState.*` -> runtime observer object
- `dive.*` -> dive state object
- `hover.*` -> hover state object
- `params.*` -> `shader.parameters.*`
- `shader.parameters.*` -> `shader.parameters.*`
- no prefix -> treated as `shader.parameters.*`

Examples:

- `observer.distance`
- `observer.orbital_inclination`
- `observer.motion`
- `look.exposure`
- `black_hole.spin`
- `accretion_mode`
- `jet.enabled`
- `camera.position.x`
- `camera.quaternion.w`
- `cameraPan.x`
- `cameraPan.y`
- `dive.currentR`
- `hover.currentR`

Important:

- If a path does not exist, it is ignored.
- Numeric targets are parsed with `parseFloat`.
- Boolean targets coerce with `!!value`.
- String enum targets should use the exact expected string.

### Captured radius tracks

`dive.currentR` and `hover.currentR` are special timeline paths used by the live capture workflow.

- `dive.currentR` only has meaning while dive mode is active.
- `hover.currentR` only has meaning while hover mode is active.

When those tracks are present and the matching mode is running, playback does not re-integrate the live solver for that segment. Instead, it seeks the dive or hover radius to the sampled value for deterministic replay.

In practice, you should usually create these tracks by recording from the `ANIMATIONS` panel rather than authoring them by hand.

## 6. Shader recompile rules

Some path changes require a shader recompile. This happens in two ways:

- automatically for known compile-sensitive paths such as `kerr_mode`, `accretion_mode`, `jet.enabled`, `grmhd.enabled`, `observer.motion`, `taa_enabled`, and `rk4_integration`
- explicitly by setting `compile: true` on a track or a `set` event

You can also insert a dedicated event:

```json
{ "t": 10.02, "action": "updateShader" }
```

## 7. Annotation notes (`action: "annotation"`)

`note` supports:

- `title` (string)
- `text` (string)
- `body` (string, fallback if `text` is omitted)
- `color` (hex `#RRGGBB`)
- `width` (number)
- `fadeIn` (number, seconds)
- `placement` (`auto`, `left`, `right`, `top`, `bottom`)
- `offset` (number)
- `boxX` / `boxY` (number, normalized `0..1`)
- `anchor` (object)

Channel support:

- use `annotationTracks` at the preset root to label annotation lanes
- use `channel` on `annotation` and `clearAnnotation` to target a specific lane

Anchor modes:

World anchor:

```json
"anchor": { "mode": "world", "target": "black_hole" }
```

Supported world targets:

- `black_hole`
- `bh`
- `center`
- `disk`
- `jet_north`
- `jet_south`
- `planet`

Custom world coordinate:

```json
"anchor": { "mode": "world", "x": 0, "y": 0, "z": 0 }
```

Screen anchor:

```json
"anchor": { "x": 0.52, "y": 0.48 }
```

## 8. Playback semantics

- `t=0` events fire when playback starts from the beginning.
- During seek, timeline state is rebuilt from start to the seek time.
- `play(true)` or replay from the end reinitializes the camera and resets dive/hover state.
- If `loop` is true, events replay correctly each cycle.
- Captured dive/hover shots replay deterministically because `dive.currentR` and `hover.currentR` are applied as radius seeks, not as free-running live integration.

## 9. Add a new preset

1. Create a JSON file in `js/app/presentation/presets/`.
2. Add it to `js/app/presentation/presets/manifest.json`:

```json
{
  "name": "My Preset",
  "file": "js/app/presentation/presets/my-preset.json"
}
```

3. Reload the page and select it in the `TIMELINE` Preset dropdown.

## 10. Practical tips

- For mode toggles such as `accretion_mode`, `kerr_mode`, or `jet.mode`, prefer `set` events plus `updateShader`.
- Keep boolean and string changes in events unless you intentionally want step behavior on a track.
- The `TEXT` workflow in the timeline editor covers common annotation editing; use raw JSON for advanced note fields and custom channels.
- If an annotation points to the wrong place in a dynamic shot, switch from a symbolic world target to an explicit screen anchor.
- Use the `ANIMATIONS` panel's `RECORD TO TIMELINE` flow for Freefall Dive or Hover shots instead of hand-authoring `startDive` / `startHover` plus radius tracks.
