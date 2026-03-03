# Presentation JSON Guide

This document explains the JSON format used by the presentation timeline system.

- Controller: `js/app/presentation/presentation-controller.js`
- Presets folder: `js/app/presentation/presets/`
- Preset manifest: `js/app/presentation/presets/manifest.json`
- UI workflow guide: `docs/presentation-editor.md`

## 1. Minimal Preset

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

## 2. Root Fields

- `name` (string): preset name shown in UI.
- `duration` (number, seconds): total timeline duration.
- `loop` (boolean): default loop flag when loaded.
- `tracks` (array): keyframed values sampled continuously.
- `events` (array): discrete actions executed at specific times.

Notes:
- If `duration` is missing or invalid, it is inferred from the largest `t` in tracks/events.
- Negative `t` values are clamped to `0`.
- Tracks and events are sorted by `t` at load time.

## 3. Track Format

Each track object:

- `path` (string): property path to animate.
- `compile` (optional boolean): force shader recompile when this track changes.
- `keys` (array): keyframes.

Each keyframe:

- `t` (number): time in seconds.
- `v` (any): value.
- `ease` (optional string): `linear` (default), `smooth`, or `smoother`.

Interpolation behavior:
- Number to number: interpolated.
- Non-number types (bool/string): stepped value (`a.v` until end of segment, then `b.v`).

## 4. Event Format

Each event object:

- `t` (number): time in seconds.
- `action` (string): one of the actions below.
- `path` / `value` / `compile` / `note`: only for actions that use them.

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

- `set`: requires `path` and `value`, optional `compile`.
- `updateShader`: no extra fields.
- dive/hover actions: no extra fields.
- `annotation`: requires `note` object.
- `clearAnnotation`: no extra fields.

## 5. Path Resolution (`path`)

Path prefixes supported by the controller:

- `cameraPan.*` -> camera pan vector (`cameraPan.x`, `cameraPan.y`)
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
- `black_hole.spin_enabled`
- `accretion_mode`
- `kerr_mode`
- `jet.enabled`
- `cameraPan.x`
- `cameraPan.y`

Important:
- If a path does not exist, it is ignored.
- Numeric targets parse with `parseFloat`.
- Boolean targets coerce with `!!value`.
- String enum targets should use exact expected strings.

## 6. Shader Recompile Rules

Some path changes require recompile. This can happen in 2 ways:

- Automatically for known compile-sensitive paths (for example `kerr_mode`, `accretion_mode`, `jet.enabled`, `grmhd.enabled`, `observer.motion`, `taa_enabled`, `rk4_integration`, etc.).
- Explicitly by setting `compile: true` on a track or `set` event.

You can also insert an explicit event:

```json
{ "t": 10.02, "action": "updateShader" }
```

## 7. Annotation Notes (`action: "annotation"`)

`note` supports:

- `title` (string)
- `text` (string)
- `body` (string, fallback if `text` omitted)
- `color` (hex `#RRGGBB`)
- `width` (number, clamped by viewport)
- `placement` (`auto`, `left`, `right`, `top`, `bottom`)
- `offset` (number, distance from anchor)
- `anchor` (object)

Anchor modes:

- World anchor:

```json
"anchor": { "mode": "world", "target": "black_hole" }
```

Supported world targets:
- `black_hole` (`bh`, `center` aliases)
- `disk`
- `jet_north`
- `jet_south`
- `planet`

Custom world coordinate:

```json
"anchor": { "mode": "world", "x": 0, "y": 0, "z": 0 }
```

Screen anchor (normalized 0..1):

```json
"anchor": { "x": 0.52, "y": 0.48 }
```

Placement behavior:

- `auto` chooses left/right based on available safe area (accounts for UI panels).
- To force position, set `placement` explicitly.
- Fine-tune with `offset`.

Example forced placement:

```json
{
  "t": 2.5,
  "action": "annotation",
  "note": {
    "title": "Photon Ring",
    "text": "Strong lensing around the shadow edge.",
    "anchor": { "mode": "world", "target": "black_hole" },
    "placement": "left",
    "offset": 56,
    "color": "#7cc5ff",
    "width": 320
  }
}
```

## 8. Playback Semantics

- `t=0` events are triggered when playback starts from the beginning.
- During seek, timeline state is rebuilt from start to seek time.
- `play(true)` or replay at end reinitializes camera and resets dive/hover states.
- If `loop` is true, events replay correctly each cycle.

## 9. Add a New Preset

1. Create a JSON file in `js/app/presentation/presets/`.
2. Add it to `js/app/presentation/presets/manifest.json`:

```json
{
  "name": "My Preset",
  "file": "js/app/presentation/presets/my-preset.json"
}
```

3. Reload the page and select it in `ANIMATIONS > PRESENTATION TIMELINE`.

## 10. Practical Tips

- For mode toggles (`accretion_mode`, `kerr_mode`, `jet.mode`), prefer `set` events + `updateShader`.
- Keep boolean/string values in events (not interpolated tracks) unless you want step behavior.
- If annotation points to the wrong region in dynamic shots, switch from symbolic target to explicit screen anchor (`anchor.x`, `anchor.y`).
- Use small event offsets like `+0.02s` after mode switches to avoid race-like visuals during heavy updates.

