---
---

# Presentation Timeline Editor Guide

This guide explains how to use the in-app timeline editor.

If you want the raw JSON schema, see `docs/presentation-json.md`.

## 1. UI overview

The interface has three sliding panels:

| Button | Location | Panel contents |
|---|---|---|
| `ANIMATIONS` | left screen edge | Freefall Dive, Hover Approach, and `RECORD TO TIMELINE` controls for those live scenarios |
| `CONTROLS` | right screen edge | Physics and rendering controls (`dat.GUI`) |
| `TIMELINE` | bottom of screen | Preset dropdown, transport, dopesheet, annotations, recording modal, import/export |

Click an edge button to open its panel.

## 2. How to open the editor

1. Click `TIMELINE` at the bottom of the screen.
2. In the Preset dropdown, select `-- new empty --` to start from a blank draft.

## 3. Timeline panel layout

The bottom panel has three main columns:

- **Track list** on the left: parameter tracks plus annotation lanes.
- **Dopesheet** in the center: keyframes, annotation bars, ruler, and scrubber.
- **Inspector** on the right: key editing or annotation editing depending on the current selection.

The transport bar includes:

- `Play / Pause / Stop`
- `current time / duration` numeric inputs
- `Preset` dropdown
- `SAVE`
- `FX`
- `REC`
- `AUTO KEY`
- `+ TRACK`
- `TEXT`
- `IMPORT`
- `EXPORT`
- `Close`

## 4. Fast workflow

This is the fastest way to animate a normal presentation:

1. Open the Timeline panel.
2. Move the playhead to the first keyframe time, usually `0.00`.
3. Click `AUTO KEY`.
4. Change controls in the right-side `CONTROLS` panel.
5. Move the playhead to the next keyframe time.
6. Click `AUTO KEY` again.
7. Repeat for more beats.
8. Press `Play` to preview.

Changes auto-apply immediately. There is no separate global apply step for ordinary timeline editing.

### Capturing Freefall Dive / Hover Approach into the timeline

The built-in `Freefall Dive` and `Hover Approach` modes can be recorded directly into the timeline. This is the intended workflow when a hand-authored keyframe edit would be too awkward or too hard to match physically.

Workflow:

1. Open `TIMELINE` and move the playhead to the shot start.
2. Open `ANIMATIONS`.
3. In either `FREEFALL DIVE` or `HOVER APPROACH`, click `RECORD TO TIMELINE`.
   - Enable `Smooth recorded camera` first if you want a light smoothing pass on the captured camera path.
4. While the capture is running, move the camera normally:
   - left drag: orbit
   - right drag: pan
   - left + right drag: roll
5. Click `STOP & SAVE` in the same section.

That capture writes all of the following into the timeline at the current playhead:

- a `startDive` or `startHover` event
- a radius track: `dive.currentR` or `hover.currentR`
- `observerState.time`
- `camera.position.x/y/z`
- `camera.quaternion.x/y/z/w`
- `cameraPan.x/y`

The smoothing option only affects the recorded camera motion (`camera.position.*`, `camera.quaternion.*`, `cameraPan.*`). It does not alter the recorded dive/hover radius track or the captured `observerState.time` shader clock.

Playback then re-seeks the dive or hover radius from those captured samples, so the replay follows the recorded shot instead of trying to re-run the live motion with slightly different timing.

## 5. Motion Functions (`FX`)

Click `FX` in the transport bar to open the Motion Functions panel. It generates blocks of keyframes from a preset motion type.

| Type | What it animates | Tracks created |
|---|---|---|
| `Intro: sky reveal` | Starts close to the hole and widens into the full scene | camera transform tracks |
| `Orbit around BH` | Camera circles the black hole at constant radius and elevation | `camera.position.*`, `camera.quaternion.*` |
| `Zoom in / out` | Observer distance change | `observer.distance` |
| `Exposure fade` | Exposure ramp | `look.exposure` |
| `Inclination sweep` | Camera latitude sweep | `observer.orbital_inclination` |

Shared parameters:

- `Start time`
- `Duration`
- `Ease`

Orbit-specific parameters:

- `Number of orbits`
- `Direction`

Press `APPLY` in the Motion Functions panel to insert the generated keys.

## 6. What Auto Key does

- First click at `0.00` or on an empty draft: prompts for `Full initial state` or `Changes only`.
- First click elsewhere: captures a baseline only.
- Next clicks: compares the current runtime state to the previous baseline.
- For each changed value, the editor writes:
  - a key at the previous baseline time
  - a key at the current time
- If nothing changed, no tracks are added.

`AUTO KEY` is best for controls-panel animation, not for Freefall Dive / Hover shots. Use `RECORD TO TIMELINE` for those.

## 7. Manual track and keyframe editing

Use the Key Inspector on the right for direct editing.

Fields:

- `Path`
- `Time`
- `Ease`
- `Value`

Inspector buttons:

- `USE TIME`
- `LIVE VALUE`
- `SET KEY`
- `DELETE KEY`

Useful workflow:

1. Click a track row.
2. Move the playhead.
3. Press `K` to key that track's current live value at the playhead.

Use `TEXT` in the transport bar to create an annotation event at the current time. Selecting an annotation bar switches the inspector into annotation-edit mode.

## 8. Interpolation rules

- Number-to-number keys interpolate.
- Boolean and string values step; there is no interpolation between them.

Practical tip:

- For mode switches such as `accretion_mode`, `kerr_mode`, or `jet.mode`, use JSON `set` events rather than numeric tracks.

## 9. Keyboard shortcuts

| Shortcut | Action |
|---|---|
| `Space` | Play / Pause |
| `Delete` or `Backspace` | Delete selected keyframes |
| `Ctrl+A` | Select all keys on the active track |
| `Ctrl+Shift+A` | Select all keys across all tracks |
| `Shift+A` | Select all keys at the current playhead time |
| Double-click keyframe | Select all keys at that same time across tracks |
| `Ctrl+C` | Copy selected keyframes |
| `Ctrl+V` | Paste copied keyframes at the current playhead |
| `K` | Key the selected track using its live value |
| `Ctrl+Z` | Undo |
| `Ctrl+Y` or `Ctrl+Shift+Z` | Redo |
| `Escape` | Close the Motion Functions panel |
| `Home` | Seek to `0.00` |
| `End` | Seek to timeline end |
| `Left` / `Right` | Nudge playhead by `0.1 s` |
| `Shift+Left` / `Shift+Right` | Nudge playhead by `1.0 s` |

## 10. Undo / Redo

Every mutating operation pushes a snapshot onto the undo stack. Use `Ctrl+Z` to undo and `Ctrl+Y` or `Ctrl+Shift+Z` to redo.

## 11. Changes are live

Normal edits update the running timeline immediately. The only apply-style action is inside the Motion Functions panel, where generated keyframes are inserted into the draft.

## 12. Panel state persistence

Closing the Timeline panel stores its state in session storage:

- selected preset
- current draft
- selected track and key selection
- playhead position
- Auto Key baseline
- recording modal choices
- annotation / parameter-HUD visibility and HUD item list

Reopening the panel restores that state for the current tab/session.

## 13. JSON import / export

Transport-bar file actions:

- `SAVE`: re-download the current draft with its linked filename
- `IMPORT`: load a `.json` timeline file into the editor
- `EXPORT`: download the current draft as `<name>.json`

Use exported JSON to edit advanced fields such as mode-switch events, compile flags, custom annotation channels, manual annotation placement, or captured dive/hover event payloads.

## 14. Playback and recording

Playback controls, preset loading, and recording settings live in the `TIMELINE` panel, not in the `ANIMATIONS` panel.

Open `TIMELINE`, then click `REC` to access:

- loop toggle
- annotation and parameter-HUD visibility
- HUD parameter list with `ADD SELECTED` / `CLEAR`
- overlay inclusion in recordings
- reset-on-record-start
- recording quality, mode, resolution, FPS, bitrate
- `PNG SNAPSHOT`
- `START REC`
- `STOP REC`

The `ANIMATIONS` panel is for the live Freefall Dive / Hover Approach modes and for capturing those runs into timeline data.

## 15. Save as a reusable preset

`-- new empty --` is a runtime editing mode. Closing the panel preserves the current draft only for the current tab/session.

To make a preset permanent in the repository:

1. Click `EXPORT`.
2. Place the downloaded file in `js/app/presentation/presets/`.
3. Add an entry in `js/app/presentation/presets/manifest.json`.
4. Reload the page.

The preset will then appear in the Timeline-panel Preset dropdown.

## 16. Troubleshooting

**Timeline panel does not open**

- Click the `TIMELINE` button at the bottom of the screen.

**My edits do not appear during playback**

- Edits auto-apply. If playback seems stale, press `Stop` and `Play` again.

**AUTO KEY reports no changes**

- You likely pressed it twice without changing any controls between captures.

**A keyframe does nothing**

- The path may be wrong or unsupported.
- Type the path in the Inspector and click `LIVE VALUE` to verify it resolves.
- Check spelling against `docs/presentation-json.md`.

**Camera rotation was not captured by AUTO KEY**

- `AUTO KEY` captures camera transforms only when you move the camera between the two captures.

**I need to animate a dive or hover shot and hand-keying is too awkward**

- Use `RECORD TO TIMELINE` in the `ANIMATIONS` panel.
- That workflow captures the live radius evolution plus camera orbit/pan/roll into timeline tracks.
