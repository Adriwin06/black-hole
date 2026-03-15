---
---

# Presentation Timeline Editor Guide

This guide explains how to use the in-app timeline editor.

If you want the raw JSON schema, see `docs/presentation-json.md`.

## 1. UI overview — three panels

The interface has three separate sliding panels:

| Button | Location | Panel contents |
|---|---|---|
| `◀ ANIMATIONS` | left screen edge | Freefall Dive, Hover Approach, **Presentation Timeline** preset picker + status |
| `▶ CONTROLS` | right screen edge | Physics / rendering controls (dat.GUI) |
| `▲ TIMELINE` | bottom of screen | **Dopesheet editor** — tracks, keyframes, annotations, transport, and the REC modal |

Click an edge button to open its panel. The `▲ TIMELINE` panel also opens automatically when you select **New Preset** in the Animations panel.

## 2. How to open the editor

**Option A — via the Animations panel:**

1. Click `◀ ANIMATIONS` on the left edge to open the Animations panel.
2. Expand the **PRESENTATION TIMELINE** section.
3. In the `Preset` dropdown, select **`New Preset`**.
   - The `▲ TIMELINE` bottom panel opens automatically with a blank draft.

**Option B — directly:**

1. Click `▲ TIMELINE` at the bottom of the screen.
2. In the Preset dropdown in the transport bar, select **`— new empty —`**.

## 3. Timeline panel layout

```
┌──────────────────────────────────────────────────────────────────────┐
│ ▶  ⏸  ■   0.00 / 12.00  [────scrubber────]  Preset ▾  SAVE  ⊕ FX  ● REC  AUTO KEY  + TRACK  ✎ TEXT  ↑ IMPORT  ↓ EXPORT  × │  ← transport bar
├────────────────┬───────────────────────────────┬─────────────────────┤
│  TRACKS        │  dopesheet lanes               │  KEY INSPECTOR      │
│  observer.di…  │  ──●──────────●──────────      │  Path               │
│  look.exposure │  ────●────────●────────        │  Time / Ease        │
│  …             │                                │  Value              │
│                │                                │  USE TIME  CAPTURE  │
│                │                                │  SET KEY  DELETE    │
└────────────────┴───────────────────────────────┴─────────────────────┘
```

**Transport bar buttons:**
- `▶ ⏸ ■` — Play / Pause / Stop
- **Time / Duration fields** — two editable number inputs showing `current / total`. Click either field and type a value to seek or change the timeline duration.
- Preset dropdown — load a saved preset or `— new empty —` to start fresh
- `SAVE` — download the current draft using its linked filename
- `⊕ FX` — open the **Motion Functions** panel (see section 5)
- `● REC` — open the recording / screenshot modal
- `AUTO KEY` — capture changed controls as keyframes at the current time
- `+ TRACK` — add a new blank track using the path typed in the Key Inspector
- `✎ TEXT` — add or edit an annotation event at the current time
- `↑ IMPORT` — import a `.json` timeline file from disk
- `↓ EXPORT` — download the current draft as a `.json` file
- `×` — close the timeline panel

**Dopesheet columns:**
- **Track list** (left) — parameter tracks plus annotation channels; click a row to select that lane
- **Lanes** (center) — keyframe dots and annotation bars; **Ctrl+click** to add/remove from a multi-selection; **double-click** a keyframe dot to select all keyframes across all tracks that share the same time; click the ruler to seek
- **Inspector** (right) — switches between the Key Inspector and the Text Event inspector depending on what you selected

## 4. Fast workflow (recommended)

This is the easiest way to animate many controls quickly.

1. Open the Timeline panel (see section 2).
2. Scrub to your first keyframe time (usually `0.00`).
3. Click **`AUTO KEY`** in the transport bar.
   - This stores a **baseline snapshot**.
4. Change settings in the **Controls panel** (observer distance, accretion mode, jet, GRMHD, look, etc.).
5. Scrub to the next keyframe time (example `5.00`).
6. Click **`AUTO KEY`** again.
   - Changed values are written as tracks with keys at both the previous and current time.
7. Repeat steps 4–6 for more keyframes.
8. Press `▶` in the transport bar to preview.

Changes **auto-apply** immediately — no separate Apply step is needed.

## 5. Motion Functions (⊕ FX)

Click `⊕ FX` in the transport bar to open the **Motion Functions** panel. It floats above the timeline and lets you generate whole blocks of keyframes from a preset motion type.

| Type | What it animates | Tracks created |
|---|---|---|
| **Intro: sky reveal** | Starts close to the hole and widens the view into the sky / full scene | Camera transform tracks |
| **Orbit around BH** | Camera circles the black hole at constant radius and elevation | `camera.position.x/y/z`, `camera.quaternion.x/y/z/w` |
| **Zoom in / out** | Observer distance change | `observer.distance` |
| **Exposure fade** | Exposure ramp | `look.exposure` |
| **Inclination sweep** | Camera latitude sweep | `observer.orbital_inclination` |

**Shared parameters for all types:**
- **Start time** — defaults to the current playhead position
- **Duration** — total seconds the motion spans
- **Ease** — easing curve applied between steps

**Orbit-specific parameters:**
- **Number of orbits** — how many full 360° circles (fractional values allowed, e.g. `0.5` = semicircle)
- **Direction** — Counter-clockwise (CCW) or Clockwise (CW) as seen from above

The orbit always uses linear easing between keyframes (32 steps/orbit), which ensures perfectly constant angular velocity with no slowdown at intermediate points. Generated quaternions are automatically sign-normalised to prevent the camera from snapping.

Press **APPLY** to write the keyframes into the draft, or press `Escape` / click outside the panel to close it without applying.

## 6. What Auto Key does exactly

- First click: captures baseline only.
- Next clicks: compares current controls vs previous baseline snapshot.
- For each changed value, writes:
  - a key at the previous snapshot time
  - a key at the current time
- If no changes are detected, no tracks are added.

This lets you animate many parameters without manually entering each path.

## 7. Manual track/keyframe editing

Use the **Key Inspector** (right column of the timeline panel) for direct control.

Fields:
- `Path` — parameter path (example `observer.distance`, `accretion_mode`, `cameraPan.x`)
- `Time` — time in seconds
- `Ease` — `linear`, `smooth`, `smoother`
- `Value` — value at that key

Buttons in the inspector:
- `USE TIME` — fills the Time field from the current playhead position
- `CAPTURE` — reads the current runtime value for the entered Path
- `SET KEY` — add or replace the key at this time
- `DELETE KEY` — remove the selected key

Click a **track row** in the left column to inspect that track.  
Click a **keyframe dot** in the dopesheet to select that specific key.  
**Ctrl+click** additional dots to build a multi-selection across one or many tracks.

Use `✎ TEXT` in the transport bar to create an annotation event at the current time. Selecting an annotation bar switches the right column into the **Text Event** inspector, where you can edit title, body, color, width, fade-in, placement, and duration/end-marker behavior.

## 8. Interpolation rules

- Number values interpolate smoothly between keys.
- Boolean/string values are stepped (switch at the key time, no interpolation).

Practical tip:
- For mode switches like `accretion_mode`, `kerr_mode`, `jet.mode` use events (`set`) in the JSON rather than numeric tracks. Text annotations are editable in the UI; other discrete events are easiest to author via **JSON import/export** (see section 13).

## 9. Keyboard shortcuts

The timeline panel responds to keyboard shortcuts when it is focused:

| Shortcut | Action |
|---|---|
| `Space` | Play / Pause |
| `Delete` or `Backspace` | Delete all selected keyframes |
| `Ctrl+A` | Select all keyframes on the active track; with `Shift` held (`Ctrl+Shift+A`), select all keyframes across all tracks |
| `Shift+A` | Select all keyframes at the current playhead time across tracks |
| Double-click keyframe | Select all keyframes across all tracks at that same time |
| `Ctrl+C` | Copy selected keyframes to clipboard |
| `Ctrl+V` | Paste copied keyframes at the current playhead time |
| `Ctrl+Z` | Undo last edit |
| `Ctrl+Y` / `Ctrl+Shift+Z` | Redo |
| `Escape` | Close the Motion Functions panel (if open) |
| `Home` | Seek to `0.00` |
| `End` | Seek to the end of the timeline (duration) |
| `←` / `→` | Nudge playhead by ±0.1 s |
| `Shift+←` / `Shift+→` | Nudge playhead by ±1.0 s |

## 10. Undo / Redo

Every mutating operation (SET KEY, DELETE KEY, AUTO KEY, + TRACK) pushes a snapshot onto an undo stack (max 40 entries). Use **Ctrl+Z** to undo and **Ctrl+Y** (or **Ctrl+Shift+Z**) to redo.

## 11. Changes are live — no APPLY needed

Every edit (SET KEY, DELETE KEY, AUTO KEY, + TRACK) immediately updates the running timeline.

The old explicit **APPLY** button no longer exists.

## 12. Panel state persistence

When you close the timeline panel, its state is saved to session storage automatically:
- The selected preset
- The current draft (all tracks and keyframes)
- The selected track and keyframe selection

Reopening the panel restores everything exactly as you left it. State persists for the browser session; it is cleared on page reload.

## 13. JSON import / export

Use the transport bar buttons to move between the UI editor and raw JSON:

- `SAVE` — re-download the current draft with its linked filename
- `↑ IMPORT` — opens a file picker; loads a `.json` timeline file into the editor
- `↓ EXPORT` — downloads the current draft as `<name>.json`

Use exported JSON to hand-edit mode-switch events, compile flags, annotation channels, manual annotation box placement, and other advanced fields not fully exposed in the dopesheet UI. See `docs/presentation-json.md` for the full schema.

## 14. Playback and recording

Playback controls and recording settings live in the **Timeline panel** (`▲ TIMELINE`), not in the Animations panel.

Open `▲ TIMELINE` → click `● REC` to access:
- Loop toggle
- Annotation and parameter-HUD visibility toggles
- Whether overlays are included in recordings
- Reset-on-record-start toggle
- Recording quality, mode, resolution, FPS, bitrate
- `PNG SNAPSHOT`, `START REC`, and `STOP REC`

The `◀ ANIMATIONS` panel only loads built-in presentation presets and shows current status.

## 15. Save as a reusable preset

`New Preset` / `— new empty —` is a runtime editing mode.  
To make a preset that persists after a page reload:

`SAVE` only downloads the current draft with its linked filename; it does not write back into the repository automatically.

1. Click `↓ EXPORT` in the transport bar.
2. Place the downloaded file in `js/app/presentation/presets/`.
3. Add an entry in `js/app/presentation/presets/manifest.json`.
4. Reload the page.

The preset now appears in the Preset dropdown in both the Animations panel and the Timeline panel.

## 16. Troubleshooting

**Timeline panel doesn't open:**
- Click the `▲ TIMELINE` button at the bottom of the screen.
- Or select `New Preset` in the Animations panel (`◀ ANIMATIONS`).

**My edits don't appear during playback:**
- Changes auto-apply; if playback seems stale, press Stop then Play again.

**AUTO KEY reports no changes:**
- You likely clicked it twice without changing any controls between clicks.
- Change controls (in the Controls panel on the right) between AUTO KEY presses.

**A keyframe does nothing:**
- The path may be wrong or unsupported.
- Type the path in the Inspector, click `CAPTURE` to verify it resolves to a value.
- Check path spelling against `docs/presentation-json.md`.

**Camera rotation was not captured by AUTO KEY:**
- AUTO KEY captures the orbit camera transform as `camera.position.*` and `camera.quaternion.*`.
- Make sure you actually moved the camera between clicks.
