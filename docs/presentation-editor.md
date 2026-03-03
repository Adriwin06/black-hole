# Presentation Timeline Editor Guide

This guide explains how to use the in-app timeline editor in:

`ANIMATIONS > PRESENTATION TIMELINE`

If you want the raw JSON schema, see `docs/presentation-json.md`.

## 1. Important: how to open the editor

The editor is only shown when **Preset = `New Preset`**.

Steps:
1. Open `ANIMATIONS > PRESENTATION TIMELINE`.
2. In the `Preset` dropdown, select **`New Preset`**.
3. Expand the `TIMELINE EDITOR` header if it is collapsed.

If you select a normal preset (like `Full Feature Tour`), the editor is hidden.

## 2. Fast workflow (recommended)

This is the easiest way to animate many controls quickly.

1. Choose `Preset = New Preset`.
2. Set `Name`, `Duration`, and optional `Loop`.
3. Move timeline `Time` to your first keyframe time (usually `0.00`).
4. Click **`AUTO KEYFRAME (FROM CONTROLS)`**.
   - This stores a **baseline snapshot**.
5. Change settings from the normal controls panel (observer distance, accretion mode, jet, GRMHD, look settings, etc.).
6. Move timeline time to the next keyframe time (example `5.00`).
7. Click **`AUTO KEYFRAME (FROM CONTROLS)`** again.
   - Changed values are added to tracks.
   - Keys are created at both previous and current times for those changed settings.
8. Repeat steps 5-7 for more keyframes.
9. Click **`APPLY`**.
10. Press `PLAY` to preview.

## 3. What Auto Keyframe does exactly

- First click: captures baseline only.
- Next clicks: compares current controls vs previous snapshot.
- For each changed value, it writes:
  - key at previous snapshot time
  - key at current time
- If no changes are detected, no tracks are added.

This lets you animate many parameters without manually entering each path.

## 4. Manual track/keyframe editing

Use **Tracks / keyframes** for direct control.

Fields:
- `Path`: parameter path (example `observer.distance`, `accretion_mode`, `cameraPan.x`)
- `Key`: time in seconds
- `Ease`: `linear`, `smooth`, `smoother`
- `Value`: value at that time

Buttons:
- `USE TIME`: fills key time from current timeline time
- `CAPTURE VALUE`: reads current runtime value for `Path`
- `ADD/UPDATE`: add key or replace key at same time
- `REMOVE`: remove nearest key for this path/time

Track list actions:
- `GO`: seek timeline to that key
- `EDIT`: load that key into editor fields
- `DEL`: delete key (or whole track with `DEL TRACK`)

## 5. Manual event editing

Use **Events** for discrete actions (toggles/mode changes/annotation popups).

Supported actions:
- `set`
- `updateShader`
- `startDive`, `pauseDive`, `resetDive`
- `startHover`, `pauseHover`, `resetHover`
- `annotation`
- `clearAnnotation`

For `set`, provide:
- `Path`
- `Value`
- optional `Force compile`

For `annotation`, provide note fields like:
- title
- text
- anchor target
- placement

Event list actions:
- `GO`, `EDIT`, `DEL`

## 6. Interpolation rules

- Number values interpolate between keys.
- Boolean/string values are stepped (switch behavior).

Practical tip:
- For mode switches like `accretion_mode`, `kerr_mode`, `jet.mode`, prefer **events** (`set`) instead of numeric tracks.

## 7. APPLY is required

The editor builds a draft timeline.  
Changes are not active until you click **`APPLY`**.

After `APPLY`:
- timeline is loaded into runtime
- `PLAY`, scrub, and recording use the updated timeline

## 8. JSON tools

Inside the editor:
- `EXPORT`: write draft JSON to the JSON text box
- `IMPORT FILE`: open file picker and import a `.json` timeline file into draft
- `IMPORT TEXT`: parse JSON text from the JSON text box into draft
- `DOWNLOAD`: save draft as a `.json` file

Use this to move between UI editing and raw JSON editing.

## 9. Save as a reusable preset

`New Preset` is a runtime/editor mode.  
To make a reusable preset:

1. Click `DOWNLOAD` in editor.
2. Put the file in `js/app/presentation/presets/`.
3. Add an entry in `js/app/presentation/presets/manifest.json`.
4. Reload page.

Now it appears in normal preset dropdown.

## 10. Troubleshooting

Editor not visible:
- Select `Preset = New Preset`.

I changed things but playback did not change:
- Click `APPLY`.

Auto keyframe says no changes:
- You likely pressed it twice without changing controls.
- Change controls between clicks.

A keyframe does nothing:
- Path may be wrong or unsupported.
- Check path spelling and try `CAPTURE VALUE`.

Camera rotation was not detected by auto keyframe:
- Newer builds capture orbit camera transform (`camera.position.*` + `camera.quaternion.*`) in AUTO KEYFRAME.
- If you still see this, make sure you actually moved the camera between clicks and then press `APPLY`.
