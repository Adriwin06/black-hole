---
---

# Black Hole Simulation

![Banner](assets/gif/banner.gif)

A real-time, GPU-accelerated browser visualization of a black hole with an accretion flow, jet models, and relativistic optical effects. Runs entirely in the browser using WebGL and [three.js](https://threejs.org).

**[Live Demo](https://adriwin06.github.io/black-hole)** â€” Chrome or Firefox on a dedicated GPU recommended.

> This is a substantially extended fork of [oseiskar/black-hole](https://github.com/oseiskar/black-hole). See [What's new](#whats-new-in-this-fork) for a summary of additions.

> Scientific scope: Schwarzschild photon geodesics follow the exact Binet equation, but rendered fidelity still depends on the integrator, step budget, and adaptive stepping. Spin, accretion, jet, and GRMHD-related options use a mix of analytic, semi-analytic, and heuristic approximations; returning radiation is not ray-traced, and jet colours use an effective-temperature proxy rather than a frequency-resolved synchrotron spectrum. See [docs/physics.html](docs/physics.html) for what is exact and what is approximate.

---

## Features

### Physics & Rendering
- **Two public spin modes** â€” `Fast (Binet lensing)` traces photons with the Schwarzschild Binet solver (exact for a = 0) plus perturbative spin heuristics; `Kerr-inspired disk velocities` keeps the same approximate photon solver but uses Kerr equatorial angular velocity to drive disk matter
- **Three accretion disk models** â€” thin disk (Shakuraâ€“Sunyaev), thick torus (ADAF/RIAF), and slim disk (super-Eddington)
- **GRMHD-inspired accretion controls** â€” plasma beta (Î²), magnetic-field strength, `R_high` electron-heating / Ti:Te proxy controls, MAD/SANE magnetic flux, MRI-inspired turbulence, and kappa-distribution electron parameters
- **Relativistic effects** â€” gravitational redshift, Doppler shift, thermal/background beaming controls (physical DÂ³ Liouville or cinematic), aberration, time dilation
- **Relativistic jets** â€” simple analytic jet or a more detailed GRMHD-inspired jet model with spine/sheath structure, reconfinement shocks, jet-corona connection, and Blandfordâ€“Znajek-inspired power scaling
- **Black-body spectrum** â€” temperature-dependent disk coloring with precomputed Planck lookup
- **Multiple tone-mapping modes** â€” ACES Filmic, AgX, and Scientific (logarithmic inferno colormap)
- **Multi-pass bloom** â€” threshold â†’ mip-chain Gaussian blur â†’ weighted composite

### Post-Processing & Quality
- **Temporal Anti-Aliasing (TAA)** â€” history accumulation with motion rejection and clip-box clamping for cleaner still frames
- **Six quality presets** â€” Mobile, Optimal, Medium, High, Ultra, and Cinematic
- **Auto GPU benchmark** â€” measures frame time on first load and keeps Optimal on capable systems or falls back to Mobile

### User Interface
- **dat.GUI control panel** â€” resizable right-side panel with collapsible folders for every parameter
- **Astrophysical presets** â€” literature-inspired starting points for M87\*, Sgr A\*, Cygnus X-1, GRS 1915+105, and more
- **Observer controls** â€” mouse orbit/pan/roll, a bottom-left observer widget with distance dial + motion toggle, and optional automatic circular orbit in the stable Schwarzschild regime (`r >= 3 r_s`)

### Presentation & Recording
- **Presentation Timeline** â€” bottom-docked dopesheet editor (inspired by Blender / After Effects) for scripted keyframe animations, preset loading, playback, and recording; supports linear, smooth, and smoother easing
- **Built-in interactive scenarios** â€” *Freefall Dive* and *Hover Approach*
- **Live scenario capture to timeline** â€” record a manual *Freefall Dive* or *Hover Approach* run, including orbit/pan camera motion and captured shader time, straight into timeline tracks, with optional camera smoothing
- **Built-in timeline presets** â€” *Full Feature Tour* (186 s) and *Orbit Showcase*
- **WebM video recording** â€” realtime MediaRecorder capture plus offline WebCodecs/WebM muxing when supported
- **Offline PNG snapshot** â€” one-click still export that forces the Cinematic (offline) quality preset before download
- **High-quality offline rendering** â€” Cinematic preset with manually boosted supersampling for publication-quality stills and video frames 

---

## Physics Documentation

See **[docs/physics.html](docs/physics.html)** for a detailed description of the models, equations, approximations, and implementation scope used in the simulation, with academic references and notes on where the renderer departs from full GRRT/GRMHD treatments.

---

## Quick Start

Clone or download this repository, then launch a local HTTP server:

```bash
python -m http.server 8000
```

Open `http://localhost:8000` in a modern browser (Chrome or Firefox recommended). A dedicated GPU is recommended for smooth rendering.

### Performance tips

| Action | Effect |
|--------|--------|
| Lower quality preset (GUI â†’ Quality) | Reduces integration steps and supersampling |
| Shrink the browser window | Fewer pixels to trace |
| Disable the planet | Removes ray-sphere intersection tests |
| Disable RK4 | Falls back to leapfrog / StÃ¶rmer-Verlet integration (faster, less accurate near the photon sphere) |
| Switch solver mode to Fast | Uses the lightweight Binet photon solver everywhere |

---

## Controls

- **Left drag** â€” orbit the camera
- **Right drag** â€” pan
- **Left + Right drag** â€” roll
- **Scroll** â€” zoom in/out
- **Controls panel** (right side, `CONTROLS â–¶`) â€” all simulation parameters
- **Animations panel** (left side, `â—€ ANIMATIONS`) â€” Freefall Dive, Hover Approach, and live capture of those runs into the timeline
- **Observer widget** (bottom-left XYZ indicator) â€” distance dial, motion toggle, camera reset
- **Timeline panel** (bottom, `â–² TIMELINE`) â€” preset loading, timeline editing, playback, recording, and export/import

### Key GUI parameters

| Parameter | Description |
|-----------|-------------|
| **a/M** | Signed black hole spin; positive and negative values mirror the spin direction |
| **solver mode** | `Fast (Binet lensing)` or `Kerr-inspired disk velocities`; the latter keeps approximate photon lensing but uses Kerr angular velocity to drive disk matter |
| **temperature (K)** | Visualized disk color temperature in Kelvin (4,500 â€“ 30,000 K) |
| **disk model** | Thin disk, thick torus (ADAF), or slim disk |
| **doppler shift (color)** | Toggle relativistic red/blue spectral shifting for black-body-based thermal emitters and the background-sky proxy; jets keep their own synchrotron-motivated transfer |
| **physical (DÂ³ Liouville)** | Use physically motivated beaming for thermal emitters and the background-sky proxy instead of the softened cinematic curve |
| **jet enabled / mode** | Toggle jets and choose simple or more detailed GRMHD-inspired shading |
| **observer motion** | Toggle automatic circular orbit around the black hole; motion is clamped to the stable Schwarzschild regime (`r >= 3 r_s`) |
| **quality preset** | Mobile / Optimal / Medium / High / Ultra / Cinematic |

### Quality preset levels

| Preset | Steps (fast / Kerr-inspired) | Supersampling | Description |
|--------|-------------------|---------------|-------------|
| Mobile | 28 / 120 | 1Ã— | 0.55Ã— resolution + TAA; fastest |
| Optimal | 100 / 400 | 1Ã— | 0.8Ã— res + TAA; recommended default |
| Medium | 100 / 400 | 1Ã— / 3Ã— | Full resolution; balanced |
| High | 320 / 520 | 4Ã— | Full resolution; GPU-intensive |
| Ultra | 600 / 1400 | 4Ã— / 6Ã— | Maximum fidelity |
| Cinematic | 600 / 1400 | 6x / 12x | Offline rendering quality |

### Built-in black hole presets

| Preset | Object | Notes |
|--------|--------|-------|
| Default | Generic BH | a/M = 0.90, thin disk |
| M87\* | Virgo A SMBH | Illustrative high-spin, MAD-like thick-torus + jet configuration inspired by EHT-era GRMHD studies |
| Sgr A\* | Milky Way centre | Illustrative moderate-spin ADAF/RIAF-style torus; EHT modeling does not determine a unique spin |
| Cygnus X-1 | X-ray binary | Near-extremal thin-disk preset inspired by continuum-fitting studies |
| GRS 1915+105 | Microquasar | Near-extremal slim-disk preset inspired by continuum-fitting and jet observations |
| Gargantua (Interstellar visuals) | Film-inspired | Warm thin disk, boosted glow, softened relativistic effects |
| Schwarzschild | Textbook case | Spin disabled, clean circular shadow |

---

## How It Works

1. Each screen pixel casts a ray from the camera into the scene.
2. The ray direction is transformed for **relativistic aberration** if the observer is moving.
3. Photon paths are traced with the **Schwarzschild Binet equation** (leapfrog / RK4). The optional `Kerr-inspired disk velocities` mode keeps the same photon solver but upgrades disk matter to a Kerr-inspired orbital-velocity model.
4. At each step, intersections with the accretion disk, GRMHD-inspired media, jets, and planet are tested and composited using Beerâ€“Lambert transmittance.
5. **Doppler shift, gravitational redshift, and beaming** are applied with model-appropriate transfer rules; jets use a separate synchrotron-motivated treatment.
6. The background sky (Milky Way panorama + star field) is rendered with optional Doppler color shifting.
7. The HDR accumulation buffer is bloom-composited and then tone-mapped to sRGB.
8. Optionally, TAA accumulates multiple jittered frames before display.

---

## Project Structure

The codebase is organized into logical modules by function and responsibility.

### **GLSL Shaders** (`shaders/raytracer/`)

```text
shaders/raytracer/
|-- core/                          # Foundational definitions
|   |-- defines.glsl               # Constants, macros, uniforms, rendering params
|   `-- math.glsl                  # Math utilities, coordinate transforms, FBM noise
|-- physics/                       # Physics models
|   |-- geodesics.glsl             # Schwarzschild Binet solver + experimental Kerr Mino-time helpers
|   |-- accretion/                 # Accretion model shards
|   |   |-- core.glsl              # Shared disk helpers and thin-disk temperature/flux logic
|   |   |-- thick-torus.glsl       # ADAF / RIAF torus model
|   |   |-- slim-disk.glsl         # Slim disk model
|   |   |-- grmhd.glsl             # GRMHD-inspired turbulence and plasma modifiers
|   |   `-- planet-irradiation.glsl # Planet temperature from disk illumination
|   |-- jet/                       # Jet model shards
|   |   |-- simple.glsl            # Simple parabolic jet
|   |   |-- physical-core.glsl     # Physical jet core model
|   |   `-- physical-grmhd.glsl    # GRMHD jet extensions
|   |-- planet.glsl                # Planet ray-sphere intersection
|   `-- background.glsl            # Galaxy/star background rendering
`-- output/                        # Rendering pipeline
    |-- tonemapping.glsl           # ACES Filmic, AgX, scientific tone-mappers
    |-- trace/                     # Core ray-marching loop shards
    |   |-- setup.glsl             # Trace setup and geodesic integration
    |   |-- accretion-thin.glsl    # Thin-disk compositing
    |   |-- accretion-thick.glsl   # Thick-torus compositing
    |   |-- accretion-slim.glsl    # Slim-disk compositing
    |   |-- jet.glsl               # Jet compositing
    |   `-- finish.glsl            # Trace finish and background compositing
    `-- main.glsl                  # GLSL main() entry point
```

### **JavaScript Modules** (`js/app/`)

```
js/app/
â”œâ”€â”€ bootstrap.js                    # Entry point: fetches GLSL shards & textures, calls init()
â”œâ”€â”€ core/                           # System core
â”‚   â”œâ”€â”€ observer.js                 # Observer state, circular-orbit kinematics, simulation time scaling
â”‚   â”œâ”€â”€ shader.js                   # Shader class, compile-time Mustache parameters
â”‚   â””â”€â”€ renderer.js                 # Three.js scene, live dive/hover capture, TAA, bloom, render loop
â”œâ”€â”€ scene/                          # Scene management
â”‚   â””â”€â”€ camera.js                   # Camera initialization & per-frame updates
â”œâ”€â”€ graphics/                       # Graphics effects
â”‚   â””â”€â”€ bloom.js                    # Multi-pass mip-chain Gaussian bloom
â”œâ”€â”€ presentation/                   # Presentation & recording system
â”‚   â”œâ”€â”€ presentation-controller.js  # Keyframe timeline engine, annotations, recording pipeline
â”‚   â”œâ”€â”€ presentation-gui.js         # Legacy presentation mini-panel helper (kept for reference)
â”‚   â”œâ”€â”€ timeline-panel.js           # Bottom dopesheet panel (transport, presets, key inspector, recording)
â”‚   â””â”€â”€ presets/                    # Built-in animation sequences (JSON)
â”‚       â”œâ”€â”€ manifest.json
â”‚       â”œâ”€â”€ full-feature-tour.json
â”‚       â””â”€â”€ orbit-showcase.json
â””â”€â”€ ui/                             # User interface
    â”œâ”€â”€ presets.js                  # Astrophysical black hole preset library
    â”œâ”€â”€ quality-presets.js          # Rendering quality preset library
    â””â”€â”€ gui.js                      # dat.GUI panel setup and parameter wiring
```

### **Other files**

```
index.html                          # Web page entry point
style.css                           # Styling (panels, timeline, controls)
three-js-monkey-patch.js            # Legacy Three.js compatibility patches
js-libs/                            # Third-party libraries (three.js, dat.GUI, webm-muxer, â€¦)
docs/
â”œâ”€â”€ physics.html                    # Comprehensive physics documentation
â”œâ”€â”€ presentation-editor.md          # Guide to using the presentation timeline editor
â””â”€â”€ presentation-json.md            # Timeline JSON schema & advanced event guide
```

---

## What's New in This Fork

Additions over the [upstream oseiskar/black-hole](https://github.com/oseiskar/black-hole):

| Feature | Details |
|---------|---------|
| WIP Kerr geodesics | Carter (1968) Mino-time integrator exists in GLSL but is not yet exposed in the UI |
| GRMHD-inspired accretion controls | plasma Î², B-field strength, R_high, MAD flux, MRI-inspired turbulence, Îº-distribution electrons |
| Presentation Timeline | Keyframe dopesheet editor with presets, transport controls, recording, and easing curves |
| Scenario capture to timeline | Records live Freefall Dive / Hover Approach runs plus camera motion into keyframed timeline data |
| WebM recording | Realtime MediaRecorder capture plus offline WebCodecs/WebMMuxer export |
| Offline PNG snapshot | One-click still image export using the Cinematic offline preset |
| Interactive observer scenarios | Freefall Dive and Hover Approach |
| Built-in timeline presets | Full Feature Tour and Orbit Showcase |
| Astrophysical BH presets | M87\*, Sgr A\*, Cygnus X-1, GRS 1915+105, Gargantua, Schwarzschild |
| Temporal Anti-Aliasing | Motion-rejection TAA for cleaner high-quality frames |
| Six quality tiers | Mobile, Optimal, Medium, High, Ultra, Cinematic |
| Three tone-mappers | ACES Filmic, AgX, Scientific (inferno colormap) |
| Three accretion models | Thin disk, thick torus (ADAF), slim disk (super-Eddington) |
| GRMHD-inspired jet model | Spine/sheath, reconfinement shocks, Blandfordâ€“Znajek-inspired power scaling |
| Resizable UI panels | Drag-to-resize controls panel and timeline |

---

## License

The source code for this fork is MIT-licensed, but the repository as distributed is a mixed-license bundle. Some bundled third-party libraries and assets use separate terms, and the shipped Milky Way panorama is not covered by the MIT code license. See [COPYRIGHT.md](COPYRIGHT.md) for the full breakdown. If you need a clean permissive redistribution, replace the restricted third-party assets first.

Originally based on [oseiskar/black-hole](https://github.com/oseiskar/black-hole) (MIT).  
Fork maintained and substantially extended by [Adriwin](https://github.com/Adriwin06).

**AI Disclaimer**: AI assisted with translating complex general-relativity equations into functional WebGL shaders and with structuring academic references. The original codebase was human-made and well-written. The main physics formulas, approximations, and external references are documented in [docs/physics.html](docs/physics.html).

---

## Demo

[![Demo video](https://img.youtube.com/vi/74kY6jmZDs8/0.jpg)](https://www.youtube.com/watch?v=74kY6jmZDs8)

