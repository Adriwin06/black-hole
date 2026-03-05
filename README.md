# Black Hole Simulation

A real-time, GPU-accelerated ray-tracing simulation of a black hole with an accretion disk, relativistic jets, and a full suite of general-relativistic optical effects. Runs entirely in the browser using WebGL and [three.js](http://threejs.org).

**[Live Demo](https://adriwin06.github.io/black-hole)** — Chrome or Firefox on a dedicated GPU recommended.

> This is a substantially extended fork of [oseiskar/black-hole](https://github.com/oseiskar/black-hole). See [What's new](#whats-new-in-this-fork) for a summary of additions.

---

## Features

### Physics & Rendering
- **Schwarzschild & Kerr black holes** — two geodesic integration modes:
  - *Schwarzschild* Binet ODE (exact for a = 0, fast)
  - *True Kerr* Carter (1968) separated equations in Mino time — correct D-shaped shadow for spinning holes
- **Three accretion disk models** — thin disk (Shakura–Sunyaev), thick torus (ADAF/RIAF), and slim disk (super-Eddington)
- **GRMHD-calibrated accretion** — magnetization (σ), electron temperature ratio (R_high), MAD/SANE magnetic flux, MRI turbulence, and kappa-distribution electrons
- **Relativistic effects** — gravitational redshift, Doppler shift, relativistic beaming (physical D³ Liouville or cinematic), aberration, time dilation
- **Relativistic jets** — simple parabolic or physical GRMHD-calibrated model with spine/sheath structure, reconfinement shocks, jet-corona connection, and Blandford–Znajek power scaling
- **Black-body spectrum** — temperature-dependent disk coloring with precomputed Planck lookup
- **Multiple tone-mapping modes** — ACES Filmic, AgX, and Scientific (logarithmic inferno colormap)
- **Multi-pass bloom** — threshold → mip-chain Gaussian blur → weighted composite

### Post-Processing & Quality
- **Temporal Anti-Aliasing (TAA)** — history accumulation with motion rejection and clip-box clamping for artifact-free still frames
- **Five quality presets** — from Mobile (28 steps, 0.55× res) to Ultra (1400 Kerr steps, 6× supersampling), plus Cinematic mode for offline rendering
- **Auto GPU benchmark** — measures frame time on first load and selects the Optimal preset automatically

### User Interface
- **dat.GUI control panel** — resizable right-side panel with collapsible folders for every parameter
- **Astrophysical presets** — one-click configurations for M87\*, Sgr A\*, Cygnus X-1, GRS 1915+105, and more, each sourced from published observations
- **Observer controls** — free orbital motion with adjustable inclination, azimuth, distance, and optional automatic circular orbit

### Presentation & Recording
- **Presentation Timeline** — bottom-docked dopesheet editor (inspired by Blender / After Effects) for scripted keyframe animations; supports linear, smooth, and smoother easing
- **Built-in animation presets** — *Full Feature Tour* (186 s), *Orbit Showcase*, *Horizon Dive*, *Hover Blueshift*
- **WebM video recording** — capture the canvas directly to a `.webm` file via MediaRecorder, rate-controlled by the timeline playback
- **High-quality offline rendering** — Cinematic preset with manually boosted supersampling for publication-quality stills and video frames 

---

## Physics Documentation

See **[docs/physics.html](docs/physics.html)** for a comprehensive description of every physics model, equation, and approximation used in the simulation, with full academic references and comparison to real observations (EHT, *Interstellar*, Luminet 1979).

---

## Quick Start

Clone or download this repository, then launch a local HTTP server:

```bash
python -m http.server 8000
```

Open `http://localhost:8000` in a modern browser (Chrome or Firefox recommended). A dedicated GPU is required for smooth rendering.

### Performance tips

| Action | Effect |
|--------|--------|
| Lower quality preset (GUI → Quality) | Reduces integration steps and supersampling |
| Shrink the browser window | Fewer pixels to trace |
| Disable the planet | Removes ray-sphere intersection tests |
| Disable RK4 | Falls back to Euler integration (faster, less accurate near photon sphere) |
| Switch to Schwarzschild mode | Binet equation is ~4× cheaper than full Kerr Carter equations |

---

## Controls

- **Left drag** — orbit the camera
- **Right drag** — pan
- **Left + Right drag** — roll
- **Scroll** — zoom in/out
- **Controls panel** (right side, ☰) — all simulation parameters
- **Timeline panel** (bottom, ▲ TIMELINE) — animation playback and recording

### Key GUI parameters

| Parameter | Description |
|-----------|-------------|
| **a/M** | Dimensionless black hole spin (0 = Schwarzschild, 1 = extremal Kerr) |
| **Kerr mode** | Switch between Schwarzschild Binet and full Carter Kerr geodesics |
| **temperature** | Accretion disk peak temperature in Kelvin (4,500 – 30,000 K) |
| **disk model** | Thin disk, thick torus (ADAF), or slim disk |
| **doppler shift** | Toggle relativistic color/brightness shifting |
| **beaming mode** | Physical (D³ Liouville) or cinematic (softened) |
| **jet model** | Off, simple, or physical (GRMHD-calibrated) |
| **observer speed** | Orbital velocity around the black hole |
| **quality** | Mobile / Optimal / Medium / High / Ultra / Cinematic |

### Quality preset levels

| Preset | Steps (std / Kerr) | Supersampling | Description |
|--------|-------------------|---------------|-------------|
| Mobile | 28 / 120 | 1× | 0.55× resolution + TAA; fastest |
| Optimal | 100 / 400 | 1× | 0.8× res + TAA; recommended default |
| Medium | 100 / 400 | 1× / 3× | Full resolution; balanced |
| High | 320 / 520 | 4× | Full resolution; GPU-intensive |
| Ultra | 600 / 1400 | 4× / 6× | Maximum fidelity |
| Cinematic | 600 / 1400 | 6x / 12x | Offline rendering quality |

### Astrophysical presets

| Preset | Object | Notes |
|--------|--------|-------|
| Default | Generic BH | a/M = 0.90, thin disk |
| M87\* | Virgo A SMBH | a/M ≈ 0.90, thick torus, physical jet; first EHT image (2019) |
| Sgr A\* | Milky Way centre | a/M ≈ 0.50, ADAF torus; EHT image (2022) |
| Cygnus X-1 | X-ray binary | a/M ≈ 0.99, thin disk (continuum-fitting spin) |
| GRS 1915+105 | Microquasar | Near-extremal spin, slim disk |

---

## How It Works

1. Each screen pixel casts a ray from the camera into the scene.
2. The ray direction is transformed for **relativistic aberration** if the observer is moving.
3. The ray is integrated using either the **Schwarzschild Binet equation** (Euler / RK4) or the **Kerr Carter equations** in Mino time, depending on the selected mode.
4. At each step, intersections with the accretion disk, GRMHD medium, jets, and planet are tested and composited using Beer–Lambert transmittance.
5. **Doppler shift, gravitational redshift, and beaming** are applied to each emission source.
6. The background sky (Milky Way panorama + star field) is rendered with optional Doppler color shifting.
7. The HDR accumulation buffer is bloom-composited and then tone-mapped to sRGB.
8. Optionally, TAA accumulates multiple jittered frames before display.

---

## Project Structure

The codebase is organized into logical modules by function and responsibility.

### **GLSL Shaders** (`shaders/raytracer/`)

```
shaders/raytracer/
├── core/                          # Foundational definitions
│   ├── defines.glsl              # Constants, macros, uniforms, rendering params
│   └── math.glsl                 # Math utilities, coordinate transforms, FBM noise
├── physics/                       # Physics models
│   ├── geodesics.glsl            # Schwarzschild Binet + Kerr Carter (Mino time)
│   ├── accretion.glsl            # Thin disk, ADAF torus, slim disk, GRMHD turbulence
│   ├── jet.glsl                  # Jet models: simple parabolic + physical GRMHD
│   ├── planet.glsl               # Planet ray-sphere intersection
│   └── background.glsl           # Galaxy/star background rendering
└── output/                        # Rendering pipeline
    ├── tonemapping.glsl          # ACES Filmic, AGX, scientific tone-mappers
    ├── trace_ray.glsl            # Core ray-marching loop with Beer-Lambert composite
    └── main.glsl                 # GLSL main() entry point
```

### **JavaScript Modules** (`js/app/`)

```
js/app/
├── bootstrap.js                    # Entry point: fetches GLSL shards & textures, calls init()
├── core/                           # System core
│   ├── observer.js                 # Observer entity, SR orbital mechanics, time dilation
│   ├── shader.js                   # Shader class, compile-time Mustache parameters
│   └── renderer.js                 # Three.js scene, init, TAA, bloom, render loop
├── scene/                          # Scene management
│   └── camera.js                   # Camera initialization & per-frame updates
├── graphics/                       # Graphics effects
│   └── bloom.js                    # Multi-pass mip-chain Gaussian bloom
├── presentation/                   # Presentation & recording system
│   ├── presentation-controller.js  # Keyframe timeline engine, MediaRecorder hooks
│   ├── presentation-gui.js         # ANIMATIONS panel UI (preset selector, status)
│   ├── timeline-panel.js           # Bottom dopesheet panel (transport, key inspector)
│   └── presets/                    # Built-in animation sequences (JSON)
│       ├── manifest.json
│       ├── full-feature-tour.json
│       ├── orbit-showcase.json
│       ├── horizon-dive.json
│       └── hover-blueshift.json
└── ui/                             # User interface
    ├── presets.js                  # Astrophysical black hole preset library
    ├── quality-presets.js          # Rendering quality preset library
    └── gui.js                      # dat.GUI panel setup and parameter wiring
```

### **Other files**

```
index.html                          # Web page entry point
style.css                           # Styling (panels, timeline, controls)
three-js-monkey-patch.js            # Three.js r71 compatibility patches
js-libs/                            # Third-party libraries (three.js, dat.GUI, webm-muxer, …)
docs/
└── physics.html                    # Comprehensive physics documentation
└── presentation-editor.md          # Guide to using the presentation timeline editor
└── presentation-json.md
```

---

## What's New in This Fork

Additions over the [upstream oseiskar/black-hole](https://github.com/oseiskar/black-hole):

| Feature | Details |
|---------|---------|
| True Kerr geodesics | Carter (1968) separated ODEs in Mino time — correct D-shaped shadow |
| GRMHD accretion model | σ, R_high, MAD flux, MRI turbulence, κ-distribution electrons |
| Presentation Timeline | Keyframe dopesheet editor with transport controls and easing curves |
| WebM recording | MediaRecorder capture synced to timeline playback |
| Built-in animation presets | Full Feature Tour, Orbit Showcase, Horizon Dive, Hover Blueshift |
| Astrophysical BH presets | M87\*, Sgr A\*, Cygnus X-1, GRS 1915+105 from published data |
| Temporal Anti-Aliasing | Motion-rejection TAA for artifact-free high-quality frames |
| Five quality tiers | Mobile through Ultra with per-mode step counts and supersampling |
| Three tone-mappers | ACES Filmic, AGX, Scientific (inferno colormap) |
| Three accretion models | Thin disk, thick torus (ADAF), slim disk (super-Eddington) |
| Physical GRMHD jet | Spine/sheath, reconfinement shocks, Blandford–Znajek power scaling |
| Resizable UI panels | Drag-to-resize controls panel and timeline |

---

## License

See [COPYRIGHT.md](COPYRIGHT.md) for full license and copyright information.

Originally based on [oseiskar/black-hole](https://github.com/oseiskar/black-hole) (MIT).  
Fork maintained and substantially extended by [Adriwin](https://github.com/Adriwin06).

**AI Disclaimer**: AI assisted with translating complex general-relativity equations into functional WebGL shaders and with structuring academic references. The original codebase was human-made and well-written. All physics references are open-access and cited in [docs/physics.html](docs/physics.html) so every equation and model can be verified against primary sources.
