# Black Hole Simulation

A real-time, GPU-accelerated ray-tracing simulation of a black hole with an accretion disk, relativistic jets, and a full suite of general-relativistic optical effects. Runs entirely in the browser using WebGL and [three.js](http://threejs.org).

## Features

- **Schwarzschild & Kerr black holes** — exact null geodesic integration (Binet equation) with perturbative Kerr frame-dragging
- **Three accretion disk models** — thin disk (Shakura–Sunyaev), thick torus (ADAF/RIAF), and slim disk (super-Eddington)
- **Relativistic effects** — gravitational redshift, Doppler shift, relativistic beaming (physical D³ or cinematic), aberration, time dilation
- **Relativistic jets** — simple parabolic or GRMHD-calibrated physical model with spine/sheath structure, reconfinement knots, and Blandford–Znajek power scaling
- **Black-body spectrum** — temperature-dependent disk coloring with precomputed Planck lookup
- **Interactive controls** — spin, temperature, observer orbit, quality presets, and many more via dat.GUI
- **Multiple tone mapping modes** — ACES Filmic, AgX, and Scientific (logarithmic inferno colormap)

## Physics Documentation

See **[docs/physics.html](docs/physics.html)** for a comprehensive description of every physics model, equation, and approximation used in the simulation, with full academic references and comparison to real observations (EHT, *Interstellar*, Luminet 1979).

## Quick Start

Clone or download this repository, then launch a local HTTP server:

```bash
python -m http.server 8000
```

Open `http://localhost:8000` in a modern browser (Chrome or Firefox recommended). A decent GPU is required for smooth rendering.

### Performance tips

| Action | Effect |
|--------|--------|
| Lower quality preset (GUI → Quality) | Reduces integration steps and supersampling |
| Shrink the browser window | Fewer pixels to trace |
| Disable the planet | Removes intersection tests per ray |
| Disable RK4 | Falls back to faster Euler integration (less accurate near photon sphere) |

## Controls

- **Drag** to rotate the camera (orbit controls)
- **Scroll** to zoom in/out
- **GUI panel** (right side) to adjust all simulation parameters

### Key GUI parameters

| Parameter | Description |
|-----------|-------------|
| **a/M** | Dimensionless black hole spin (0 = Schwarzschild, 1 = extremal Kerr) |
| **temperature** | Accretion disk peak temperature in Kelvin (4,500 – 30,000 K) |
| **disk model** | Thin disk, thick torus (ADAF), or slim disk |
| **doppler shift** | Toggle relativistic color/brightness shifting |
| **beaming mode** | Physical (D³ Liouville) or cinematic (softened) |
| **jet model** | Off, simple, or physical (GRMHD-calibrated) |
| **observer speed** | Orbital velocity around the black hole |
| **quality** | Fast / Medium / High presets |

## Project Structure

```
index.html                Web page and shader loading
main.js                   Application logic, GUI, observer physics
raytracer.glsl            GPU ray-tracer (all physics equations)
style.css                 Styling
three-js-monkey-patch.js  Three.js compatibility patch
docs/physics.html         Comprehensive physics documentation
js-libs/                  Third-party libraries (three.js, dat.GUI, etc.)
```

## How It Works

1. Each screen pixel casts a ray from the camera into the scene.
2. The ray direction is transformed for **relativistic aberration** if the observer is moving.
3. The ray is traced by integrating the **Binet equation** (Schwarzschild geodesic ODE) using RK4 or Euler, with an optional Kerr frame-dragging correction.
4. At each integration step, intersections with the accretion disk, jets, and planet are tested.
5. **Doppler shift, gravitational redshift, and beaming** are applied to the emission.
6. The background sky (Milky Way panorama + stars) is rendered with optional Doppler color shifting.
7. The final HDR color is tone-mapped for display.

## License

See [COPYRIGHT.md](COPYRIGHT.md) for license and copyright information.

Originally based on [oseiskar/black-hole](https://github.com/oseiskar/black-hole).

**AI Disclaimer**: AI helped a lot in the code, like for translating the complex general relativity equations into functional WebGL shaders, as well as structuring the academic references. Fortunately, the original code was well-written and was a really good, human-made base. Physics research papers are all open access and properly cited in the documentation, so you can verify every equation and model against the primary sources.