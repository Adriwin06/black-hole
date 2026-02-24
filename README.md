---
---

# Ray-traced simulation of a black hole

In this simulation, the light ray paths are computed by integrating an ODE describing the Schwarzschild geodesics using GLSL on the GPU, leveraging WebGL and [three.js](http://threejs.org). This should result to a fairly physically accurate gravitational lensing effect. Various other relativistic effects have also been added and their contributions can be toggled from the GUI.
The simulation has normalized units such that the Schwarzschild radius of the black hole is one and the speed of light is one length unit per second (unless changed using the "time scale" parameter).

See **[this page](https://oseiskar.github.io/black-hole/docs/physics.html)** ([PDF version](https://oseiskar.github.io/black-hole/docs/physics.pdf)) for a more detailed description of the physics of the simulation.

### Interstellar-style controls

The GUI now includes controls matching the visual effects discussed in the Interstellar black-hole breakdown:

* **Accretion disk > temperature (K):** directly adjusts disk black-body color over the 4,500 K to 30,000 K range (orange to blue-white).
 * **Black hole > rotating shadow:** enables/disables spin-induced asymmetry of the black-hole shadow.
 * **Black hole > a/M:** sets the dimensionless spin parameter.
 * **Black hole > shadow squeeze:** scales the spin deformation strength (use lower values for a film-like reduced effect).
 * **Relativistic effects > doppler shift / beaming:** controls blue-shifted brightening of approaching plasma and red-shifted dimming of receding plasma.

### System requirements

The simulation needs a decent GPU and a recent variant of Chrome or Firefox to run smoothly. In addition to changing simulation quality from the GUI, frame rate can be increased by shrinking the browser window and/or reducing screen resolution. Disabling the planet from the GUI also increases frame rate.

Example: runs 30+ fps at resolution 1920 x 1080 in Chrome 48 on a Linux desktop with GeForce GTX 750 Ti and "high" simulation quality

### Known artefacts

 * The accretion disk is stylized and procedurally textured for visual turbulence; it is not a full GRMHD/volumetric simulation.
 * The spectrum used in modeling the Doppler shift of the Milky Way background image is quite arbitrary (not based on real spectral data) and consequently the Doppler-shifted background colors may be wrong.
 * The lighting model of the planet is based on a point-like light source and a quite unphysical ambient component.
 * In the "medium" quality mode, the planet deforms unphysically when it travels between the camera and the black hole.
 * The light paths bend a bit more than they should due to low ODE solver step counts (see [numeric tests](https://github.com/oseiskar/black-hole/blob/numeric-notebooks/numeric_tests.ipynb)), but this seems to happen in a systematic way so that the image looks very similar in comparison to a more accurate simulation.
 * Lorentz contraction causes jagged looks in the planet when simultaneously enabled with "light travel time" and the planet is close to the black hole.
 * Texture sampling issues cause unintended star blinking.

_see **[COPYRIGHT.md](https://github.com/oseiskar/black-hole/blob/master/COPYRIGHT.md)** for license and copyright info_
