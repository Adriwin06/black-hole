---
---

# Ray-traced simulation of a black hole

In this simulation, the light ray paths are computed by integrating an ODE describing the Schwarzschild geodesics using GLSL on the GPU, leveraging WebGL and [three.js](http://threejs.org). This should result to a fairly physically accurate gravitational lensing effect. Various other relativistic effects have also been added and their contributions can be toggled from the GUI.
The simulation has normalized units such that the Schwarzschild radius of the black hole is one and the speed of light is one length unit per second (unless changed using the "time scale" parameter).

See **[this page](https://oseiskar.github.io/black-hole/docs/physics.html)** ([PDF version](https://oseiskar.github.io/black-hole/docs/physics.pdf)) for a more detailed description of the physics of the simulation.

### Physics accuracy

The simulation implements several physically accurate relativistic effects:

* **Schwarzschild geodesics:** The core light-bending equation $\ddot{u} = -u(1 - \frac{3}{2}u^2)$ is exact for null geodesics in Schwarzschild spacetime.
* **Gravitational redshift:** Correctly implemented as $z = \sqrt{(1 - r_s/r_{emit})/(1 - r_s/r_{obs})}$.
* **Relativistic Doppler effect:** Uses the full formula $f_r = f_s / [\gamma(1 + \vec{v} \cdot \hat{n})]$.
* **Relativistic beaming:** Two modes available:
  - *Physical (D³ Liouville):* Implements proper $D^3$ intensity scaling per Liouville's theorem ($I/\nu^3$ is Lorentz invariant). Creates dramatic left-right brightness asymmetry matching theoretical predictions.
  - *Cinematic:* Softened beaming (~$D^{2.15}$ with clamped Doppler factor) for artistic rendering similar to movie visualizations.
* **Accretion disk temperature profile:** Follows the Shakura-Sunyaev thin disk model $T \propto r^{-3/4}(1 - \sqrt{r_{in}/r})^{1/4}$.
* **ISCO (Innermost Stable Circular Orbit):** Now dynamically calculated using the Bardeen-Press-Teukolsky formula based on black hole spin. For Schwarzschild ($a=0$), ISCO = 3 $r_s$; varies from 0.5 $r_s$ (prograde, maximal spin) to 4.5 $r_s$ (retrograde).
* **Keplerian disk velocity:** $v_\phi = 1/\sqrt{2(r-1)}$ in Schwarzschild units.

The Kerr (spinning black hole) frame-dragging effect is approximated rather than fully derived from the Kerr metric, which would require implementing the Carter constant and Boyer-Lindquist coordinate geodesics.

### Interstellar-style controls

The GUI now includes controls matching the visual effects discussed in the Interstellar black-hole breakdown:

* **Accretion disk > temperature (K):** directly adjusts disk black-body color over the 4,500 K to 30,000 K range (orange to blue-white).
 * **Black hole > rotating shadow:** enables/disables spin-induced asymmetry of the black-hole shadow.
 * **Black hole > a/M:** sets the dimensionless spin parameter.
 * **Black hole > shadow squeeze:** scales the spin deformation strength (use lower values for a film-like reduced effect).
 * **Relativistic effects > doppler shift / beaming:** controls blue-shifted brightening of approaching plasma and red-shifted dimming of receding plasma.
 * **Relativistic effects > physical (D³ Liouville):** toggles between physically accurate D³ beaming (dramatic asymmetry) and softened cinematic beaming.

### System requirements

The simulation needs a decent GPU and a recent variant of Chrome or Firefox to run smoothly. In addition to changing simulation quality from the GUI, frame rate can be increased by shrinking the browser window and/or reducing screen resolution. Disabling the planet from the GUI also increases frame rate.

Example: runs 30+ fps at resolution 1920 x 1080 in Chrome 48 on a Linux desktop with GeForce GTX 750 Ti and "high" simulation quality

### Known artefacts

 * The accretion disk is stylized and procedurally textured for visual turbulence; it is not a full GRMHD/volumetric simulation.
 * The spectrum used in modeling the Doppler shift of the Milky Way background image is quite arbitrary (not based on real spectral data) and consequently the Doppler-shifted background colors may be wrong.
 * The lighting model of the planet is based on a point-like light source and a quite unphysical ambient component.
 * In the "medium" quality mode, the planet deforms unphysically when it travels between the camera and the black hole.
 * Adaptive step sizing reduces light path bending errors near the photon sphere (r = 1.5 $r_s$), but some inaccuracy remains at lower quality settings.
 * The Kerr frame-dragging is an approximation; true Kerr geodesics would require the Carter constant formalism.
 * Lorentz contraction causes jagged looks in the planet when simultaneously enabled with "light travel time" and the planet is close to the black hole.
 * Texture sampling issues cause unintended star blinking.

_see **[COPYRIGHT.md](https://github.com/oseiskar/black-hole/blob/master/COPYRIGHT.md)** for license and copyright info_
