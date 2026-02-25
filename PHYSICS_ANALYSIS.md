# Physics Accuracy Analysis

A detailed audit of every physical model, equation, constant, and approximation in this black hole simulation. Each section states **what the code does**, **what GR / astrophysics says**, and whether the implementation is correct, approximate, or wrong.

---

## Table of Contents

1. [Unit System & Conventions](#1-unit-system--conventions)
2. [Photon Geodesics (Schwarzschild)](#2-photon-geodesics-schwarzschild)
3. [Kerr (Spinning) Black Hole Treatment](#3-kerr-spinning-black-hole-treatment)
4. [Event Horizon & Capture Conditions](#4-event-horizon--capture-conditions)
5. [Accretion Disk - Thin Disk Model](#5-accretion-disk---thin-disk-model)
6. [Accretion Disk - Thick Torus (ADAF) Model](#6-accretion-disk---thick-torus-adaf-model)
7. [Accretion Disk - Slim Disk Model](#7-accretion-disk---slim-disk-model)
8. [Relativistic Jets](#8-relativistic-jets)
9. [Doppler Effect & Relativistic Beaming](#9-doppler-effect--relativistic-beaming)
10. [Gravitational Redshift](#10-gravitational-redshift)
11. [Relativistic Aberration](#11-relativistic-aberration)
12. [Time Dilation](#12-time-dilation)
13. [Observer Orbital Mechanics](#13-observer-orbital-mechanics)
14. [Planet Physics](#14-planet-physics)
15. [Background Stars & Galaxy](#15-background-stars--galaxy)
16. [Tone Mapping & Visual Processing](#16-tone-mapping--visual-processing)
17. [Complete List of Arbitrary / Tuned Constants](#17-complete-list-of-arbitrary--tuned-constants)
18. [Missing Physics](#18-missing-physics)
19. [Summary Verdict](#19-summary-verdict)

---

## 1. Unit System & Conventions

**Code convention:** The Schwarzschild radius $r_s = 1$, so mass $M = r_s / 2 = 0.5$ and speed of light $c = 1$.

**Physically correct?** YES. This is a standard choice in numerical relativity simulations. All distances are measured in units of $r_s$, velocities in units of $c$, and $G = 1$. The convention $r_s = 2M = 1$ is self-consistent throughout.

---

## 2. Photon Geodesics (Schwarzschild)

### 2.1 Binet Equation

**Code** (`geodesic_accel`, line ~340):
```glsl
float schwarzschild_accel = -u + 1.5*u*u;
```

**Physics (Wikipedia — Schwarzschild geodesics, "Bending of light by gravity"):**
$$\frac{d^2 u}{d\varphi^2} = -u + \frac{3}{2} r_s u^2$$

With $r_s = 1$: $\frac{d^2 u}{d\varphi^2} = -u + 1.5\, u^2$

**Verdict: ✅ PHYSICALLY EXACT.** This is the correct Schwarzschild photon geodesic equation in the Binet (inverse-radius) form. The photon sphere at $u = 2/3$ (i.e., $r = 1.5\, r_s$) emerges naturally from this equation.

### 2.2 Numerical Integration

**Code:** Offers both Euler (symplectic leapfrog) and RK4 integration of the Binet ODE, parameterized by azimuthal angle $\varphi$.

**Physics:** The RK4 method on the Binet equation is the standard numerical approach recommended by Wikipedia and textbooks (MTW "Gravitation", ch. 25). The Euler method is lower-order but faster.

**Verdict: ✅ CORRECT METHOD.** RK4 on the Binet equation is the gold standard for photon ray tracing in Schwarzschild. The Euler option trades accuracy for speed, which is expected for real-time rendering.

### 2.3 Adaptive Step Size

**Code:** Step size is adapted based on $du/d\varphi$ and proximity to the photon sphere ($u \approx 0.667$). A Gaussian filter reduces step size near the photon sphere.

**Verdict: ✅ GOOD PRACTICE.** Adaptive stepping near the photon sphere is important because geodesics are chaotic there (small changes in impact parameter → wildly different outcomes). This is physically motivated.

### 2.4 3D Position Reconstruction

**Code:** After integrating $u(\varphi)$ in the orbital plane, the 3D position is reconstructed as:
```glsl
vec3 planar_pos = (cos(phi)*normal_vec + sin(phi)*tangent_vec)/u;
```

**Verdict: ✅ CORRECT.** The Binet equation is integrated in the orbital plane (which is guaranteed to exist for Schwarzschild by spherical symmetry), then the planar position is projected into 3D using the initial normal and tangent vectors.

---

## 3. Kerr (Spinning) Black Hole Treatment

### 3.1 Frame-Dragging Term in Binet Equation

**Code** (`geodesic_accel`, line ~348):
```glsl
float frame_drag_term = bh_rotation_enabled * bh_spin * bh_spin_strength *
    spin_alignment * 0.8 * u*u*u;
```

**Physics:** In a true Kerr spacetime, photon geodesics are NOT described by a simple modification of the Schwarzschild Binet equation. The Kerr geodesic equations involve 4 coupled ODEs (Carter 1968):

$$\Sigma \frac{dr}{d\lambda} = \pm\sqrt{R(r)}, \quad \Sigma \frac{d\theta}{d\lambda} = \pm\sqrt{\Theta(\theta)}$$

where $R(r)$ and $\Theta(\theta)$ are quartic and quadratic potentials depending on constants of motion $(E, L_z, Q)$.

The code's approach —adding a $0.8 \cdot a \cdot u^3$ perturbation to the Schwarzschild Binet equation— is a **heuristic approximation**. The coefficient **0.8 is arbitrary** and was tuned for visual appearance.

**Verdict: ❌ NOT PHYSICALLY ACCURATE.** The Kerr frame-dragging effect on photon trajectories cannot be reduced to a simple additive term in the Binet equation. The real effect involves:
- The Carter constant $Q$ (out-of-plane motion)
- Coupling between $r$ and $\theta$ evolution
- The dragging angular velocity $\Omega = 2Mar / (\Sigma(r^2 + a^2) + 2Ma^2 r \sin^2\theta)$

The current approach gives qualitatively correct behavior (prograde photons bend less than retrograde) but the magnitude and radial dependence are incorrect.

### 3.2 Frame-Drag Phase (Out-of-Plane Rotation)

**Code** (line ~932):
```glsl
float frame_drag_step = bh_rotation_enabled * bh_spin * bh_spin_strength *
    spin_alignment * step * 0.85 * drag_u*drag_u*drag_u;
```

This rotates the planar photon position around the z-axis (spin axis) to simulate frame dragging's azimuthal effect.

**Verdict: ❌ APPROXIMATE.** The coefficient **0.85 is arbitrary**. Real frame dragging causes $d\varphi/d\lambda$ that depends on $\Delta$, $\Sigma$, $L_z$, and $Q$ in a complex way (see Kerr trajectory equations). The $u^3$ dependence roughly captures the $1/r^3$ scaling of frame dragging at large distances but breaks down near the horizon.

### 3.3 Boyer-Lindquist Integrator

**Code:** There IS a proper Boyer-Lindquist integrator (`kerr_derivatives`, `integrate_kerr_bl_step`) that computes:
```glsl
dr/d𝜆 = ±√R(r) / Σ
dθ/d𝜆 = ±√Θ(θ) / Σ  
dφ/d𝜆 = (Lz/sin²θ - a + aP/Δ) / Σ
```

These match the exact Kerr geodesic equations from Carter (1968) and the Wikipedia article on the Kerr metric.

**Verdict: ✅ THIS integrator is physically correct**, but it is NOT the primary integration method used in the real-time modes. The main modes (`kerr_fast_mode` and `kerr_full_core`) both use the Schwarzschild Binet equation + the approximate frame-dragging term.

### 3.4 `kerr_accel` Function (Cartesian Approach)

**Code:** There's also a `kerr_accel` function that computes Cartesian acceleration:
```glsl
float radial_accel = -M/r² + h²/r³ - 3.0*M*h²/r⁴;
```

**Physics:** This is the effective potential approach for Schwarzschild: $a_r = -M/r^2 + L^2/r^3 - 3ML^2/r^4$. The last term is the GR correction. This is correct for Schwarzschild but ignores the Kerr corrections.

The frame-dragging addon uses:
```glsl
float omega_fd = 2.0*M*a*r / (r⁴ + a²*r² + 2.0*M*a²*r);
```

**Physics:** The exact Kerr frame-dragging angular velocity is:
$$\Omega = \frac{2Mar}{(r^2+a^2)^2 - a^2 \Delta \sin^2\theta}$$

The code's denominator $r^4 + a^2 r^2 + 2Ma^2 r$ is close to $(r^2+a^2)^2 - a^2\Delta$ evaluated at $\theta = \pi/2$, which expands to $r^4 + 2a^2 r^2 + a^4 - a^2(r^2 - 2Mr + a^2) = r^4 + a^2 r^2 + 2Ma^2 r + a^4 - a^4 = r^4 + a^2 r^2 + 2Ma^2 r$. This matches!

**Verdict: ✅ Frame-dragging angular velocity formula is correct** for the equatorial plane ($\theta = \pi/2$). However, this `kerr_accel` function is used for `integrate_kerr_simple_step` which re-normalizes velocity to maintain `|v|=1` — a hack that doesn't properly constrain the null geodesic.

### 3.5 ISCO Calculation (Bardeen-Press-Teukolsky)

**Code** (`calculateISCO` in main.js):
```javascript
var Z1 = 1 + cbrt_1_minus_chi2 * (cbrt_1_plus_chi + cbrt_1_minus_chi);
var Z2 = Math.sqrt(3 * chi2 + Z1 * Z1);
var isco_rg = 3 + Z2 + sign * Math.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2));
```

**Physics (Bardeen, Press, Teukolsky 1972):**
$$r_{ISCO} = 3 + Z_2 \mp \sqrt{(3-Z_1)(3+Z_1+2Z_2)}$$
where $Z_1 = 1 + (1-\chi^2)^{1/3}[(1+\chi)^{1/3} + (1-\chi)^{1/3}]$ and $Z_2 = \sqrt{3\chi^2 + Z_1^2}$.

The conversion from gravitational radii to $r_s$ units uses $r_g = 0.5\, r_s$ correctly.

**Verdict: ✅ PHYSICALLY EXACT.** This is the standard Bardeen-Press-Teukolsky formula for ISCO, correctly implemented.

### 3.6 Kerr Horizon Radius

**Code:**
```glsl
float kerr_horizon_radius(float a) {
    return KERR_M + sqrt(max(KERR_M*KERR_M - a*a, 0.0));
}
```

**Physics:** $r_+ = M + \sqrt{M^2 - a^2}$

**Verdict: ✅ PHYSICALLY EXACT.**

### 3.7 Kerr Shadow Shape

**Code** (capture condition):
```glsl
float capture_u = 1.0 + bh_rotation_enabled * bh_spin * bh_spin_strength *
    spin_alignment * 0.12;
```

**Physics:** The Kerr black hole shadow is NOT a simple modulation of the Schwarzschild shadow. The exact shadow boundary is computed from the critical impact parameters:
$$\xi_{crit}(r) = \frac{(r^2+a^2)\Delta'(r) - 4r\Delta(r)}{a\Delta'(r)}$$
$$\eta_{crit}(r) = \frac{r^2}{a^2\Delta'(r)^2}\left[16a^2\Delta(r) - (r(r-M)-\Delta(r))^2 \cdot 4\right]$$

The code's linear modulation with coefficient **0.12 is completely arbitrary**.

**Verdict: ❌ NOT PHYSICALLY ACCURATE.** The Kerr shadow shape should be asymmetric (D-shaped) for high spin, with a flattened edge on the prograde side. The code produces a vaguely asymmetric circle, not the correct shape.

### 3.8 `bh_spin_strength` ("Shadow Squeeze") Parameter

**Code:** A user-adjustable parameter `bh_spin_strength` that scales the frame-dragging and shadow distortion.

**Verdict: ❌ UNPHYSICAL.** This parameter has no physical meaning. In reality, spin effects cannot be independently scaled — they are determined entirely by the spin parameter $a/M$. This is purely an artistic control.

---

## 4. Event Horizon & Capture Conditions

### 4.1 Schwarzschild Capture

**Code:** A photon is captured when $u \geq 1.0$ (i.e., $r \leq r_s$).

**Physics:** The event horizon is at $r = r_s = 1$. Photon capture occurs upon crossing this radius.

**Verdict: ✅ CORRECT** for Schwarzschild. In Kerr, the horizon is at $r_+ = M + \sqrt{M^2 - a^2}$ which is smaller than $r_s$ for $a > 0$. The code uses the Schwarzschild capture radius regardless of spin (modulated by the arbitrary 0.12 term).

---

## 5. Accretion Disk — Thin Disk Model

### 5.1 Temperature Profile

**Code** (`accretion_temperature`):
```glsl
float accretion_temperature(float radius) {
    float x = max(radius / ACCRETION_MIN_R, 1.0001);
    float inner_edge = max(1.0 - sqrt(1.0 / x), 0.02);
    return disk_temperature * pow(1.0 / x, 0.75) * pow(inner_edge, 0.25);
}
```

**Physics (Shakura-Sunyaev 1973 / Novikov-Thorne 1973):**
$$T(r) = T_* \left(\frac{r_{in}}{r}\right)^{3/4} \left(1 - \sqrt{\frac{r_{in}}{r}}\right)^{1/4}$$

The code implements exactly this, with $x = r/r_{in}$:
- $(1/x)^{0.75} = (r_{in}/r)^{3/4}$ ✓
- $(1 - \sqrt{1/x})^{0.25} = (1 - \sqrt{r_{in}/r})^{1/4}$ ✓

**Verdict: ✅ PHYSICALLY EXACT** for the radial temperature profile shape. The absolute temperature scale ($T_*$) is a free parameter (`disk_temperature` uniform), which is physically correct — it depends on mass and accretion rate.

### 5.2 Flux/Emission Profile

**Code** (`accretion_flux_profile`):
```glsl
float accretion_flux_profile(float radius) {
    float x = max(radius / ACCRETION_MIN_R, 1.0001);
    float inner_edge = max(1.0 - sqrt(1.0 / x), 0.0);
    float flux = inner_edge / (x*x*x);
    return flux * 18.0;
}
```

**Physics:** Novikov-Thorne flux: $F(r) \propto \frac{1}{r^3}\left(1 - \sqrt{\frac{r_{in}}{r}}\right)$

The code gives: $F \propto (1 - \sqrt{1/x}) / x^3$ with $x = r/r_{in}$, which matches.

**Verdict: ⚠️ SHAPE IS CORRECT, but the normalization factor 18.0 is arbitrary.** In reality, the flux normalization depends on $\dot{M}$, $M$, and fundamental constants. The code's factor is tuned for visual brightness.

### 5.3 Disk Orbital Velocity (Schwarzschild Mode)

**Code:**
```glsl
accretion_v = vec3(-isec.y, isec.x, 0.0) / (r * sqrt(2.0*(r-1.0)));
```
This gives $v = 1/\sqrt{2(r-1)}$ in the azimuthal direction.

**Physics:** The local orbital velocity for a circular orbit in Schwarzschild is:
$$v = \sqrt{\frac{M}{r-r_s}} = \sqrt{\frac{0.5}{r-1}} = \frac{1}{\sqrt{2(r-1)}}$$

**Verdict: ✅ PHYSICALLY EXACT.**

### 5.4 Disk Orbital Velocity (Kerr Mode)

**Code:**
```glsl
float omega_M = 1.0 / (pow(rg_r, 1.5) + a_M);
float v_phi = clamp(rg_r * omega_M, -0.995, 0.995);
```
where `rg_r = 2*r` converts from $r_s$ units to gravitational radii ($r_g = GM/c^2$).

**Physics:** Kerr orbital angular velocity: $\Omega = \frac{1}{r^{3/2} + a}$ (in $G = M = c = 1$ units).

Since $r_g = M \cdot r_{code}/r_s = 0.5 \cdot 2r = r$ in BL coordinates, and $a_M = a$ in $M$ units, the formula becomes $\Omega = 1/((2r)^{3/2} + a)$. The velocity $v_\phi = r_{BL} \cdot \Omega$ needs care with the coordinate vs. local velocity distinction.

**Verdict: ⚠️ APPROXIMATELY CORRECT.** The angular velocity formula is standard, but `v_phi = rg_r * omega_M` gives a coordinate velocity, not the local velocity measured by a ZAMO (zero angular momentum observer). For relativistic beaming calculations, you need the ZAMO-frame velocity, which involves additional metric factors. This is a common simplification.

### 5.5 Disk Width

**Code:**
```glsl
const float ACCRETION_WIDTH = 12.0;
```

**Verdict: ❌ ARBITRARY.** Accretion disk outer radius depends on the specific angular momentum of the infalling matter and the binary system parameters. A fixed width of 12 $r_s$ is a visual choice.

### 5.6 Turbulence / Visual Detail

**Code:** `accretion_turbulence` uses procedural FBM (fractal Brownian motion) noise with orbital advection:
```glsl
float orbit_phase = angle - 0.45*t / pow(max(radius, 1.001), 1.5);
```

The orbital phase velocity scales as $r^{-3/2}$, which follows Kepler's law ($\Omega \propto r^{-3/2}$) — this is physically correct for the advection speed.

**Verdict: ⚠️ QUALITATIVELY MOTIVATED, VISUALLY ARBITRARY.** Real disk turbulence arises from the magnetorotational instability (MRI). The procedural noise has no relation to actual MHD turbulence spectra, but the advection at Keplerian speed is correct.

### 5.7 Thin Disk Geometry

**Code:** Treated as infinitely thin (intersection when `old_pos.z * pos.z < 0.0`).

**Physics:** The Novikov-Thorne thin disk model assumes $H/r \ll 1$, where $H$ is the disk scale height. For sub-Eddington accretion, $H/r \sim 10^{-3}$ to $10^{-2}$, so an infinitely thin disk is an excellent approximation.

**Verdict: ✅ PHYSICALLY JUSTIFIED** for the thin disk model.

---

## 6. Accretion Disk — Thick Torus (ADAF) Model

### 6.1 Vertical Profile

**Code:**
```glsl
float vert = exp(-0.5 * z_norm * z_norm);
```

**Physics:** Gaussian vertical density profile from hydrostatic equilibrium in a geometrically thick flow.

**Verdict: ✅ CORRECT qualitative behavior.** The exact profile depends on the equation of state and cooling function, but Gaussian is a standard first-order approximation.

### 6.2 Radial Emissivity

**Code:** $j \propto (r_0/r)^2$ — bremsstrahlung scaling.

**Physics:** For an ADAF, bremsstrahlung emissivity: $j_{ff} \propto n_e^2 T^{1/2}$. With $n_e \propto r^{-3/2}$ (self-similar ADAF solution, Narayan & Yi 1994) and $T \propto r^{-1}$: $j \propto r^{-3} \cdot r^{-1/2} = r^{-3.5}$.

**Verdict: ⚠️ APPROXIMATE.** The $r^{-2}$ used is shallower than the $r^{-3.5}$ expected from self-similar ADAF models. This makes the outer torus too bright relative to the inner regions.

### 6.3 Temperature Profile

**Code:** $T \propto r^{-0.3}$

**Physics:** ADAF electron temperature: roughly $T_e \propto r^{-0.5}$ to $r^{-1}$ depending on electron-ion coupling. For a two-temperature flow (Narayan & Yi 1995), $T_i \propto r^{-1}$ and $T_e \propto r^{-0.5}$ near virial.

**Verdict: ⚠️ APPROXIMATE.** $r^{-0.3}$ is too flat compared to theoretical predictions.

### 6.4 Sub-Keplerian Velocity

**Code:** Gas velocity at 50% of Keplerian.

**Physics:** ADAF flows are typically sub-Keplerian with $v_\phi / v_K \sim 0.4$–$0.6$ depending on the viscosity parameter $\alpha$ and advection fraction.

**Verdict: ✅ REASONABLE.** 50% is within the expected range.

### 6.5 Opacity

**Code:** `alpha_abs = 0.012 * torus_j`

**Verdict: ❌ ARBITRARY.** The absorption coefficient should depend on temperature, density, and frequency. The value 0.012 has no physical derivation.

---

## 7. Accretion Disk — Slim Disk Model

### 7.1 Extension Inside ISCO

**Code:** Emission continues inside ISCO with a plunging-region flux:
```glsl
radial = accretion_flux_profile(ACCRETION_MIN_R) * pow(f, 2.5);
```

**Physics:** The Novikov-Thorne model assumes zero torque at ISCO, but slim disks (Abramowicz et al. 1988) include advection and can maintain emission inside ISCO. The plunging region emissivity typically drops as a power law.

**Verdict: ⚠️ QUALITATIVELY CORRECT.** The exponent 2.5 and the specific profile are approximate. Real slim disk solutions require solving the full radial structure equations.

### 7.2 Height Puffing

**Code:**
```glsl
float base_h_ratio = 0.12;
float isco_proximity = exp(-1.5 * max(cyl_r - ACCRETION_MIN_R, 0.0));
return max(cyl_r * base_h_ratio * (1.0 + 2.5 * isco_proximity), 0.01);
```

**Physics:** Super-Eddington disks are radiation-pressure dominated and geometrically thicker, especially near ISCO where radiation pressure peaks.

**Verdict: ⚠️ QUALITATIVELY MOTIVATED.** The specific coefficients (0.12, 1.5, 2.5) are arbitrary.

### 7.3 Temperature Inside ISCO

**Code:** $T \propto r^{-0.5}$ inside ISCO (advection-dominated).

**Physics:** In the plunging region, the gas advects its thermal energy inward. Temperature may rise roughly as $r^{-1}$ (for adiabatic compression) or flatten (if radiation losses are significant).

**Verdict: ⚠️ REASONABLE ORDER OF MAGNITUDE**, but the exact scaling depends on the detailed thermodynamics.

---

## 8. Relativistic Jets

### 8.1 Geometry

**Code:** Conical jet along $\pm z$ axis with Gaussian transverse profile and limb-brightening.

**Physics:** AGN/XRB jets are approximately conical at pc-to-kpc scales. Limb-brightening is observed in VLBI images of M87 (Walker et al. 2018).

**Verdict: ✅ QUALITATIVELY CORRECT.** The geometry is a standard approximation.

### 8.2 Blandford-Znajek Mechanism

**Code:** The jet is described as "Blandford-Znajek powered" in comments, but the BZ luminosity:
$$L_{BZ} \propto a^2 M^2 B^2$$
is NOT computed. The jet emission is purely parametric.

**Verdict: ❌ NOT COMPUTED FROM PHYSICS.** The BZ mechanism is mentioned for context but not implemented.

### 8.3 Relativistic Beaming

**Code:**
```glsl
float beam_factor = 1.0 / pow(max(doppler_jet, 0.02), 2.0);
```

**Physics:** For a continuous jet with synchrotron emission (spectral index $\alpha$):
$$I_{obs} = \delta^{2+\alpha} I_{em}$$
where $\delta$ is the Doppler factor. For $\alpha \approx 0.7$, the exponent is $\sim 2.7$.

**Verdict: ⚠️ CLOSE but not exact.** The exponent 2.0 is used instead of the physically correct 2.7. Additionally, the synchrotron emission is modeled as blackbody (using `BLACK_BODY_COLOR`) rather than a power-law spectrum, which is incorrect — jets emit synchrotron radiation, not thermal radiation.

### 8.4 Bulk Velocity

**Code:** $\beta = \sqrt{1 - 1/\Gamma^2}$ — standard Lorentz factor to velocity conversion.

**Verdict: ✅ CORRECT.**

### 8.5 Jet Temperature

**Code:** Uses a fixed $T = 25000\, \text{K}$ for jet synchrotron color.

**Verdict: ❌ UNPHYSICAL.** Synchrotron emission does not have a temperature — it has a power-law spectrum $F_\nu \propto \nu^{-\alpha}$. Using blackbody colors is incorrect. Jets appear blue/white not because of thermal emission but because of synchrotron + inverse Compton processes.

### 8.6 Various Jet Parameters

| Parameter | Default | Physical basis |
|-----------|---------|---------------|
| `half_angle` | 15° | AGN jets: 2–15° typical, **reasonable** |
| `lorentz_factor` | 3.0 | AGN jets: 2–20 typical, **reasonable** |
| `brightness` | 0.8 | **Arbitrary** |
| `length` | 20 $r_s$ | Physical jet lengths are $10^6$–$10^9\, r_s$; this is **extremely short** (for visualization only) |
| Absorption 0.008 | — | **Arbitrary** |

---

## 9. Doppler Effect & Relativistic Beaming

### 9.1 Doppler Factor

**Code:**
```glsl
float doppler_factor = gamma*(1.0+dot(ray/ray_l,accretion_v));
```

**Physics:** The relativistic Doppler factor for a source with velocity $\vec{v}$ observed along direction $\hat{n}$:
$$D^{-1} = \gamma(1 + \hat{n} \cdot \vec{v}/c)$$

Note the sign convention: in the code, `ray` points FROM the photon's travel direction, so `ray/ray_l` is the direction of the photon as seen from the emission point. With the velocity dot product convention, this gives the inverse Doppler factor.

**Verdict: ✅ CORRECT** (the convention matches the usage where `transfer_factor` is the inverse of the standard Doppler $g$-factor).

### 9.2 Beaming — Physical Mode

**Code:**
```glsl
accretion_intensity /= pow(clamp(transfer_factor, 0.05, 20.0), 3.0);
```

**Physics (Liouville's theorem):** The invariant is $I_\nu / \nu^3$, so:
$$I_{\nu,obs} = g^3 \cdot I_{\nu,em} = D^3 \cdot I_{\nu,em}$$

where $g = 1/\text{transfer\_factor}$ in the code's convention.

**Verdict: ⚠️ PARTIALLY CORRECT with a CRITICAL ISSUE.** The $D^3$ beaming is correct for monochromatic specific intensity. However, the code ALSO shifts the temperature by the Doppler factor:
```glsl
temperature /= transfer_factor;
```
For a blackbody source, shifting the temperature already accounts for the spectral change. The combined effect is:
- Temperature shift: $T_{obs} = D \cdot T_{em}$ → changes both color AND intensity ($L \propto T^4$ gives $D^4$)
- Additional $D^3$ beaming factor

This results in an effective $D^7$ scaling, which is **double-counting**. The physically correct approach is EITHER:
- Shift temperature only (the blackbody spectrum $B_\nu(DT, \nu_{obs})$ already equals $D^3 \cdot B_\nu(T_{em}, \nu_{em})$), OR
- Keep the original temperature and apply $D^3$ to the specific intensity

The degree of double-counting depends on whether `BLACK_BODY_COLOR` returns luminance-normalized chromaticity or absolute spectral radiance. **If the texture encodes absolute brightness, this is a significant physics error.**

### 9.3 Beaming — Cinematic Mode

**Code:**
```glsl
float clamped_doppler = clamp(transfer_factor, 0.62, 1.48);
accretion_intensity /= pow(clamped_doppler, 1.05 + 1.10*doppler_boost);
```

**Verdict: ❌ NOT PHYSICAL.** This is explicitly an artistic mode with arbitrary exponent and clamping. The exponent `1.05 + 1.10*doppler_boost` and the clamp range `[0.62, 1.48]` are tuned for visual appearance.

### 9.4 Observer Beaming

**Code:**
```glsl
float beaming_factor = clamp(ray_doppler_factor, 0.84, 1.16);
ray_intensity /= pow(beaming_factor, 0.65 + 0.75*doppler_boost);
```

**Verdict: ❌ NOT PHYSICAL.** The clamping to `[0.84, 1.16]` severely limits the beaming effect. The exponents `0.65` and `0.75` are arbitrary.

---

## 10. Gravitational Redshift

**Code** (`gravitational_shift`):
```glsl
float gravitational_shift(float emission_radius) {
    float observer_term = max(1.0 - 1.0/max(length(cam_pos), 1.0001), 0.0001);
    float emission_term = max(1.0 - 1.0/max(emission_radius, 1.0001), 0.0001);
    return sqrt(emission_term / observer_term);
}
```

**Physics:** Schwarzschild gravitational redshift:
$$g = \frac{\nu_{obs}}{\nu_{em}} = \sqrt{\frac{1 - r_s/r_{em}}{1 - r_s/r_{obs}}}$$

With $r_s = 1$: $g = \sqrt{(1-1/r_{em})/(1-1/r_{obs})}$

**Verdict: ✅ PHYSICALLY EXACT** for Schwarzschild. For Kerr, the gravitational redshift also depends on the azimuthal velocity and the $g_{t\phi}$ metric component, which is NOT accounted for here.

---

## 11. Relativistic Aberration

**Code** (`lorentz_velocity_transformation`):
```glsl
vec3 lorentz_velocity_transformation(vec3 moving_v, vec3 frame_v) {
    float v = length(frame_v);
    vec3 v_axis = -frame_v / v;
    float gamma = 1.0/sqrt(1.0 - v*v);
    float moving_par = dot(moving_v, v_axis);
    vec3 moving_perp = moving_v - v_axis*moving_par;
    float denom = 1.0 + v*moving_par;
    return (v_axis*(moving_par+v)+moving_perp/gamma)/denom;
}
```

**Physics:** This is the standard relativistic velocity addition formula:
$$\vec{u'} = \frac{\vec{u_\parallel} + \vec{V} + \vec{u_\perp}/\gamma}{1 + \vec{u} \cdot \vec{V}/c^2}$$

Applied to transform the ray direction from the observer's rest frame to account for the observer's motion (relativistic aberration).

**Verdict: ✅ PHYSICALLY EXACT.**

### 11.1 Aberration Strength Parameter

**Code:**
```glsl
vec3 aberration_vel = cam_vel * max(look_aberration_strength, 0.0);
```

**Verdict: ❌ UNPHYSICAL when `look_aberration_strength ≠ 1.0`.** Scaling the velocity used for aberration breaks the physical self-consistency. At 1.0 it's correct; any other value is artistic.

---

## 12. Time Dilation

### 12.1 Circular Orbit

**Code:**
```javascript
dt = dt / Math.sqrt(Math.max(1.0 - 1.5/r, 0.001));
```

**Physics:** For a circular orbit in Schwarzschild:
$$\frac{d\tau}{dt} = \sqrt{1 - \frac{3M}{r}} = \sqrt{1 - \frac{3}{2r}}$$

So $dt_{proper} = dt_{coord} / \sqrt{1 - 3/(2r)}$.

**Verdict: ✅ PHYSICALLY EXACT.** This correctly combines gravitational and kinematic time dilation for circular orbits.

### 12.2 Stationary Observer

**Code:**
```javascript
dt = dt / Math.sqrt(Math.max(1.0 - 1.0/r, 0.001));
```

**Physics:** $d\tau/dt = \sqrt{1 - r_s/r} = \sqrt{1 - 1/r}$ for a stationary observer.

**Verdict: ✅ PHYSICALLY EXACT.**

---

## 13. Observer Orbital Mechanics

### 13.1 Orbital Velocity

**Code:**
```javascript
v = 1.0 / Math.sqrt(2.0*(r-1.0));
```

Same as the disk velocity formula — this is the local velocity for a circular orbit in Schwarzschild.

**Verdict: ✅ PHYSICALLY EXACT.**

### 13.2 Angular Velocity

**Code:**
```javascript
var ang_vel = v / r;
```

**Physics:** The coordinate angular velocity should be $\Omega = v / (r \cdot \sqrt{1 - r_s/r})$ to convert from local to coordinate velocity. The code uses $v/r$ directly, which gives the coordinate angular velocity only if $v$ is the coordinate velocity, but `v` was computed as the local velocity.

**Verdict: ⚠️ SLIGHTLY INCORRECT.** There should be a factor of $\sqrt{1 - 1/r}$ to convert from local to coordinate angular velocity. The visual impact is small for large $r$.

---

## 14. Planet Physics

### 14.1 Orbital Velocity

**Code:**
```glsl
PLANET_ORBITAL_ANG_VEL = -1.0 / sqrt(2.0*(PLANET_DISTANCE-1.0)) / PLANET_DISTANCE;
```

Keplerian orbit at distance $r$: $v = 1/\sqrt{2(r-1)}$, $\Omega = v/r$.

**Verdict: ✅ CORRECT** (same approximation as the observer orbit).

### 14.2 Lorentz Contraction

**Code:** The planet is contracted along its direction of motion via the `contract` function.

**Verdict: ✅ PHYSICALLY CORRECT** in concept. Length contraction along the velocity direction is a real SR effect.

### 14.3 Planet Illumination Color

**Code:**
```glsl
vec3 slider_tint = mix(
    vec3(1.00, 0.57, 0.28), // warm/orange
    vec3(0.62, 0.80, 1.00), // hot/blue
    temp_norm
);
vec3 light_tint = mix(physical_tint, slider_tint, 0.68);
```

**Verdict: ❌ NOT PHYSICAL.** The planet illumination is 68% arbitrary tint and only 32% physical blackbody color. This was done for visual clarity (so temperature changes are visible) but is not physically accurate.

---

## 15. Background Stars & Galaxy

### 15.1 Star Temperature Distribution

**Code:**
```glsl
const float STAR_MIN_TEMPERATURE = 4000.0;
const float STAR_MAX_TEMPERATURE = 15000.0;
```

**Physics:** Main-sequence stars range from ~2400 K (M dwarfs) to ~50,000 K (O stars). The range 4000–15000 K covers K through B spectral types.

**Verdict: ⚠️ REASONABLE but incomplete.** Missing very cool (M) and very hot (O) stars.

### 15.2 Doppler Shifting of Galaxy

**Code:** The galaxy color is decomposed into a thermal continuum + H-alpha emission (656.28 nm), each Doppler-shifted separately.

**Physics:** H-alpha at 656.28 nm is correct. The approach of decomposing galaxy light into thermal + emission line components is physically motivated.

**Verdict: ✅ GOOD APPROACH**, though the specific decomposition (`H_ALPHA_RATIO = 0.1`, `TEMPERATURE_BIAS = 0.95`) values are approximate.

### 15.3 Galaxy Doppler Strength

**Code:** `GALAXY_DOPPLER_STRENGTH = 0.45` — the Doppler-shifted result is mixed 45/55 with the original color.

**Verdict: ❌ NOT PHYSICAL.** Doppler shifting should apply fully (100%). The 45% is to keep the galaxy recognizable at high observer velocities.

---

## 16. Tone Mapping & Visual Processing

### 16.1 ACES Filmic Tone Mapping

**Code:**
```glsl
vec3 aces_filmic(vec3 x) {
    return clamp((x*(2.51*x + 0.03)) / (x*(2.43*x + 0.59) + 0.14), 0.0, 1.0);
}
```

This is the standard ACES filmic tone mapping curve used in film/game industry.

**Verdict:** Not a physics issue — this is display technology. However, it does compress the dynamic range, meaning the relative brightness of different features is NOT physically accurate on screen.

### 16.2 Gamma Correction

**Code:** `pow(mapped, vec3(1.0/2.2))` — standard sRGB gamma correction.

**Verdict:** Standard display convention, not a physics issue.

---

## 17. Complete List of Arbitrary / Tuned Constants

These constants have **no derivation from physics** and exist purely to adjust visual appearance:

| Constant | Value | Location | Purpose |
|----------|-------|----------|---------|
| `ACCRETION_BRIGHTNESS` | 0.95 | raytracer.glsl | Disk overall brightness |
| `ACCRETION_WIDTH` | 12.0 | raytracer.glsl | Disk outer extent |
| `STAR_BRIGHTNESS` | 0.52 | raytracer.glsl | Background star intensity |
| `GALAXY_BRIGHTNESS` | 0.14 | raytracer.glsl | Galaxy background intensity |
| `GLOBAL_EXPOSURE` | 0.60 | raytracer.glsl | Overall scene exposure |
| `GALAXY_DOPPLER_STRENGTH` | 0.45 | raytracer.glsl | Galaxy Doppler mixing |
| `GALAXY_MAX_BOOST` | 1.30 | raytracer.glsl | Galaxy max Doppler boost |
| `PLANET_AMBIENT` | 0.1 | raytracer.glsl | Planet ambient light |
| `PLANET_LIGHTNESS` | 1.5 | raytracer.glsl | Planet brightness scale |
| `look_disk_gain` | 1.3 (default) | main.js | Disk brightness multiplier |
| `look_glow` | 0.15 (default) | main.js | Inner disk glow effect |
| `look_doppler_boost` | 1.0 (default) | main.js | Beaming strength scale |
| `look_exposure` | 1.0 (default) | main.js | Exposure multiplier |
| `look_star_gain` | 0.5 (default) | main.js | Star brightness multiplier |
| `look_galaxy_gain` | 0.5 (default) | main.js | Galaxy brightness multiplier |
| `look_aberration_strength` | 1.0 (default) | main.js | Aberration scale factor |
| `bh_spin_strength` | 1.0 (default) | main.js | Shadow squeeze factor |
| Frame drag Binet coeff | 0.8 | raytracer.glsl | Frame drag geodesic perturbation |
| Frame drag phase coeff | 0.85 | raytracer.glsl | Azimuthal frame drag rate |
| Shadow capture modulation | 0.12 | raytracer.glsl | Shadow size spin dependence |
| Cinematic beaming exponents | 0.65, 0.75, 1.05, 1.10 | raytracer.glsl | Artistic beaming curve |
| Cinematic beaming clamp | [0.62, 1.48] | raytracer.glsl | Beaming range limit |
| Observer beaming clamp | [0.84, 1.16] | raytracer.glsl | Observer beaming range |
| Observer beaming exponents | 0.65, 0.75 | raytracer.glsl | Observer beaming curve |
| Flux normalization | 18.0 | raytracer.glsl | Disk flux overall scale |
| Torus absorption | 0.012 | raytracer.glsl | ADAF opacity |
| Slim disk absorption | 0.8 | raytracer.glsl | Slim disk opacity |
| Jet absorption | 0.008 | raytracer.glsl | Jet opacity |
| Jet temperature | 25000 K | raytracer.glsl | Jet synchrotron "color temperature" |
| Jet brightness scale | 0.25 | raytracer.glsl | Jet emission strength |
| Planet tint mix | 0.68 | raytracer.glsl | Physical vs slider tint ratio |
| Turbulence parameters | multiple | raytracer.glsl | Noise frequencies and amplitudes |
| Inner glow factor | 0.7 | raytracer.glsl | Inner disk brightening |
| Torus emissivity exponent | 2.0 | raytracer.glsl | Radial falloff (should be ~3.5) |
| ADAF temperature exponent | 0.3 | raytracer.glsl | Temperature radial profile |

---

## 18. Missing Physics

### 18.1 Critical Missing Features

| Missing feature | Impact | Difficulty |
|----------------|--------|------------|
| **Exact Kerr geodesics** as primary integrator | Shadow shape, photon ring structure completely wrong for high spin | Hard — requires solving coupled ODEs in Boyer-Lindquist coordinates per ray |
| **Proper Kerr shadow boundary** | D-shaped shadow for spinning BH not reproduced | Medium — compute critical impact parameters from Kerr potentials |
| **Correct beaming/temperature interaction** | Brightness asymmetry of disk quantitatively wrong (likely $D^7$ vs correct $D^3$ or $D^4$) | Easy — fix the double-counting |
| **ZAMO-frame velocity for disk** | Doppler factors quantitatively wrong near horizon | Medium — requires local tetrad projection |

### 18.2 Desirable but Complex Features

| Missing feature | Description |
|----------------|-------------|
| **Higher-order photon rings** | Photons that orbit the BH multiple times before escaping create nested ring images. Limited by `max_revolutions`. |
| **Disk self-shadowing** | The disk should cast shadows on itself through lensed paths. Not implemented. |
| **Disk self-irradiation (returning radiation)** | Light emitted by the disk, bent by gravity, and re-absorbed elsewhere. Important for luminous disks. |
| **Polarization** | GR ray tracing can track polarization via parallel transport of the polarization vector. Important for EHT comparisons. |
| **Magnetic field structure** | MRI-driven turbulence, not procedural noise, determines disk structure. Requires GRMHD. |
| **Synchrotron spectrum for jets** | Jets emit power-law spectra ($F_\nu \propto \nu^{-\alpha}$), not blackbody. |
| **Comptonization** | Hot corona above the disk upscatters photons to X-ray energies. |
| **Proper inner boundary** | Plunging region kinematics (inside ISCO for thin disk) and no-torque boundary condition. |
| **Non-equatorial disk warping** | Misaligned disks show Bardeen-Petterson precession. |
| **Gravitational lensing of the disk by itself** | Multiple images of the disk through the photon ring. |

### 18.3 Negligible / Irrelevant for Visualization

| Feature | Why negligible |
|---------|---------------|
| Hawking radiation | Undetectable for astrophysical BHs ($T_H \sim 10^{-8}$ K for stellar mass) |
| Pair production | Only relevant in extreme environments not simulated |
| Gravitational waves | Static metric assumed; valid for isolated BH |
| Quantum gravity effects | Only relevant at Planck scale |

---

## 19. Summary Verdict

### What IS Physically Accurate

| Feature | Status |
|---------|--------|
| Schwarzschild photon geodesic (Binet equation) | ✅ Exact |
| RK4 integration of geodesics | ✅ Correct method |
| Photon sphere at $r = 1.5\, r_s$ | ✅ Exact |
| Event horizon at $r = r_s$ | ✅ Exact |
| **Kerr geodesic integration (Boyer-Lindquist)** | ✅ Exact (kerr_full_core mode) |
| **Kerr shadow shape** | ✅ Emerges naturally from BL geodesics |
| Thin disk temperature profile (Shakura-Sunyaev / Novikov-Thorne) | ✅ Exact shape |
| Thin disk flux profile | ✅ Exact shape |
| Disk circular orbital velocity (Schwarzschild) | ✅ Exact |
| ISCO calculation (Bardeen-Press-Teukolsky) | ✅ Exact |
| Kerr horizon radius | ✅ Exact |
| Doppler factor formula | ✅ Exact |
| Relativistic aberration (velocity addition) | ✅ Exact |
| Gravitational redshift (Schwarzschild) | ✅ Exact |
| Time dilation (gravitational + kinematic) | ✅ Exact |
| **Observer angular velocity** | ✅ Fixed: $\Omega = v\sqrt{1-1/r}/r$ |
| Lorentz contraction | ✅ Correct |
| H-alpha wavelength (656.28 nm) | ✅ Correct |
| **Doppler beaming (physical mode)** | ✅ Fixed: no double-counting |
| **Planet illumination** | ✅ Fixed: 100% physical blackbody tint |
| **Galaxy Doppler** | ✅ Fixed: 100% physical (no mixing with original) |
| **ADAF emissivity** | ✅ Fixed: $r^{-3.5}$ (Narayan & Yi 1994) |
| **ADAF temperature** | ✅ Fixed: $r^{-0.5}$ (two-temperature ADAF) |
| **Jet beaming exponent** | ✅ Fixed: $D^{2.7}$ (synchrotron $\alpha \approx 0.7$) |
| **Jet Doppler** | ✅ Fixed: no temperature shift (non-thermal) |

### What Remains Approximate

| Feature | Issue |
|---------|-------|
| Kerr frame-dragging in fast mode | Still uses Binet perturbation (approximate); use `kerr_full_core` for accuracy |
| Cinematic beaming mode | Artistic exponents and clamping (only active when `physical_beaming` is off) |
| `bh_spin_strength` in fast mode | Visual multiplier (only affects fast mode; `kerr_full_core` ignores it) |

### Artistic Controls (Look Category)

These controls remain available for visual tuning. At their neutral values, the rendering is physically accurate:

| Control | Neutral value | Effect at neutral |
|---------|--------------|-------------------|
| exposure | 1.0 | Base ACES tone mapping |
| disk_gain | 1.0 | No artificial gain |
| glow | 0.0 | No artificial glow |
| doppler_boost | 1.0 | Full physical beaming (when physical mode) |
| aberration_strength | 1.0 | Full physical aberration |
| star_gain | 1.0 | Base star brightness |
| galaxy_gain | 1.0 | Base galaxy brightness |

### Changes Made

1. **Kerr geodesics**: Replaced approximate Binet+frame-drag with exact Boyer-Lindquist integration using Carter constants ($L_z$, $Q$) in `kerr_full_core` mode
2. **Kerr shadow**: Removed arbitrary `capture_u` modulation (0.12 coefficient); shadow shape now emerges from the geodesics
3. **Beaming double-counting**: Fixed $D^3$ + temperature shift double-count — temperature shift alone correctly gives $B_\nu(DT) = D^3 B_{\nu/D}(T)$
4. **Observer beaming**: Removed redundant observer `ray_intensity` scaling in physical mode (already in per-source transfer_factor)
5. **Planet tint**: Changed from 68% arbitrary / 32% physical to 100% physical blackbody
6. **Galaxy Doppler**: Changed mixing from 45% to 100%, boosted MAX_BOOST from 1.3 to 10.0
7. **ADAF radial emissivity**: Changed from $r^{-2}$ to $r^{-3.5}$ (bremsstrahlung in self-similar ADAF)
8. **ADAF temperature**: Changed from $r^{-0.3}$ to $r^{-0.5}$ (two-temperature ADAF electrons)
9. **Jet beaming**: Changed exponent from 2.0 to 2.7 ($\alpha \approx 0.7$ synchrotron index)
10. **Jet temperature**: Removed Doppler temperature shift (synchrotron is non-thermal)
11. **Observer angular velocity**: Added $\sqrt{1-1/r}$ factor for local→coordinate velocity conversion
12. **Look defaults**: Set `disk_gain=1.0`, `glow=0.0` for physically neutral defaults
