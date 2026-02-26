// Role: Black hole preset library — physically motivated parameter sets for
//       well-known astrophysical objects. Each preset is applied wholesale by
//       applyBlackHolePreset() in gui.js to override all relevant shader and
//       observer parameters in one click.

/*global BH_PRESETS:true */
var BH_PRESETS = {
    'Default': {
        // Simulation defaults — use this to restore the starting configuration.
        spin_enabled: true, spin: 0.90, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'thin_disk',
        disk_temperature: 5000,
        torus: { r0: 4.0, h_ratio: 0.45 },
        jet: { enabled: false, mode: 'simple', half_angle: 5.0,
               lorentz_factor: 3.0, brightness: 1.2, length: 30.0,
               magnetization: 10.0, knot_spacing: 6.0, corona_brightness: 1.5 },
        observer: { distance: 11.0 },
        beaming: true,
        physical_beaming: true,
        doppler_shift: true,
        disk_gain: 1.0,
        glow: 0.0,
        tonemap_mode: 1,
        bloom: { enabled: true, strength: 0.35, threshold: 0.65, radius: 0.85 }
    },
    'M87*': {
        // Supermassive BH in M87 (Virgo A), M = 6.5×10⁹ M☉, first EHT image (2019).
        // Spin: a/M ≈ 0.90 ± 0.1 (Tamburini et al. 2019, twisted photon OAM).
        // Confirmed by Daly 2019 (outflow method): a/M = 1.00 ± 0.15.
        // Accretion: ADAF/RIAF thick torus — low-luminosity AGN, sub-mm synchrotron.
        // Prominent relativistic jet (HST-1, superluminal knots).
        spin_enabled: true, spin: 0.90, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'thick_torus',
        disk_temperature: 20000,
        torus: { r0: 4.0, h_ratio: 0.45 },
        jet: { enabled: true, mode: 'physical', half_angle: 4.0,
               lorentz_factor: 5.0, brightness: 1.2, length: 35.0,
               magnetization: 15.0, knot_spacing: 7.0, corona_brightness: 2.0 },
        observer: { distance: 11.0 },
        beaming: true,
        physical_beaming: true,
        doppler_shift: true,
        disk_gain: 1.0,
        glow: 0.0,
        tonemap_mode: 1,
        bloom: { enabled: true, strength: 0.40, threshold: 0.55, radius: 0.90 }
    },
    'Sgr A*': {
        // Milky Way centre SMBH, M = 4.297 ± 0.012 × 10⁶ M☉, EHT image (2022).
        // Spin: highly debated — a* < 0.1 (Fragione & Loeb 2020) to
        //   a* = 0.90 ± 0.06 (Daly et al. 2023). Mid-range a* ≈ 0.50 adopted.
        // Accretion: RIAF/ADAF; quiescent low-luminosity state.
        // No persistent jet detected observationally.
        spin_enabled: true, spin: 0.50, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'thick_torus',
        disk_temperature: 15000,
        torus: { r0: 3.5, h_ratio: 0.50 },
        jet: { enabled: false, mode: 'physical', half_angle: 4.0,
               lorentz_factor: 3.0, brightness: 0.6, length: 20.0,
               magnetization: 10.0, knot_spacing: 5.0, corona_brightness: 1.0 },
        observer: { distance: 11.0 },
        beaming: true,
        physical_beaming: true,
        doppler_shift: true,
        disk_gain: 1.0,
        glow: 0.0,
        tonemap_mode: 1,
        bloom: { enabled: true, strength: 0.35, threshold: 0.60, radius: 0.85 }
    },
    'Cygnus X-1': {
        // Stellar-mass BH X-ray binary with HDE 226868 (blue supergiant).
        // M ≈ 21.2 M☉ (Miller-Jones et al. 2021).
        // Spin: a/M > 0.983 at 3σ (Gou et al. 2011, continuum fitting) — capped 0.99.
        // Accretion: geometrically thin Shakura-Sunyaev disk (soft state).
        // Inner disk temperature ~10⁷ K (X-ray); visual proxy shown here.
        // Microquasar with transient jets; disabled here for canonical soft state.
        spin_enabled: true, spin: 0.99, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'thin_disk',
        disk_temperature: 12000,
        torus: { r0: 4.0, h_ratio: 0.45 },
        jet: { enabled: false, mode: 'simple', half_angle: 5.0,
               lorentz_factor: 2.0, brightness: 0.8, length: 20.0,
               magnetization: 10.0, knot_spacing: 6.0, corona_brightness: 1.0 },
        observer: { distance: 8.0 },
        beaming: true,
        physical_beaming: true,
        doppler_shift: true,
        disk_gain: 1.0,
        glow: 0.0,
        tonemap_mode: 1,
        bloom: { enabled: true, strength: 0.30, threshold: 0.60, radius: 0.80 }
    },
    'GRS 1915+105': {
        // Stellar-mass BH microquasar, M = 12.4 ± 2 M☉.
        // First galactic object with superluminal jets (v_app ~ 1.25c).
        // Spin: a/M > 0.98, near-extreme Kerr (McClintock et al. 2006).
        // Rotates ≥ 950 rev/s. Super-Eddington accretion episodes → slim disk.
        // Jets at ~0.9c; Lorentz factor Γ ≈ 2–5 (apparent superluminal).
        spin_enabled: true, spin: 0.98, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'slim_disk',
        disk_temperature: 22000,
        torus: { r0: 4.0, h_ratio: 0.45 },
        jet: { enabled: true, mode: 'simple', half_angle: 3.0,
               lorentz_factor: 4.0, brightness: 1.0, length: 25.0,
               magnetization: 10.0, knot_spacing: 6.0, corona_brightness: 1.5 },
        observer: { distance: 9.0 },
        beaming: true,
        physical_beaming: true,
        doppler_shift: true,
        disk_gain: 1.0,
        glow: 0.0,
        tonemap_mode: 1,
        bloom: { enabled: true, strength: 0.35, threshold: 0.55, radius: 0.85 }
    },
    'Gargantua (Interstellar visuals)': {
        spin_enabled: true, spin: 0.7, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'thin_disk',
        disk_temperature: 5800,
        torus: { r0: 4.0, h_ratio: 0.45 },
        jet: { enabled: false, mode: 'simple', half_angle: 5.0,
               lorentz_factor: 3.0, brightness: 1.0, length: 30.0,
               magnetization: 10.0, knot_spacing: 6.0, corona_brightness: 1.5 },
        observer: { distance: 11.0 },
        beaming: false,
        physical_beaming: false,
        doppler_shift: false,
        disk_gain: 2.0,
        glow: 1.0,
        tonemap_mode: 0,
        bloom: { enabled: true, strength: 0.65, threshold: 0.45, radius: 0.92 }
    },
    'Schwarzschild': {
        // Idealised non-rotating black hole (a/M = 0).
        // Classical textbook case: symmetric circular shadow, no frame dragging.
        spin_enabled: false, spin: 0.0, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'thin_disk',
        disk_temperature: 5000,
        torus: { r0: 4.0, h_ratio: 0.45 },
        jet: { enabled: false, mode: 'simple', half_angle: 5.0,
               lorentz_factor: 3.0, brightness: 1.0, length: 30.0,
               magnetization: 10.0, knot_spacing: 6.0, corona_brightness: 1.5 },
        observer: { distance: 11.0 },
        beaming: true,
        physical_beaming: true,
        doppler_shift: true,
        disk_gain: 1.0,
        glow: 0.0,
        tonemap_mode: 1,
        bloom: { enabled: true, strength: 0.30, threshold: 0.65, radius: 0.80 }
    }
};
