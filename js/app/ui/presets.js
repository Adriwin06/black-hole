// Role: Black hole preset library — literature-inspired parameter sets for
//       well-known astrophysical systems. These are illustrative starting
//       points rather than definitive observational fits. Each preset is
//       applied wholesale by applyBlackHolePreset() in gui.js.

/*global BH_PRESETS:true */
var BH_PRESETS = {
    'Default': {
        // Simulation defaults — use this to restore the starting configuration.
        spin_enabled: true, spin: 0.90, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'thin_disk', disk_self_irradiation: true,
        disk_temperature: 5000,
        torus: { r0: 4.0, h_ratio: 0.45, radial_falloff: 2.5, opacity: 0.015, outer_radius: 3.5 },
        slim: { h_ratio: 0.15, opacity: 0.6, puff_factor: 2.5 },
        jet: { enabled: false, mode: 'simple', half_angle: 5.0,
               lorentz_factor: 3.0, brightness: 1.2, length: 30.0,
               magnetization: 10.0, knot_spacing: 6.0, corona_brightness: 1.5,
               base_width: 0.4, corona_extent: 0.5 },
        grmhd: { enabled: true, r_high: 40.0, magnetic_beta: 10.0, mad_flux: 0.0,
                 density_scale: 1.0, turbulence_amp: 1.0, electron_kappa: 5.0,
                 magnetic_field_str: 1.0 },
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
        // Illustrative M87*-inspired preset:
        // high-spin, jet-producing, MAD-like thick torus motivated by EHT-era
        // GRMHD studies. The exact spin and plasma parameters remain model dependent.
        spin_enabled: true, spin: 0.90, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'thick_torus', disk_self_irradiation: true,
        disk_temperature: 20000,
        torus: { r0: 4.0, h_ratio: 0.45, radial_falloff: 2.5, opacity: 0.02, outer_radius: 3.5 },
        slim: { h_ratio: 0.15, opacity: 0.6, puff_factor: 2.5 },
        jet: { enabled: true, mode: 'physical', half_angle: 4.0,
               lorentz_factor: 5.0, brightness: 1.2, length: 35.0,
               magnetization: 15.0, knot_spacing: 7.0, corona_brightness: 1.5,
               base_width: 0.35, corona_extent: 0.6 },
        // EHT Paper V-style MAD-like morphology: high R_high, strong magnetic flux,
        // bright crescent, and an active jet.
        grmhd: { enabled: true, r_high: 80.0, magnetic_beta: 5.0, mad_flux: 0.8,
                 density_scale: 1.2, turbulence_amp: 1.5, electron_kappa: 4.5,
                 magnetic_field_str: 1.5 },
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
        // Illustrative Sgr A*-inspired preset:
        // moderate spin, quiescent ADAF/RIAF-like torus, no persistent jet.
        // EHT modeling does not isolate a unique spin, so this is a representative
        // visualization preset rather than a measurement claim.
        spin_enabled: true, spin: 0.50, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'thick_torus', disk_self_irradiation: true,
        disk_temperature: 15000,
        torus: { r0: 3.5, h_ratio: 0.50, radial_falloff: 2.8, opacity: 0.012, outer_radius: 3.0 },
        slim: { h_ratio: 0.15, opacity: 0.6, puff_factor: 2.5 },
        jet: { enabled: false, mode: 'physical', half_angle: 4.0,
               lorentz_factor: 3.0, brightness: 0.6, length: 20.0,
               magnetization: 10.0, knot_spacing: 5.0, corona_brightness: 1.0,
               base_width: 0.4, corona_extent: 0.5 },
        // Lower R_high and weaker MAD flux than the M87* preset, consistent with
        // a less jet-dominated and more weakly magnetized visualization.
        grmhd: { enabled: true, r_high: 20.0, magnetic_beta: 15.0, mad_flux: 0.2,
                 density_scale: 0.8, turbulence_amp: 1.8, electron_kappa: 5.0,
                 magnetic_field_str: 1.0 },
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
        // Near-extremal thin-disk X-ray-binary preset inspired by continuum-fitting
        // work on Cygnus X-1. The visible-light temperature is a rendering proxy,
        // not the true X-ray disk temperature.
        spin_enabled: true, spin: 0.99, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'thin_disk', disk_self_irradiation: true,
        disk_temperature: 12000,
        torus: { r0: 4.0, h_ratio: 0.45, radial_falloff: 2.5, opacity: 0.015, outer_radius: 3.5 },
        slim: { h_ratio: 0.15, opacity: 0.6, puff_factor: 2.5 },
        jet: { enabled: false, mode: 'simple', half_angle: 5.0,
               lorentz_factor: 2.0, brightness: 0.8, length: 20.0,
               magnetization: 10.0, knot_spacing: 6.0, corona_brightness: 1.0,
               base_width: 0.4, corona_extent: 0.5 },
        // Thin-disk-biased GRMHD-inspired morphology with moderate magnetization.
        grmhd: { enabled: true, r_high: 10.0, magnetic_beta: 20.0, mad_flux: 0.1,
                 density_scale: 1.5, turbulence_amp: 0.8, electron_kappa: 6.0,
                 magnetic_field_str: 0.8 },
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
        // Near-extremal slim-disk microquasar preset inspired by continuum-fitting
        // and superluminal-jet observations of GRS 1915+105.
        spin_enabled: true, spin: 0.98, spin_strength: 1.0,
        accretion_disk: true, accretion_mode: 'slim_disk', disk_self_irradiation: true,
        disk_temperature: 22000,
        torus: { r0: 4.0, h_ratio: 0.45, radial_falloff: 2.5, opacity: 0.015, outer_radius: 3.5 },
        slim: { h_ratio: 0.18, opacity: 0.7, puff_factor: 3.0 },
        jet: { enabled: true, mode: 'simple', half_angle: 3.0,
               lorentz_factor: 4.0, brightness: 1.0, length: 25.0,
               magnetization: 10.0, knot_spacing: 6.0, corona_brightness: 1.5,
               base_width: 0.4, corona_extent: 0.5 },
        // Strongly magnetized, super-Eddington-leaning morphology for a bright microquasar.
        grmhd: { enabled: true, r_high: 5.0, magnetic_beta: 3.0, mad_flux: 0.6,
                 density_scale: 2.0, turbulence_amp: 2.0, electron_kappa: 3.5,
                 magnetic_field_str: 2.0 },
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
        accretion_disk: true, accretion_mode: 'thin_disk', disk_self_irradiation: true,
        disk_temperature: 4500,
        torus: { r0: 4.0, h_ratio: 0.45, radial_falloff: 2.5, opacity: 0.015, outer_radius: 3.5 },
        slim: { h_ratio: 0.15, opacity: 0.6, puff_factor: 2.5 },
        jet: { enabled: false, mode: 'simple', half_angle: 5.0,
               lorentz_factor: 3.0, brightness: 1.0, length: 30.0,
               magnetization: 10.0, knot_spacing: 6.0, corona_brightness: 1.5,
               base_width: 0.4, corona_extent: 0.5 },
        grmhd: { enabled: true, r_high: 40.0, magnetic_beta: 10.0, mad_flux: 0.0,
                 density_scale: 1.0, turbulence_amp: 1.0, electron_kappa: 5.0,
                 magnetic_field_str: 1.0 },
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
        accretion_disk: true, accretion_mode: 'thin_disk', disk_self_irradiation: true,
        disk_temperature: 5000,
        torus: { r0: 4.0, h_ratio: 0.45, radial_falloff: 2.5, opacity: 0.015, outer_radius: 3.5 },
        slim: { h_ratio: 0.15, opacity: 0.6, puff_factor: 2.5 },
        jet: { enabled: false, mode: 'simple', half_angle: 5.0,
               lorentz_factor: 3.0, brightness: 1.0, length: 30.0,
               magnetization: 10.0, knot_spacing: 6.0, corona_brightness: 1.5,
               base_width: 0.4, corona_extent: 0.5 },
        // Schwarzschild: GRMHD off for clean textbook appearance
        grmhd: { enabled: false, r_high: 40.0, magnetic_beta: 10.0, mad_flux: 0.0,
                 density_scale: 1.0, turbulence_amp: 1.0, electron_kappa: 5.0,
                 magnetic_field_str: 1.0 },
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
