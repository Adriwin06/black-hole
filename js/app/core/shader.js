// Role: Shader configuration and compilation. Holds all compile-time parameters
//       that control the Mustache-templated GLSL (quality, modes, enabled features)
//       and exposes compile() to render the template into a ready-to-use fragment
//       shader string. Also defines the degToRad helper.

import { Mustache } from '../vendor.js';
import { diveState, hoverState } from './scenario-state.js';

export function degToRad(a) {
    return Math.PI * a / 180.0;
}

export function Shader(mustacheTemplate) {
    this.parameters = {
        n_steps: 100,
        sample_count: 1,
        resolution_scale: 1.0,
        taa_enabled: false,
        taa: {
            history_weight: 0.88,
            clip_box: 0.06,
            motion_rejection: 8.0,
            max_camera_delta: 0.08,
            motion_clip_scale: 0.6
        },
        max_revolutions: 2.0,
        rk4_integration: false,
        cinematic_tonemap: true,
        quality: 'high',
        kerr_mode: 'kerr_inspired_disk_velocity',
        accretion_disk: true,
        accretion_mode: 'thin_disk',
        disk_self_irradiation: true,
        disk_temperature: 5000.0,
        torus: {
            r0: 4.0,
            h_ratio: 0.45,
            radial_falloff: 2.5,
            opacity: 0.015,
            outer_radius: 3.5
        },
        slim: {
            h_ratio: 0.15,
            opacity: 0.6,
            puff_factor: 2.5
        },
        jet: {
            enabled: false,
            mode: 'simple',
            half_angle: 5.0,
            lorentz_factor: 3.0,
            brightness: 1.2,
            length: 30.0,
            magnetization: 10.0,
            knot_spacing: 6.0,
            corona_brightness: 1.5,
            base_width: 0.4,
            corona_extent: 0.5
        },
        grmhd: {
            enabled: true,
            r_high: 40.0,
            magnetic_beta: 10.0,
            mad_flux: 0.0,
            density_scale: 1.0,
            turbulence_amp: 1.0,
            electron_kappa: 5.0,
            magnetic_field_str: 1.0
        },
        black_hole: {
            spin_enabled: true,
            spin: 0.90,
            spin_strength: 1.0
        },
        look: {
            tonemap_mode: 1,
            exposure: 1.0,
            disk_gain: 1.0,
            glow: 0.0,
            doppler_boost: 1.0,
            aberration_strength: 1.0,
            star_gain: 0.4,
            galaxy_gain: 0.4
        },
        bloom: {
            enabled: true,
            strength: 0.35,
            threshold: 0.65,
            radius: 0.85
        },
        planet: {
            enabled: true,
            distance: 14.0,
            radius: 0.4
        },
        lorentz_contraction: true,
        gravitational_time_dilation: true,
        aberration: true,
        beaming: true,
        physical_beaming: true,
        doppler_shift: true,
        light_travel_time: true,
        time_scale: 1.0,
        turbulence_loop_enabled: false,
        turbulence_loop_seconds: 20.0,
        observer: {
            motion: false,
            distance: 11.0,
            orbital_inclination: -10
        },
        dive: {
            speed: 1.0,
            autoOrient: true
        },
        hover: {
            speed: 0.3
        },
        planetEnabled: function() {
            return this.planet.enabled && this.quality !== 'mobile';
        },
        observerMotion: function() {
            return this.observer.motion;
        }
    };

    var that = this;
    this.needsUpdate = false;

    this.hasMovingParts = function() {
        return that.parameters.accretion_disk ||
            that.parameters.jet.enabled ||
            that.parameters.planet.enabled ||
            that.parameters.observer.motion ||
            (diveState.active && !diveState.paused) ||
            (hoverState.active && !hoverState.paused);
    };

    this.compile = function() {
        if (that.parameters.kerr_mode === 'realtime_full_kerr_core') {
            that.parameters.kerr_mode = 'kerr_inspired_disk_velocity';
        }

        var kerrInspiredMode =
            (that.parameters.kerr_mode === 'kerr_inspired_disk_velocity');

        that.parameters.kerr_fast_mode = (that.parameters.kerr_mode === 'fast');
        that.parameters.kerr_inspired_mode = kerrInspiredMode;
        that.parameters.kerr_full_geodesic = false;
        that.parameters.kerr_inspired_velocity = kerrInspiredMode;

        var accMode = that.parameters.accretion_mode;
        var diskOn = that.parameters.accretion_disk;
        that.parameters.accretion_thin_disk = diskOn && (accMode === 'thin_disk');
        that.parameters.accretion_thick_torus = diskOn && (accMode === 'thick_torus');
        that.parameters.accretion_slim_disk = diskOn && (accMode === 'slim_disk');
        that.parameters.jet_enabled = that.parameters.jet.enabled;
        that.parameters.jet_simple =
            that.parameters.jet.enabled && (that.parameters.jet.mode === 'simple');
        that.parameters.jet_physical =
            that.parameters.jet.enabled && (that.parameters.jet.mode === 'physical');
        that.parameters.grmhd_enabled = !!that.parameters.grmhd.enabled;
        that.parameters.disk_self_irradiation_enabled = !!that.parameters.disk_self_irradiation;

        return Mustache.render(mustacheTemplate, that.parameters);
    };
}
