export function setupAccretionLookFolders(options) {
    options = options || {};

    var gui = options.gui;
    var p = options.parameters;
    var observer = options.observer;
    var addControl = options.addControl;
    var setControlVisible = options.setControlVisible;
    var setControlsVisible = options.setControlsVisible;
    var getUpdateDependentVisibility = options.getUpdateDependentVisibility || function() { return null; };
    var updateShader = options.updateShader || function() {};
    var updateUniformsLive = options.updateUniformsLive || function() {};
    var markShaderDirty = options.markShaderDirty || function() {};
    var diskTemperatureMin = options.diskTemperatureMin;
    var diskTemperatureMax = options.diskTemperatureMax;

    function syncDependentVisibility() {
        var updateDependentVisibility = getUpdateDependentVisibility();
        if (typeof updateDependentVisibility === 'function') {
            updateDependentVisibility();
        }
    }

    var diskFolder = gui.addFolder('Accretion disk');
    addControl(diskFolder, p, 'accretion_disk', {
        name: 'enabled',
        trackBlackHolePreset: true,
        onChange: function() {
            syncDependentVisibility();
            updateShader();
        },
        help: 'Toggles thermal emission from the accretion disk.'
    });
    var accretionModeCtrl = addControl(diskFolder, p, 'accretion_mode', {
        name: 'accretion type',
        options: ['thin_disk', 'thick_torus', 'slim_disk'],
        trackBlackHolePreset: true,
        onChange: function() {
            syncDependentVisibility();
            updateShader();
            observer.turbulenceTimeOffset = -observer.time;
            markShaderDirty();
        },
        help: 'Thin disk: zero-torque Shakura-Sunyaev / Novikov-Thorne-style proxy (quasars/XRBs). Thick torus: ADAF/RIAF (M87*/Sgr A*). Slim disk: super-Eddington.'
    });
    var diskSelfIrradiationCtrl = addControl(diskFolder, p, 'disk_self_irradiation', {
        name: 'self-irradiation',
        trackBlackHolePreset: true,
        onChange: function() { updateShader(); },
        help: 'Heuristic Cunningham-inspired inner-disk brightening. Increases local flux near the ISCO with a spin-scaled falloff; this is not recursive ray-traced returning radiation.'
    });
    var diskTempCtrl = addControl(diskFolder, p, 'disk_temperature', {
        min: diskTemperatureMin,
        max: diskTemperatureMax,
        step: 1,
        name: 'temperature (K)',
        trackBlackHolePreset: true,
        onChange: updateUniformsLive,
        help: 'Rest-frame disk color temperature before relativistic shifts.'
    });

    var torusCenterCtrl = addControl(diskFolder, p.torus, 'r0', { min: 1.5, max: 10.0, step: 0.1, name: 'torus center r', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Center radius of the torus in r_s units (thick torus mode only).' });
    var torusHRCtrl = addControl(diskFolder, p.torus, 'h_ratio', { min: 0.1, max: 1.0, step: 0.01, name: 'torus H/R', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Height-to-radius ratio of the torus cross-section (thick torus mode only).' });
    var torusFalloffCtrl = addControl(diskFolder, p.torus, 'radial_falloff', { min: 0.5, max: 5.0, step: 0.1, name: 'radial falloff', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Power-law index for radial emissivity decay. Physical ADAF: 2-4. Lower = flatter profile, higher = more concentrated.' });
    var torusOpacityCtrl = addControl(diskFolder, p.torus, 'opacity', { min: 0.001, max: 0.15, step: 0.001, name: 'opacity', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Absorption coefficient. Low = optically thin (ADAF). Higher = more self-shielding and defined torus surface.' });
    var torusOuterCtrl = addControl(diskFolder, p.torus, 'outer_radius', { min: 1.5, max: 8.0, step: 0.1, name: 'outer extent', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Outer edge multiplier relative to torus center r0. Controls how far the torus extends.' });

    var slimHRCtrl = addControl(diskFolder, p.slim, 'h_ratio', { min: 0.05, max: 0.5, step: 0.01, name: 'slim H/R', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Base height-to-radius ratio of the slim disk. Higher = geometrically thicker.' });
    var slimOpacityCtrl = addControl(diskFolder, p.slim, 'opacity', { min: 0.1, max: 2.0, step: 0.01, name: 'slim opacity', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Absorption coefficient. Higher = more opaque, surface-like. Lower = more translucent, volumetric.' });
    var slimPuffCtrl = addControl(diskFolder, p.slim, 'puff_factor', { min: 0.0, max: 6.0, step: 0.1, name: 'ISCO puff', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'How much the disk puffs up near the ISCO due to radiation pressure. Higher = more pronounced thickening.' });

    var torusRows = [torusCenterCtrl, torusHRCtrl, torusFalloffCtrl, torusOpacityCtrl, torusOuterCtrl];
    var slimRows = [slimHRCtrl, slimOpacityCtrl, slimPuffCtrl];

    var grmhdEnabledCtrl = addControl(diskFolder, p.grmhd, 'enabled', {
        name: 'GRMHD-inspired',
        trackBlackHolePreset: true,
        onChange: function() {
            syncDependentVisibility();
            updateShader();
        },
        help: 'Enables GRMHD-inspired morphology controls: two-temperature plasma weighting, synchrotron-inspired emissivity corrections, MRI-inspired turbulence, and magnetic-field effects. These are semi-analytic visual models, not full GRMHD evolution.'
    });
    var grmhdRHighCtrl = addControl(diskFolder, p.grmhd, 'r_high', { min: 1.0, max: 160.0, step: 1.0, name: 'R_high', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Renderer-side Ti:Te / electron-heating proxy. Higher R_high suppresses high-beta disk-body emission most strongly in the GRMHD torus branch, while thin/slim disks keep the effect deliberately milder so the RGB disk remains visible.' });
    var grmhdBetaCtrl = addControl(diskFolder, p.grmhd, 'magnetic_beta', { min: 0.01, max: 100.0, step: 0.1, name: 'plasma beta', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Ratio of gas-to-magnetic pressure in the disk midplane. Low beta = strongly magnetized, high beta = gas dominated. GRMHD simulations: beta ~ 1-30.' });
    var grmhdMADCtrl = addControl(diskFolder, p.grmhd, 'mad_flux', { min: 0.0, max: 1.0, step: 0.01, name: 'MAD flux', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Magnetically Arrested Disk saturation. 0 = SANE (standard), 1 = full MAD. MAD state creates stronger B-fields, prominent m=1 spiral arms, and more powerful jets (Tchekhovskoy+ 2011).' });
    var grmhdDensityCtrl = addControl(diskFolder, p.grmhd, 'density_scale', { min: 0.1, max: 5.0, step: 0.1, name: 'density scale', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Normalization of the GRMHD density profile. Higher values increase overall emissivity and absorption.' });
    var grmhdTurbCtrl = addControl(diskFolder, p.grmhd, 'turbulence_amp', { min: 0.0, max: 3.0, step: 0.1, name: 'MRI turbulence', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Amplitude of MRI-driven turbulent density fluctuations. 0 = smooth, 1 = typical GRMHD, 3 = strongly turbulent. Controls log-normal density PDF and spiral arm contrast.' });
    var grmhdKappaCtrl = addControl(diskFolder, p.grmhd, 'electron_kappa', { min: 2.5, max: 8.0, step: 0.1, name: 'kappa (non-thermal)', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Kappa-distribution index for non-thermal electrons from magnetic reconnection. Kappa ~ 3.5 = strong non-thermal tail, Kappa ~ 8 = nearly thermal.' });
    var grmhdBFieldCtrl = addControl(diskFolder, p.grmhd, 'magnetic_field_str', { min: 0.1, max: 5.0, step: 0.1, name: 'B-field strength', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Overall magnetic field strength scaling. Affects synchrotron emissivity, self-absorption, and electron temperature.' });
    var grmhdRows = [grmhdRHighCtrl, grmhdBetaCtrl, grmhdMADCtrl, grmhdDensityCtrl, grmhdTurbCtrl, grmhdKappaCtrl, grmhdBFieldCtrl];

    setControlVisible(accretionModeCtrl, !!p.accretion_disk);
    setControlVisible(diskSelfIrradiationCtrl, !!p.accretion_disk);
    setControlVisible(diskTempCtrl, !!p.accretion_disk);
    setControlsVisible(torusRows, !!p.accretion_disk && p.accretion_mode === 'thick_torus');
    setControlsVisible(slimRows, !!p.accretion_disk && p.accretion_mode === 'slim_disk');
    setControlVisible(grmhdEnabledCtrl, !!p.accretion_disk);
    setControlsVisible(grmhdRows, !!p.accretion_disk && !!p.grmhd.enabled);
    diskFolder.open();

    var jetFolder = gui.addFolder('Relativistic jets');
    addControl(jetFolder, p.jet, 'enabled', { name: 'enabled', trackBlackHolePreset: true, onChange: function() { syncDependentVisibility(); updateShader(); }, help: 'Bipolar analytic jets aligned with the spin axis. The detailed mode adds GRMHD-inspired structure and Blandford-Znajek-like power scaling.' });
    var jetModeCtrl = addControl(jetFolder, p.jet, 'mode', { name: 'jet model', options: ['simple', 'physical'], trackBlackHolePreset: true, onChange: function() { syncDependentVisibility(); updateShader(); }, help: 'Simple: smooth analytic jet. Physical: more detailed GRMHD-inspired jet model with spine/sheath, reconfinement knots, corona base, and disk occultation.' });
    var jetAngleCtrl = addControl(jetFolder, p.jet, 'half_angle', { min: 1.0, max: 25.0, step: 0.5, name: 'half-angle (deg)', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Opening half-angle of the jet. Typical AGN jets: 2-7 deg. Parabolic collimation narrows the beam with distance.' });
    var jetLorentzCtrl = addControl(jetFolder, p.jet, 'lorentz_factor', { min: 1.1, max: 20.0, step: 0.1, name: 'Lorentz Gamma', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Bulk Lorentz factor for relativistic beaming.' });
    var jetBrightCtrl = addControl(jetFolder, p.jet, 'brightness', { min: 0.05, max: 3.0, step: 0.01, name: 'brightness', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Overall jet synchrotron emission strength.' });
    var jetLengthCtrl = addControl(jetFolder, p.jet, 'length', { min: 5.0, max: 60.0, step: 1.0, name: 'length (r_s)', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Visible jet length in Schwarzschild radii.' });
    var jetMagCtrl = addControl(jetFolder, p.jet, 'magnetization', { min: 1.0, max: 50.0, step: 0.5, name: 'sigma (magnetization)', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Plasma magnetization at jet base.' });
    var jetKnotCtrl = addControl(jetFolder, p.jet, 'knot_spacing', { min: 2.0, max: 15.0, step: 0.5, name: 'knot spacing', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Spacing of reconfinement shock knots.' });
    var jetCoronaCtrl = addControl(jetFolder, p.jet, 'corona_brightness', { min: 0.0, max: 5.0, step: 0.1, name: 'corona glow', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Brightness of the jet-corona connection at the base.' });
    var jetBaseWidthCtrl = addControl(jetFolder, p.jet, 'base_width', { min: 0.0, max: 1.0, step: 0.01, name: 'base width', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Controls how wide the jet funnel is at the base.' });
    var jetCoronaExtentCtrl = addControl(jetFolder, p.jet, 'corona_extent', { min: 0.0, max: 1.0, step: 0.01, name: 'corona extent', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Radial spread of the corona base emission.' });
    var jetCommonCtrls = [jetModeCtrl, jetAngleCtrl, jetLorentzCtrl, jetBrightCtrl, jetLengthCtrl];
    var jetPhysicalCtrls = [jetMagCtrl, jetKnotCtrl, jetCoronaCtrl, jetBaseWidthCtrl, jetCoronaExtentCtrl];

    var spinFolder = gui.addFolder('Black hole');
    addControl(spinFolder, p.black_hole, 'spin_enabled', { name: 'rotation enabled', trackBlackHolePreset: true, onChange: function() { syncDependentVisibility(); updateUniformsLive(); }, help: 'Enables Kerr-like rotation effects. Disable for a Schwarzschild-style shadow.' });
    var spinCtrl = addControl(spinFolder, p.black_hole, 'spin', { min: -0.99, max: 0.99, step: 0.01, name: 'a/M', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Signed Kerr spin.' });
    var spinStrengthCtrl = addControl(spinFolder, p.black_hole, 'spin_strength', { min: 0.0, max: 1.4, step: 0.01, name: 'shadow squeeze', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Visual strength multiplier for rotation-induced asymmetry in the shadow.' });
    spinFolder.open();

    var lookFolder = gui.addFolder('Look');
    addControl(lookFolder, p.look, 'tonemap_mode', { options: { 'ACES Filmic': 0, 'AgX': 1, 'Scientific (Log)': 2 }, name: 'tonemapper', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'ACES: classic cinematic. AgX: better saturation handling for extreme HDR. Scientific: logarithmic false-color like EHT papers.' });
    addControl(lookFolder, p.look, 'exposure', { min: 0.6, max: 2.5, step: 0.01, name: 'exposure', onChange: updateUniformsLive, help: 'Overall tone-mapped brightness.' });
    var lookDiskGainCtrl = addControl(lookFolder, p.look, 'disk_gain', { min: 0.4, max: 4.0, step: 0.01, name: 'disk intensity', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Brightness multiplier applied to the accretion disk emission.' });
    var lookGlowCtrl = addControl(lookFolder, p.look, 'glow', { min: 0.0, max: 2.0, step: 0.01, name: 'inner glow', trackBlackHolePreset: true, onChange: updateUniformsLive, help: 'Extra bloom-like emphasis near the hotter inner disk region.' });
    var lookDopplerBoostCtrl = addControl(lookFolder, p.look, 'doppler_boost', { min: 0.0, max: 2.5, step: 0.01, name: 'doppler boost', onChange: updateUniformsLive, help: 'Controls visual strength of relativistic beaming contrast.' });
    var lookAberrationStrengthCtrl = addControl(lookFolder, p.look, 'aberration_strength', { min: 0.0, max: 3.0, step: 0.01, name: 'aberration strength', onChange: updateUniformsLive, help: 'Scales apparent directional warping from observer motion.' });
    addControl(lookFolder, p.look, 'star_gain', { min: 0.0, max: 2.5, step: 0.01, name: 'star gain', onChange: updateUniformsLive, help: 'Brightness multiplier for background stars.' });
    addControl(lookFolder, p.look, 'galaxy_gain', { min: 0.0, max: 2.5, step: 0.01, name: 'galaxy gain', onChange: updateUniformsLive, help: 'Brightness multiplier for the background galaxy map.' });
    lookFolder.open();

    var ppFolder = gui.addFolder('Post-processing');
    addControl(ppFolder, p.bloom, 'enabled', { name: 'bloom', trackBlackHolePreset: true, onChange: function() { syncDependentVisibility(); markShaderDirty(); }, help: 'Screen-space Gaussian bloom approximating lens glow and optical spill.' });
    var bloomStrengthCtrl = addControl(ppFolder, p.bloom, 'strength', { min: 0.0, max: 2.0, step: 0.01, name: 'bloom strength', trackBlackHolePreset: true, onChange: markShaderDirty, help: 'Overall intensity of the bloom glow added to the image.' });
    var bloomThresholdCtrl = addControl(ppFolder, p.bloom, 'threshold', { min: 0.0, max: 1.0, step: 0.01, name: 'bloom threshold', trackBlackHolePreset: true, onChange: markShaderDirty, help: 'Minimum pixel brightness that contributes to bloom.' });
    var bloomRadiusCtrl = addControl(ppFolder, p.bloom, 'radius', { min: 0.0, max: 1.0, step: 0.01, name: 'bloom radius', trackBlackHolePreset: true, onChange: markShaderDirty, help: 'Controls the width of the bloom halo.' });
    ppFolder.open();

    return {
        accretionModeCtrl: accretionModeCtrl,
        diskSelfIrradiationCtrl: diskSelfIrradiationCtrl,
        diskTempCtrl: diskTempCtrl,
        torusRows: torusRows,
        slimRows: slimRows,
        grmhdEnabledCtrl: grmhdEnabledCtrl,
        grmhdRows: grmhdRows,
        jetCommonCtrls: jetCommonCtrls,
        jetPhysicalCtrls: jetPhysicalCtrls,
        spinCtrl: spinCtrl,
        spinStrengthCtrl: spinStrengthCtrl,
        lookDiskGainCtrl: lookDiskGainCtrl,
        lookGlowCtrl: lookGlowCtrl,
        lookDopplerBoostCtrl: lookDopplerBoostCtrl,
        lookAberrationStrengthCtrl: lookAberrationStrengthCtrl,
        bloomStrengthCtrl: bloomStrengthCtrl,
        bloomThresholdCtrl: bloomThresholdCtrl,
        bloomRadiusCtrl: bloomRadiusCtrl
    };
}
