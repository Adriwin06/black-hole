import {
    camera,
    observer,
    shader,
    cameraControls,
    cameraPan,
    distanceController,
    refreshAllControllersGlobal
} from '../../core/runtime/runtime-state.js';
import { diveState, hoverState } from '../../core/scenarios/scenario-state.js';
import { updateCamera } from '../../scene/camera.js';
import { seekDive } from '../../core/scenarios/dive.js';
import { seekHover } from '../../core/scenarios/hover.js';
import { getBlackHoleUiBinding } from '../../core/runtime/runtime-registry.js';

function clonePresentationPathValue(value) {
    return JSON.parse(JSON.stringify(value));
}

export function resolvePresentationPath(path) {
    if (!shader || !path || typeof path !== 'string') return null;

    var clean = path.trim();
    var root = null;
    var parts = [];

    if (clean.indexOf('cameraPan.') === 0) {
        root = cameraPan;
        parts = clean.substring('cameraPan.'.length).split('.');
    } else if (clean.indexOf('camera.') === 0) {
        root = camera;
        parts = clean.substring('camera.'.length).split('.');
    } else if (clean.indexOf('observerState.') === 0) {
        root = observer;
        parts = clean.substring('observerState.'.length).split('.');
    } else if (clean.indexOf('dive.') === 0) {
        root = diveState;
        parts = clean.substring('dive.'.length).split('.');
    } else if (clean.indexOf('hover.') === 0) {
        root = hoverState;
        parts = clean.substring('hover.'.length).split('.');
    } else if (clean.indexOf('params.') === 0) {
        root = shader.parameters;
        parts = clean.substring('params.'.length).split('.');
    } else if (clean.indexOf('shader.parameters.') === 0) {
        root = shader.parameters;
        parts = clean.substring('shader.parameters.'.length).split('.');
    } else {
        root = shader.parameters;
        parts = clean.split('.');
    }

    if (!parts.length) return null;
    return { root: root, parts: parts, originalPath: clean };
}

export function presentationPathNeedsCompile(path) {
    if (typeof path !== 'string') return false;
    var clean = path.trim();
    if (!clean) return false;

    if (clean.indexOf('params.') === 0) clean = clean.substring('params.'.length);
    if (clean.indexOf('shader.parameters.') === 0) {
        clean = clean.substring('shader.parameters.'.length);
    }

    switch (clean) {
        case 'kerr_mode':
        case 'accretion_disk':
        case 'accretion_mode':
        case 'disk_self_irradiation':
        case 'jet.enabled':
        case 'jet.mode':
        case 'grmhd.enabled':
        case 'planet.enabled':
        case 'aberration':
        case 'beaming':
        case 'physical_beaming':
        case 'doppler_shift':
        case 'light_travel_time':
        case 'gravitational_time_dilation':
        case 'lorentz_contraction':
        case 'cinematic_tonemap':
        case 'observer.motion':
        case 'quality':
        case 'taa_enabled':
        case 'rk4_integration':
            return true;
        default:
            return false;
    }
}

export function refreshPresentationUiBindings() {
    var refreshControllers = (typeof getBlackHoleUiBinding === 'function')
        ? getBlackHoleUiBinding('refreshControllers')
        : refreshAllControllersGlobal;
    if (typeof refreshControllers === 'function') {
        refreshControllers();
    }

    var observerDistanceBinding = (typeof getBlackHoleUiBinding === 'function')
        ? getBlackHoleUiBinding('observerDistance')
        : (typeof distanceController !== 'undefined' ? distanceController : null);
    if (observerDistanceBinding &&
        typeof observerDistanceBinding.updateDisplay === 'function') {
        observerDistanceBinding.updateDisplay();
    }
}

export function setPresentationInteractionLock(locked) {
    if (typeof cameraControls !== 'undefined' && cameraControls) {
        cameraControls.enabled = !locked;
    }
}

export function setPresentationPathValue(path, value) {
    var cleanPath = (typeof path === 'string') ? path.trim() : '';
    if (cleanPath === 'dive.currentR') {
        var diveRadius = parseFloat(value);
        if (!isFinite(diveRadius) || (!diveState.active && !diveState.reachedSingularity)) {
            return false;
        }
        var diveChanged = Math.abs(diveState.currentR - diveRadius) > 1e-8 ||
            !diveState.timelineDriven || !diveState.paused;
        seekDive(diveRadius, { timelineDriven: true });
        return diveChanged;
    }
    if (cleanPath === 'hover.currentR') {
        var hoverRadius = parseFloat(value);
        if (!isFinite(hoverRadius) || !hoverState.active) {
            return false;
        }
        var hoverChanged = Math.abs(hoverState.currentR - hoverRadius) > 1e-8 ||
            !hoverState.timelineDriven || !hoverState.paused;
        seekHover(hoverRadius, { timelineDriven: true });
        return hoverChanged;
    }

    var resolved = resolvePresentationPath(path);
    if (!resolved) return false;

    var obj = resolved.root;
    var parts = resolved.parts;
    for (var i = 0; i < parts.length - 1; i++) {
        var key = parts[i];
        if (!obj || typeof obj !== 'object' || !(key in obj)) return false;
        obj = obj[key];
    }
    if (!obj || typeof obj !== 'object') return false;

    var leaf = parts[parts.length - 1];
    if (!(leaf in obj)) return false;

    var current = obj[leaf];
    var next = value;
    if (typeof current === 'number') {
        var numeric = parseFloat(next);
        if (!isFinite(numeric)) return false;
        next = numeric;
    } else if (typeof current === 'boolean') {
        next = !!next;
    }

    var changed;
    if (typeof current === 'number' && typeof next === 'number') {
        changed = Math.abs(current - next) > 1e-8;
    } else {
        changed = current !== next;
    }
    if (!changed) return false;

    obj[leaf] = next;

    var isCameraOrObserverPath =
        resolved.originalPath.indexOf('cameraPan.') === 0 ||
        resolved.originalPath.indexOf('camera.') === 0 ||
        resolved.originalPath.indexOf('observerState.') === 0 ||
        resolved.originalPath.indexOf('observer.') === 0 ||
        resolved.originalPath.indexOf('params.observer.') === 0 ||
        resolved.originalPath.indexOf('shader.parameters.observer.') === 0;

    if (resolved.originalPath.indexOf('camera.quaternion.') === 0 &&
        camera && camera.quaternion &&
        typeof camera.quaternion.normalize === 'function') {
        camera.quaternion.normalize();
    }
    if (resolved.originalPath.indexOf('camera.') === 0 &&
        camera && typeof camera.updateMatrixWorld === 'function') {
        camera.updateMatrixWorld(true);
    }

    if (isCameraOrObserverPath && camera && typeof updateCamera === 'function') {
        updateCamera();
    }

    shader.needsUpdate = true;
    return true;
}

export function getPresentationPathValue(path) {
    var resolved = resolvePresentationPath(path);
    if (!resolved) return undefined;

    var obj = resolved.root;
    var parts = resolved.parts;
    for (var i = 0; i < parts.length - 1; i++) {
        var key = parts[i];
        if (!obj || typeof obj !== 'object' || !(key in obj)) return undefined;
        obj = obj[key];
    }
    if (!obj || typeof obj !== 'object') return undefined;

    var leaf = parts[parts.length - 1];
    if (!(leaf in obj)) return undefined;

    var value = obj[leaf];
    if (value && typeof value === 'object') {
        return clonePresentationPathValue(value);
    }
    return value;
}


