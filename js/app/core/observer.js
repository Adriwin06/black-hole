// Role: Observer entity - tracks position, velocity, and orientation of the
//       in-simulation camera/observer. Handles circular orbital motion and the
//       simulation-time scaling used for orbit/hover modes. Also exports
//       formatThousands.

import { THREE } from '../vendor.js';
import { degToRad } from './shader.js';

var boundShader = null;

export function bindObserverShader(shaderInstance) {
    boundShader = shaderInstance || null;
}

export function formatThousands(value) {
    return Math.round(value).toString().replace(/\B(?=(\d{3})+(?!\d))/g, " ");
}

// Close stationary views remain available, but circular observer motion is
// clamped to the stable Schwarzschild-orbit regime for massive observers.
export var OBSERVER_DISTANCE_MIN = 1.5;
export var OBSERVER_ORBIT_MIN = 3.0;
export var OBSERVER_DISTANCE_MAX = 30.0;
export var PLANET_ORBIT_MIN = 3.0;

export function clampObserverDistance(distance, motionEnabled) {
    var minDistance = motionEnabled ? OBSERVER_ORBIT_MIN : OBSERVER_DISTANCE_MIN;
    return Math.max(minDistance, Math.min(OBSERVER_DISTANCE_MAX, distance));
}

export function clampPlanetOrbitDistance(distance) {
    return Math.max(PLANET_ORBIT_MIN, distance);
}

export function Observer() {
    this.position = new THREE.Vector3(10, 0, 0);
    this.velocity = new THREE.Vector3(0, 1, 0);
    this.orientation = new THREE.Matrix3();
    this.time = 0.0;
}

Observer.prototype.orbitalFrame = function() {
    var orbitalY = (new THREE.Vector3())
        .subVectors(this.velocity.clone().normalize().multiplyScalar(4.0), this.position)
        .normalize();

    var orbitalZ = (new THREE.Vector3())
        .crossVectors(this.position, orbitalY).normalize();
    var orbitalX = (new THREE.Vector3()).crossVectors(orbitalY, orbitalZ);

    return (new THREE.Matrix4()).makeBasis(
        orbitalX,
        orbitalY,
        orbitalZ
    ).linearPart();
};

Observer.prototype.move = function(dt) {
    if (!boundShader) return;
    dt *= boundShader.parameters.time_scale;

    var r;
    var v = 0;

    if (boundShader.parameters.observer.motion) {
        r = clampObserverDistance(boundShader.parameters.observer.distance, true);
        boundShader.parameters.observer.distance = r;
        v = 1.0 / Math.sqrt(2.0 * (r - 1.0));

        var angVel = v * Math.sqrt(1.0 - 1.0 / r) / r;
        var angle = this.time * angVel;
        var s = Math.sin(angle);
        var c = Math.cos(angle);

        this.position.set(c * r, s * r, 0);
        this.velocity.set(-s * v, c * v, 0);

        var alpha = degToRad(boundShader.parameters.observer.orbital_inclination);
        var orbitCoords = (new THREE.Matrix4()).makeRotationY(alpha);

        this.position.applyMatrix4(orbitCoords);
        this.velocity.applyMatrix4(orbitCoords);
    } else {
        r = this.position.length();
    }

    if (boundShader.parameters.gravitational_time_dilation) {
        if (v > 0) {
            dt = dt / Math.sqrt(Math.max(1.0 - 1.5 / r, 0.001));
        } else {
            dt = dt / Math.sqrt(Math.max(1.0 - 1.0 / r, 0.001));
        }
    }

    this.time += dt;
};


