// Role: Camera management - initializes the Three.js orbit camera at a default
//       viewing angle, synchronizes it with the observer's orientation on each
//       frame, and provides a Frobenius-norm matrix distance helper used by the
//       animate loop to skip redundant renders.

import { THREE } from '../vendor.js';
import { degToRad } from '../core/shader.js';
import { diveState, hoverState } from '../core/scenarios/scenario-state.js';
import { camera, observer, shader } from '../core/runtime/runtime-state.js';

export function initializeCamera(targetCamera) {
    var pitchAngle = 3.0;
    var yawAngle = 0.0;

    targetCamera.matrixWorldInverse.makeRotationX(degToRad(-pitchAngle));
    targetCamera.matrixWorldInverse.multiply(
        new THREE.Matrix4().makeRotationY(degToRad(-yawAngle))
    );

    var m = targetCamera.matrixWorldInverse.elements;
    targetCamera.position.set(m[2], m[6], m[10]);
}

export function updateCamera() {
    camera.updateMatrixWorld();
    camera.matrixWorldInverse.getInverse(camera.matrixWorld);

    var m = camera.matrixWorldInverse.elements;
    var cameraMatrix;

    if (shader.parameters.observer.motion) {
        cameraMatrix = new THREE.Matrix3();
    } else {
        cameraMatrix = observer.orientation;
    }

    cameraMatrix.set(
        m[0], m[1], m[2],
        m[8], m[9], m[10],
        m[4], m[5], m[6]
    );

    if (shader.parameters.observer.motion) {
        observer.orientation = observer.orbitalFrame().multiply(cameraMatrix);
    } else if (diveState.active) {
        var orbitRot = new THREE.Matrix3();
        orbitRot.set(m[0], m[1], m[2], m[8], m[9], m[10], m[4], m[5], m[6]);

        var inward = diveState.direction.clone().negate();
        var upHint = Math.abs(inward.z) < 0.99
            ? new THREE.Vector3(0, 0, 1)
            : new THREE.Vector3(0, 1, 0);
        var right = (new THREE.Vector3()).crossVectors(inward, upHint).normalize();
        var up = (new THREE.Vector3()).crossVectors(right, inward).normalize();

        var diveFrame = (new THREE.Matrix4()).makeBasis(right, inward, up).linearPart();
        observer.orientation = diveFrame.multiply(orbitRot);
        shader.needsUpdate = true;
    } else if (hoverState.active) {
        var hoverOrbitRot = new THREE.Matrix3();
        hoverOrbitRot.set(m[0], m[1], m[2], m[8], m[9], m[10], m[4], m[5], m[6]);

        var hoverInward = hoverState.direction.clone().negate();
        var hoverUpHint = Math.abs(hoverInward.z) < 0.99
            ? new THREE.Vector3(0, 0, 1)
            : new THREE.Vector3(0, 1, 0);
        var hoverRight = (new THREE.Vector3())
            .crossVectors(hoverInward, hoverUpHint)
            .normalize();
        var hoverUp = (new THREE.Vector3())
            .crossVectors(hoverRight, hoverInward)
            .normalize();

        var hoverFrame = (new THREE.Matrix4())
            .makeBasis(hoverRight, hoverInward, hoverUp)
            .linearPart();
        observer.orientation = hoverFrame.multiply(hoverOrbitRot);
        shader.needsUpdate = true;
    } else {
        var p = new THREE.Vector3(
            cameraMatrix.elements[6],
            cameraMatrix.elements[7],
            cameraMatrix.elements[8]
        );

        var dist = shader.parameters.observer.distance;
        observer.position.set(-p.x * dist, -p.y * dist, -p.z * dist);
        observer.velocity.set(0, 0, 0);
    }
}

export function frobeniusDistance(matrix1, matrix2) {
    var sum = 0.0;
    for (var i in matrix1.elements) {
        var diff = matrix1.elements[i] - matrix2.elements[i];
        sum += diff * diff;
    }
    return Math.sqrt(sum);
}
