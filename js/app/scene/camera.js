// Role: Camera management — initialises the Three.js orbit camera at a default
//       viewing angle, synchronises it with the observer's orientation on each
//       frame, and provides a Frobenius-norm matrix distance helper used by the
//       animate loop to skip redundant renders.

function initializeCamera(camera) {

    var pitchAngle = 3.0, yawAngle = 0.0;

    // there are nicely named methods such as "lookAt" in the camera object
    // but there do not do a thing to the projection matrix due to an internal
    // representation of the camera coordinates using a quaternion (nice)
    camera.matrixWorldInverse.makeRotationX(degToRad(-pitchAngle));
    camera.matrixWorldInverse.multiply(new THREE.Matrix4().makeRotationY(degToRad(-yawAngle)));

    var m = camera.matrixWorldInverse.elements;

    camera.position.set(m[2], m[6], m[10]);
}

function updateCamera( event ) {

    // Keep matrices current when called from OrbitControls change events.
    camera.updateMatrixWorld();
    camera.matrixWorldInverse.getInverse(camera.matrixWorld);

    var zoom_dist = camera.position.length();
    var m = camera.matrixWorldInverse.elements;
    var camera_matrix;

    if (shader.parameters.observer.motion) {
        camera_matrix = new THREE.Matrix3();
    }
    else {
        camera_matrix = observer.orientation;
    }

    camera_matrix.set(
        // row-major, not the same as .elements (nice)
        // y and z swapped for a nicer coordinate system
        m[0], m[1], m[2],
        m[8], m[9], m[10],
        m[4], m[5], m[6]
    );

    if (shader.parameters.observer.motion) {

        observer.orientation = observer.orbitalFrame().multiply(camera_matrix);

    } else if (diveState && diveState.active) {
        // Dive frame: inward direction in column 1 → cam_z (look direction)
        // through the y/z swap, matching how orbitalFrame places orbital_y
        // (the look direction) in column 1.
        // Use a fresh matrix for orbit-control rotation to avoid stale refs.
        var orbitRot = new THREE.Matrix3();
        orbitRot.set(m[0], m[1], m[2], m[8], m[9], m[10], m[4], m[5], m[6]);

        var inward = diveState.direction.clone().negate();
        var up_hint = Math.abs(inward.z) < 0.99
            ? new THREE.Vector3(0, 0, 1) : new THREE.Vector3(0, 1, 0);
        var right = (new THREE.Vector3()).crossVectors(inward, up_hint).normalize();
        var up = (new THREE.Vector3()).crossVectors(right, inward).normalize();

        var diveFrame = (new THREE.Matrix4()).makeBasis(right, inward, up).linearPart();
        observer.orientation = diveFrame.multiply(orbitRot);
        shader.needsUpdate = true;
    } else {

        var p = new THREE.Vector3(
            camera_matrix.elements[6],
            camera_matrix.elements[7],
            camera_matrix.elements[8]);

        var dist = shader.parameters.observer.distance;
        observer.position.set(-p.x*dist, -p.y*dist, -p.z*dist);
        observer.velocity.set(0,0,0);
    }
}

function frobeniusDistance(matrix1, matrix2) {
    var sum = 0.0;
    for (var i in matrix1.elements) {
        var diff = matrix1.elements[i] - matrix2.elements[i];
        sum += diff*diff;
    }
    return Math.sqrt(sum);
}
