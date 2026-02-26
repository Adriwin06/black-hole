// Role: Observer entity — tracks position, velocity, and orientation of the
//       in-simulation camera/observer. Handles circular orbital motion with
//       full special-relativistic time dilation. Also exports formatThousands.

function formatThousands(value) {
    return Math.round(value).toString().replace(/\B(?=(\d{3})+(?!\d))/g, " ");
}

function Observer() {
    this.position = new THREE.Vector3(10,0,0);
    this.velocity = new THREE.Vector3(0,1,0);
    this.orientation = new THREE.Matrix3();
    this.time = 0.0;
}

Observer.prototype.orbitalFrame = function() {

    //var orbital_y = observer.velocity.clone().normalize();
    var orbital_y = (new THREE.Vector3())
        .subVectors(observer.velocity.clone().normalize().multiplyScalar(4.0),
            observer.position).normalize();

    var orbital_z = (new THREE.Vector3())
        .crossVectors(observer.position, orbital_y).normalize();
    var orbital_x = (new THREE.Vector3()).crossVectors(orbital_y, orbital_z);


    return (new THREE.Matrix4()).makeBasis(
        orbital_x,
        orbital_y,
        orbital_z
    ).linearPart();
};

Observer.prototype.move = function(dt) {

    dt *= shader.parameters.time_scale;

    var r;
    var v = 0;

    // motion on a pre-defined cirular orbit
    if (shader.parameters.observer.motion) {

        r = shader.parameters.observer.distance;
        v =  1.0 / Math.sqrt(2.0*(r-1.0));
        // Convert local velocity to coordinate angular velocity
        // Ω = v·sqrt(1-r_s/r)/r = 1/sqrt(2r³) for circular Schwarzschild orbit
        var ang_vel = v * Math.sqrt(1.0 - 1.0/r) / r;
        var angle = this.time * ang_vel;

        var s = Math.sin(angle), c = Math.cos(angle);

        this.position.set(c*r, s*r, 0);
        this.velocity.set(-s*v, c*v, 0);

        var alpha = degToRad(shader.parameters.observer.orbital_inclination);
        var orbit_coords = (new THREE.Matrix4()).makeRotationY(alpha);

        this.position.applyMatrix4(orbit_coords);
        this.velocity.applyMatrix4(orbit_coords);
    }
    else {
        r = this.position.length();
    }

    if (shader.parameters.gravitational_time_dilation) {
        if (v > 0) {
            // Circular orbit: combined gravitational + kinematic time dilation
            // dτ/dt = sqrt(1 - 3M/r) = sqrt(1 - 3/(2r)) for r_s = 1, M = 0.5
            dt = dt / Math.sqrt(Math.max(1.0 - 1.5/r, 0.001));
        } else {
            // Stationary observer: gravitational time dilation only
            // dτ/dt = sqrt(1 - r_s/r) = sqrt(1 - 1/r)
            dt = dt / Math.sqrt(Math.max(1.0 - 1.0/r, 0.001));
        }
    }

    this.time += dt;
};
