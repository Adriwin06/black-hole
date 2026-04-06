export function lookAtOriginQuat(px, py, pz) {
    var r = Math.sqrt(px * px + py * py + pz * pz);
    if (r < 1e-8) return { x: 0, y: 0, z: 0, w: 1 };
    var zx = px / r;
    var zy = py / r;
    var zz = pz / r;
    var upx = 0;
    var upy = 1;
    var upz = 0;
    if (Math.abs(zy) > 0.999) {
        upx = 0;
        upy = 0;
        upz = 1;
    }
    var rx = upy * zz - upz * zy;
    var ry = upz * zx - upx * zz;
    var rz = upx * zy - upy * zx;
    var rlen = Math.sqrt(rx * rx + ry * ry + rz * rz);
    if (rlen < 1e-8) return { x: 0, y: 0, z: 0, w: 1 };
    rx /= rlen;
    ry /= rlen;
    rz /= rlen;
    var ux = zy * rz - zz * ry;
    var uy = zz * rx - zx * rz;
    var uz = zx * ry - zy * rx;
    var m00 = rx;
    var m01 = ux;
    var m02 = zx;
    var m10 = ry;
    var m11 = uy;
    var m12 = zy;
    var m20 = rz;
    var m21 = uz;
    var m22 = zz;
    var trace = m00 + m11 + m22;
    var qx;
    var qy;
    var qz;
    var qw;
    if (trace > 0) {
        var s = 0.5 / Math.sqrt(trace + 1.0);
        qw = 0.25 / s;
        qx = (m21 - m12) * s;
        qy = (m02 - m20) * s;
        qz = (m10 - m01) * s;
    } else if (m00 > m11 && m00 > m22) {
        var s2 = 2.0 * Math.sqrt(1.0 + m00 - m11 - m22);
        qw = (m21 - m12) / s2;
        qx = 0.25 * s2;
        qy = (m01 + m10) / s2;
        qz = (m02 + m20) / s2;
    } else if (m11 > m22) {
        var s3 = 2.0 * Math.sqrt(1.0 + m11 - m00 - m22);
        qw = (m02 - m20) / s3;
        qx = (m01 + m10) / s3;
        qy = 0.25 * s3;
        qz = (m12 + m21) / s3;
    } else {
        var s4 = 2.0 * Math.sqrt(1.0 + m22 - m00 - m11);
        qw = (m10 - m01) / s4;
        qx = (m02 + m20) / s4;
        qy = (m12 + m21) / s4;
        qz = 0.25 * s4;
    }
    return { x: qx, y: qy, z: qz, w: qw };
}

export function quatMul(q1, q2) {
    return {
        w: q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
        x: q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
        y: q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
        z: q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w
    };
}

export function quatSlerp(q1, q2, t) {
    var dot = q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
    if (dot < 0) {
        q2 = { x: -q2.x, y: -q2.y, z: -q2.z, w: -q2.w };
        dot = -dot;
    }
    if (dot > 0.9995) {
        var rx = q1.x + t * (q2.x - q1.x);
        var ry = q1.y + t * (q2.y - q1.y);
        var rz = q1.z + t * (q2.z - q1.z);
        var rw = q1.w + t * (q2.w - q1.w);
        var rlen = Math.sqrt(rx * rx + ry * ry + rz * rz + rw * rw);
        return { x: rx / rlen, y: ry / rlen, z: rz / rlen, w: rw / rlen };
    }
    var theta0 = Math.acos(dot);
    var sinTheta0 = Math.sin(theta0);
    var theta = theta0 * t;
    var s1 = Math.cos(theta) - dot * Math.sin(theta) / sinTheta0;
    var s2 = Math.sin(theta) / sinTheta0;
    return {
        x: s1 * q1.x + s2 * q2.x,
        y: s1 * q1.y + s2 * q2.y,
        z: s1 * q1.z + s2 * q2.z,
        w: s1 * q1.w + s2 * q2.w
    };
}

export function sampleTrackDraft(track, t) {
    var keys = track.keys;
    if (!keys || !keys.length) return 0;
    if (t <= keys[0].t) return +keys[0].v;
    if (t >= keys[keys.length - 1].t) return +keys[keys.length - 1].v;
    for (var i = 0; i < keys.length - 1; i++) {
        if (t >= keys[i].t && t <= keys[i + 1].t) {
            var dt = keys[i + 1].t - keys[i].t;
            if (dt < 1e-8) return +keys[i].v;
            return +keys[i].v + (+keys[i + 1].v - +keys[i].v) * ((t - keys[i].t) / dt);
        }
    }
    return +keys[keys.length - 1].v;
}
