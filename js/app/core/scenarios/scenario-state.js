import { THREE } from '../../vendor.js';

export const diveState = {
    active: false,
    paused: false,
    speed: 1.0,
    cinematic: false,
    autoOrient: true,
    timelineDriven: false,
    currentR: 11.0,
    direction: new THREE.Vector3(1, 0, 0),
    startPosition: new THREE.Vector3(10, 0, 0),
    startVelocity: new THREE.Vector3(0, 1, 0),
    startRenderSettings: null,
    prevMotionState: true,
    prevDistance: 11.0,
    reachedSingularity: false
};

export const hoverState = {
    active: false,
    paused: false,
    speed: 0.3,
    timelineDriven: false,
    currentR: 11.0,
    direction: new THREE.Vector3(1, 0, 0),
    startPosition: new THREE.Vector3(10, 0, 0),
    startVelocity: new THREE.Vector3(0, 1, 0),
    prevMotionState: true,
    prevDistance: 11.0,
    minR: 1.0002
};

export const animationTimelineCaptureState = {
    active: false,
    mode: '',
    startedAtMs: 0,
    lastElapsed: 0,
    nextSampleTime: 0,
    sampleInterval: 1.0 / 30.0,
    cameraSmoothingEnabled: (function() {
        try {
            return localStorage.getItem('black-hole.anim-capture.camera-smoothing') === '1';
        } catch (e) {
            return false;
        }
    })(),
    samples: [],
    startPosition: null,
    startVelocity: null,
    prevMotionState: true,
    prevDistance: 11.0,
    startObserverTime: 0.0,
    feedback: {
        dive: { text: 'Idle', tone: '' },
        hover: { text: 'Idle', tone: '' }
    }
};


