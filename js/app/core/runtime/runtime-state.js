import { THREE } from '../../vendor.js';
import { Observer, bindObserverShader } from '../observer.js';

export const DISK_TEMPERATURE_MIN = 4500.0;
export const DISK_TEMPERATURE_MAX = 30000.0;

export let container = null;
export let stats = null;
export let camera = null;
export let scene = null;
export let renderer = null;
export let cameraControls = null;
export let shader = null;
export const observer = new Observer();
export const cameraPan = new THREE.Vector2(0, 0);
export let distanceController = null;
export let refreshAllControllersGlobal = null;
export let bloomPass = null;
export let taaPass = null;
export let shaderUniforms = null;
export const baseDevicePixelRatio = Math.max(window.devicePixelRatio || 1.0, 1.0);
export let isMobileClient = false;
export let rendererContextLost = false;
export const lastTaaCameraMat = new THREE.Matrix4().identity();

export let updateUniforms = function() {};
export let applyRenderScaleFromSettings = function() {};
export let resetTemporalAAHistory = function() {};

export function setContainer(value) { container = value; }
export function setStats(value) { stats = value; }
export function setCamera(value) { camera = value; }
export function setScene(value) { scene = value; }
export function setRenderer(value) { renderer = value; }
export function setCameraControls(value) { cameraControls = value; }
export function setShader(value) {
    shader = value;
    bindObserverShader(value);
}
export function setDistanceController(value) { distanceController = value; }
export function setRefreshAllControllersGlobal(value) { refreshAllControllersGlobal = value; }
export function setBloomPass(value) { bloomPass = value; }
export function setTaaPass(value) { taaPass = value; }
export function setShaderUniforms(value) { shaderUniforms = value; }
export function setIsMobileClient(value) { isMobileClient = !!value; }
export function setRendererContextLost(value) { rendererContextLost = !!value; }
export function setUpdateUniforms(value) { updateUniforms = value || function() {}; }
export function setApplyRenderScaleFromSettings(value) {
    applyRenderScaleFromSettings = value || function() {};
}
export function setResetTemporalAAHistory(value) {
    resetTemporalAAHistory = value || function() {};
}


