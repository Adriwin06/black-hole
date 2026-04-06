"use strict";

var app = window.BLACK_HOLE_APP || (window.BLACK_HOLE_APP = {});

app.runtimeApis = app.runtimeApis || {};
app.uiBindings = app.uiBindings || {};

export function registerBlackHoleRuntimeApi(name, api) {
    if (!name) return api || null;
    app.runtimeApis[name] = api || null;
    return app.runtimeApis[name];
}

export function getBlackHoleRuntimeApi(name) {
    if (!name) return null;
    return app.runtimeApis[name] || null;
}

export function registerBlackHoleUiBinding(name, binding) {
    if (!name) return binding || null;
    app.uiBindings[name] = binding || null;
    return app.uiBindings[name];
}

export function getBlackHoleUiBinding(name) {
    if (!name) return null;
    return app.uiBindings[name] || null;
}

window.registerBlackHoleRuntimeApi = registerBlackHoleRuntimeApi;
window.getBlackHoleRuntimeApi = getBlackHoleRuntimeApi;
window.registerBlackHoleUiBinding = registerBlackHoleUiBinding;
window.getBlackHoleUiBinding = getBlackHoleUiBinding;
