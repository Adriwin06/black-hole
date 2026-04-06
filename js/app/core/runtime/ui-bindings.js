import { getBlackHoleUiBinding } from './runtime-registry.js';
import {
    distanceController,
    refreshAllControllersGlobal
} from './runtime-state.js';

export function getTimelinePanelBinding() {
    return getBlackHoleUiBinding('timelinePanel') || null;
}

export function getRefreshControllersBinding() {
    var registeredRefresh = getBlackHoleUiBinding('refreshControllers');
    if (typeof registeredRefresh === 'function') return registeredRefresh;
    return refreshAllControllersGlobal;
}

export function refreshRendererUiBindings() {
    var refresh = getRefreshControllersBinding();
    if (typeof refresh === 'function') refresh();
}

export function getObserverDistanceBinding() {
    return getBlackHoleUiBinding('observerDistance') || distanceController;
}

export function updateObserverDistanceBinding() {
    var binding = getObserverDistanceBinding();
    if (binding && typeof binding.updateDisplay === 'function') {
        binding.updateDisplay();
    }
}


