export function setupObserverOrbitWidget(options) {
    var parameters = options.parameters;
    var observerOrbitMin = options.observerOrbitMin;
    var observerDistanceMin = options.observerDistanceMin;
    var observerDistanceMax = options.observerDistanceMax;
    var clampObserverDistance = options.clampObserverDistance;
    var updateCamera = options.updateCamera;
    var markShaderDirty = options.markShaderDirty;
    var updateShader = options.updateShader || markShaderDirty;
    var resetObserverCamera = options.resetObserverCamera;
    var showObserverControlHint = options.showObserverControlHint;
    var registerBlackHoleUiBinding = options.registerBlackHoleUiBinding;
    var setDistanceController = options.setDistanceController;

    var gizmo = document.createElement('div');
    gizmo.id = 'axes-gizmo-container';
    gizmo.innerHTML =
        '<div id="observer-orbit-shell" class="observer-orbit-shell">' +
            '<svg id="observer-distance-arc" class="observer-distance-arc" viewBox="0 0 140 140" aria-label="Observer distance dial">' +
                '<path id="observer-distance-track" class="observer-distance-track"></path>' +
                '<path id="observer-distance-fill" class="observer-distance-fill"></path>' +
                '<path id="observer-distance-hit" class="observer-distance-hit"></path>' +
                '<circle id="observer-distance-knob" class="observer-distance-knob" r="6"></circle>' +
            '</svg>' +
            '<button id="observer-orbit-motion" type="button" class="observer-orbit-btn observer-orbit-btn-motion" aria-pressed="false">MOTION</button>' +
            '<button id="observer-orbit-reset" type="button" class="observer-orbit-btn observer-orbit-btn-reset">RESET</button>' +
            '<div id="observer-orbit-distance" class="observer-orbit-distance">11.0 r_s</div>' +
            '<div class="observer-orbit-help">L orbit | R pan | L+R roll</div>' +
            '<canvas id="axes-gizmo" width="80" height="80" aria-label="XYZ orientation indicator"></canvas>' +
        '</div>' +
        '<div class="observer-orbit-tag">OBSERVER</div>';
    document.body.appendChild(gizmo);

    var orbitShell = document.getElementById('observer-orbit-shell');
    var motionBtn = document.getElementById('observer-orbit-motion');
    var resetBtn = document.getElementById('observer-orbit-reset');
    var distanceReadout = document.getElementById('observer-orbit-distance');
    var distanceArcSvg = document.getElementById('observer-distance-arc');
    var distanceTrack = document.getElementById('observer-distance-track');
    var distanceFill = document.getElementById('observer-distance-fill');
    var distanceHit = document.getElementById('observer-distance-hit');
    var distanceKnob = document.getElementById('observer-distance-knob');

    var ARC = { cx: 70, cy: 70, radius: 56, startDeg: 60, endDeg: 300 };
    var draggingDistance = false;

    function clamp(value, minValue, maxValue) {
        return Math.max(minValue, Math.min(maxValue, value));
    }

    function setOrbitMenuOpen(isOpen) {
        orbitShell.classList.toggle('is-open', !!isOpen);
    }

    function getObserverDistanceMin() {
        return parameters.observer.motion ? observerOrbitMin : observerDistanceMin;
    }

    function distanceToProgress(distance) {
        var minDistance = getObserverDistanceMin();
        return (distance - minDistance) / (observerDistanceMax - minDistance);
    }

    function progressToDistance(progress) {
        var minDistance = getObserverDistanceMin();
        return minDistance + clamp(progress, 0.0, 1.0) * (observerDistanceMax - minDistance);
    }

    function arcPoint(angleDeg) {
        var radians = (angleDeg - 90.0) * Math.PI / 180.0;
        return {
            x: ARC.cx + ARC.radius * Math.cos(radians),
            y: ARC.cy + ARC.radius * Math.sin(radians)
        };
    }

    function describeArc(startDeg, endDeg) {
        var start = arcPoint(startDeg);
        var end = arcPoint(endDeg);
        var largeArcFlag = Math.abs(endDeg - startDeg) > 180 ? 1 : 0;
        var sweepFlag = endDeg >= startDeg ? 1 : 0;
        return 'M ' + start.x.toFixed(3) + ' ' + start.y.toFixed(3) +
            ' A ' + ARC.radius + ' ' + ARC.radius + ' 0 ' + largeArcFlag + ' ' + sweepFlag +
            ' ' + end.x.toFixed(3) + ' ' + end.y.toFixed(3);
    }

    function updateDistanceArc(distance) {
        var progress = clamp(distanceToProgress(distance), 0.0, 1.0);
        var angle = ARC.startDeg + progress * (ARC.endDeg - ARC.startDeg);
        var fillEnd = angle;
        if (Math.abs(fillEnd - ARC.startDeg) < 0.001) fillEnd += 0.001;

        distanceTrack.setAttribute('d', describeArc(ARC.startDeg, ARC.endDeg));
        distanceFill.setAttribute('d', describeArc(ARC.startDeg, fillEnd));
        distanceHit.setAttribute('d', describeArc(ARC.startDeg, ARC.endDeg));

        var knobPos = arcPoint(angle);
        distanceKnob.setAttribute('cx', knobPos.x.toFixed(3));
        distanceKnob.setAttribute('cy', knobPos.y.toFixed(3));
    }

    function setDistanceFromPointerEvent(event) {
        var rect = distanceArcSvg.getBoundingClientRect();
        if (!rect.width || !rect.height) return;

        var px = (event.clientX - rect.left) * (140.0 / rect.width);
        var py = (event.clientY - rect.top) * (140.0 / rect.height);
        var angle = (Math.atan2(py - ARC.cy, px - ARC.cx) * 180.0 / Math.PI + 90.0 + 360.0) % 360.0;
        angle = clamp(angle, ARC.startDeg, ARC.endDeg);

        var progress = (angle - ARC.startDeg) / (ARC.endDeg - ARC.startDeg);
        parameters.observer.distance = progressToDistance(progress);
        updateCamera();
        markShaderDirty();
        updateObserverWidget();
    }

    function updateObserverWidget() {
        if (!motionBtn || !distanceReadout) return;
        var motionEnabled = !!parameters.observer.motion;
        motionBtn.classList.toggle('is-active', motionEnabled);
        motionBtn.setAttribute('aria-pressed', motionEnabled ? 'true' : 'false');
        distanceReadout.textContent = parameters.observer.distance.toFixed(1) + ' r_s';
        updateDistanceArc(parameters.observer.distance);
    }

    var observerDistanceBinding = { updateDisplay: updateObserverWidget };
    setDistanceController(observerDistanceBinding);
    registerBlackHoleUiBinding('observerDistance', observerDistanceBinding);
    updateObserverWidget();

    function toggleOrbitMenuFromEvent(event) {
        if (event) {
            event.stopPropagation();
            if (event.preventDefault) event.preventDefault();
        }
        setOrbitMenuOpen(!orbitShell.classList.contains('is-open'));
    }

    if (window.PointerEvent) {
        document.getElementById('axes-gizmo').addEventListener('pointerdown', toggleOrbitMenuFromEvent);
    } else {
        document.getElementById('axes-gizmo').addEventListener('touchstart', toggleOrbitMenuFromEvent, { passive: false });
        document.getElementById('axes-gizmo').addEventListener('click', toggleOrbitMenuFromEvent);
    }

    orbitShell.addEventListener('pointerdown', function(event) {
        event.stopPropagation();
    });

    document.addEventListener('pointerdown', function(event) {
        if (!orbitShell.contains(event.target)) setOrbitMenuOpen(false);
    });

    motionBtn.addEventListener('click', function() {
        parameters.observer.motion = !parameters.observer.motion;
        parameters.observer.distance = clampObserverDistance(
            parameters.observer.distance,
            parameters.observer.motion
        );
        updateCamera();
        updateShader();
        showObserverControlHint(
            parameters.observer.motion ? 'Stable circular orbit enabled.' : 'Stationary observer.'
        );
        updateObserverWidget();
        setOrbitMenuOpen(true);
    });

    distanceHit.addEventListener('pointerdown', function(event) {
        draggingDistance = true;
        setOrbitMenuOpen(true);
        setDistanceFromPointerEvent(event);
        if (distanceHit.setPointerCapture) distanceHit.setPointerCapture(event.pointerId);
        event.preventDefault();
    });
    distanceHit.addEventListener('pointermove', function(event) {
        if (!draggingDistance) return;
        setDistanceFromPointerEvent(event);
        event.preventDefault();
    });
    distanceHit.addEventListener('pointerup', function(event) {
        draggingDistance = false;
        if (distanceHit.releasePointerCapture) distanceHit.releasePointerCapture(event.pointerId);
    });
    distanceHit.addEventListener('pointercancel', function() {
        draggingDistance = false;
    });

    resetBtn.addEventListener('click', function() {
        resetObserverCamera();
        updateObserverWidget();
        setOrbitMenuOpen(true);
    });

    return {
        syncControls: updateObserverWidget
    };
}
