export function setupAnimationsPanel(options) {
    var initialObserverDistance = options.initialObserverDistance;
    var diveState = options.diveState;
    var hoverState = options.hoverState;
    var startDive = options.startDive;
    var resetDive = options.resetDive;
    var seekDive = options.seekDive;
    var startHover = options.startHover;
    var resetHover = options.resetHover;
    var seekHover = options.seekHover;
    var toggleAnimationTimelineCapture = options.toggleAnimationTimelineCapture;
    var setAnimationTimelineCaptureCameraSmoothingEnabled =
        options.setAnimationTimelineCaptureCameraSmoothingEnabled;
    var updateAnimationTimelineCaptureUi = options.updateAnimationTimelineCaptureUi;
    var buildTimelinePanel = options.buildTimelinePanel;
    var getTimelinePanelBinding = options.getTimelinePanelBinding;

    var ANIM_PANEL_WIDTH_STORAGE_KEY = 'black-hole.anim-panel.width';
    var ANIM_PANEL_DEFAULT_WIDTH = 420;
    var ANIM_PANEL_MIN_WIDTH = 340;
    var ANIM_PANEL_MAX_WIDTH = 680;

    var panel = document.createElement('div');
    panel.id = 'anim-panel';
    panel.className = 'sp-panel sp-panel--left sp-panel--collapsed';
    panel.innerHTML =
        '<div class="sp-resize sp-resize--right"></div>' +
        '<div class="sp-header">' +
            '<span class="sp-title">ANIMATIONS</span>' +
            '<button id="anim-close-btn" class="sp-close" type="button" ' +
                'aria-label="Close Animations panel">&times;</button>' +
        '</div>' +
        '<div class="sp-content">' +
        '<div class="anim-section" id="dive-section">' +
            '<button class="anim-section-toggle" id="dive-section-toggle" ' +
                'type="button" aria-expanded="false" aria-controls="dive-section-body">' +
                '<span class="anim-section-arrow">&#9654;</span> FREEFALL DIVE' +
            '</button>' +
            '<div class="anim-section-body" id="dive-section-body">' +
                '<div class="dive-desc">Radial plunge into the black hole interior. ' +
                'Uses the Schwarzschild geodesic solver through the event horizon, ' +
                'with additional rendering simplifications for interior presentation.</div>' +
                '<button id="dive-start-btn" class="dive-btn dive-btn-start">' +
                    '\u25b6 START DIVE</button>' +
                '<div class="dive-control-row">' +
                    '<label>Fall speed</label>' +
                    '<input type="range" id="dive-speed" min="0.01" max="5.0" ' +
                        'step="0.01" value="1.0">' +
                    '<span id="dive-speed-val">1.0\u00d7</span>' +
                '</div>' +
                '<div class="dive-control-row dive-cinematic-row">' +
                    '<label for="dive-cinematic">Auto-speed</label>' +
                    '<input type="checkbox" id="dive-cinematic">' +
                    '<span class="dive-cinematic-hint">Slow near photon sphere &amp; horizon</span>' +
                '</div>' +
                '<div id="dive-horizon-track" class="dive-horizon-track">' +
                    '<div id="dive-horizon-bar" class="dive-horizon-fill outside"></div>' +
                    '<div class="dive-horizon-label">Event Horizon</div>' +
                '</div>' +
                '<div class="dive-readout">' +
                    '<div id="dive-radius" class="dive-metric">' +
                        'r = ' + initialObserverDistance.toFixed(2) + ' r<sub>s</sub></div>' +
                    '<div id="dive-velocity" class="dive-metric">v = 0.000 c</div>' +
                    '<div id="dive-status" class="dive-status ready">Ready</div>' +
                '</div>' +
                '<div class="anim-capture-block anim-capture-block--dive">' +
                    '<div class="anim-capture-head">' +
                        '<span class="anim-capture-label">Timeline capture</span>' +
                        '<span id="dive-capture-status" class="anim-capture-status">Idle</span>' +
                    '</div>' +
                    '<div class="anim-capture-hint">Record the live dive plus orbit/pan camera moves into the bottom timeline at the current playhead.</div>' +
                    '<label class="anim-capture-check" for="dive-capture-smooth">' +
                        '<input type="checkbox" id="dive-capture-smooth">' +
                        '<span>Smooth recorded camera</span>' +
                    '</label>' +
                    '<button id="dive-capture-btn" class="dive-btn dive-btn-capture">' +
                        '&#9679; RECORD TO TIMELINE</button>' +
                '</div>' +
                '<button id="dive-reset-btn" class="dive-btn dive-btn-reset" disabled>\u21ba RESET</button>' +
            '</div>' +
        '</div>' +
        '<div class="anim-section" id="hover-section">' +
            '<button class="anim-section-toggle" id="hover-section-toggle" ' +
                'type="button" aria-expanded="false" aria-controls="hover-section-body">' +
                '<span class="anim-section-arrow">&#9654;</span> HOVER APPROACH' +
            '</button>' +
            '<div class="anim-section-body" id="hover-section-body">' +
                '<div class="hover-desc">Stationary observer firing thrusters to hover ' +
                'at a fixed radius. Zero velocity, pure gravitational blueshift. ' +
                'Light from infinity gains energy falling into the potential well.</div>' +
                '<button id="hover-start-btn" class="hover-btn hover-btn-start">' +
                    '\u25b6 START HOVER</button>' +
                '<div class="hover-control-row">' +
                    '<label>Descent speed</label>' +
                    '<input type="range" id="hover-speed" min="0.01" max="2.0" ' +
                        'step="0.01" value="0.3">' +
                    '<span id="hover-speed-val">0.3\u00d7</span>' +
                '</div>' +
                '<div id="hover-horizon-track" class="hover-horizon-track">' +
                    '<div id="hover-horizon-bar" class="hover-horizon-fill normal"></div>' +
                    '<div class="hover-horizon-label">Event Horizon</div>' +
                '</div>' +
                '<div class="hover-readout">' +
                    '<div id="hover-radius" class="hover-metric">' +
                        'r = ' + initialObserverDistance.toFixed(2) + ' r<sub>s</sub></div>' +
                    '<div id="hover-blueshift" class="hover-metric">D<sub>grav</sub> = 1.00\u00d7</div>' +
                    '<div id="hover-accel" class="hover-metric">a = 0.00 c\u00b2/r<sub>s</sub></div>' +
                    '<div id="hover-status" class="hover-status ready">Ready</div>' +
                '</div>' +
                '<div class="anim-capture-block anim-capture-block--hover">' +
                    '<div class="anim-capture-head">' +
                        '<span class="anim-capture-label">Timeline capture</span>' +
                        '<span id="hover-capture-status" class="anim-capture-status">Idle</span>' +
                    '</div>' +
                    '<div class="anim-capture-hint">Record the live hover descent plus orbit/pan camera moves into the bottom timeline at the current playhead.</div>' +
                    '<label class="anim-capture-check" for="hover-capture-smooth">' +
                        '<input type="checkbox" id="hover-capture-smooth">' +
                        '<span>Smooth recorded camera</span>' +
                    '</label>' +
                    '<button id="hover-capture-btn" class="hover-btn hover-btn-capture">' +
                        '&#9679; RECORD TO TIMELINE</button>' +
                '</div>' +
                '<button id="hover-reset-btn" class="hover-btn hover-btn-reset" disabled>\u21ba RESET</button>' +
            '</div>' +
        '</div>' +
        '</div>';
    document.body.appendChild(panel);

    var animResizer = panel.querySelector('.sp-resize');
    var animResizing = false;
    var animStartX = 0;
    var animStartW = ANIM_PANEL_DEFAULT_WIDTH;

    function readAnimWidth() {
        try {
            var w = parseFloat(localStorage.getItem(ANIM_PANEL_WIDTH_STORAGE_KEY));
            return isFinite(w) ? w : null;
        } catch (e) {
            return null;
        }
    }

    function saveAnimWidth(w) {
        try {
            localStorage.setItem(ANIM_PANEL_WIDTH_STORAGE_KEY, String(Math.round(w)));
        } catch (e) {}
    }

    function applyAnimWidth(w, persist) {
        w = Math.round(Math.max(ANIM_PANEL_MIN_WIDTH, Math.min(ANIM_PANEL_MAX_WIDTH, w)));
        panel.style.width = w + 'px';
        if (persist) saveAnimWidth(w);
    }

    applyAnimWidth(readAnimWidth() || ANIM_PANEL_DEFAULT_WIDTH, false);

    animResizer.addEventListener('pointerdown', function(e) {
        if (e.button !== 0) return;
        animResizing = true;
        animStartX = e.clientX;
        animStartW = panel.getBoundingClientRect().width || ANIM_PANEL_DEFAULT_WIDTH;
        document.body.classList.add('sp-resizing');
        if (animResizer.setPointerCapture) animResizer.setPointerCapture(e.pointerId);
        document.addEventListener('pointermove', onAnimResizeMove);
        document.addEventListener('pointerup', onAnimResizeEnd);
        e.preventDefault();
    });

    function onAnimResizeMove(e) {
        if (!animResizing) return;
        applyAnimWidth(animStartW + (e.clientX - animStartX), false);
        e.preventDefault();
    }

    function onAnimResizeEnd() {
        if (!animResizing) return;
        animResizing = false;
        document.body.classList.remove('sp-resizing');
        applyAnimWidth(parseFloat(panel.style.width) || ANIM_PANEL_DEFAULT_WIDTH, true);
        document.removeEventListener('pointermove', onAnimResizeMove);
        document.removeEventListener('pointerup', onAnimResizeEnd);
    }

    var openBtn = document.createElement('button');
    openBtn.id = 'anim-open-btn';
    openBtn.className = 'sp-open-btn sp-open-btn--left';
    openBtn.type = 'button';
    openBtn.innerHTML = '&#9664; ANIMATIONS';
    openBtn.title = 'Open Animations panel';
    document.body.appendChild(openBtn);

    function setAnimPanelOpen(isOpen) {
        panel.classList.toggle('sp-panel--collapsed', !isOpen);
        openBtn.classList.toggle('sp-open-btn--hidden', isOpen);
        document.body.classList.toggle('has-anim-panel', isOpen);
    }

    setAnimPanelOpen(false);

    openBtn.addEventListener('click', function() { setAnimPanelOpen(true); });
    document.getElementById('anim-close-btn').addEventListener('click', function() {
        setAnimPanelOpen(false);
    });

    function setupSectionToggle(sectionId) {
        var section = document.getElementById(sectionId);
        var toggle = section.querySelector('.anim-section-toggle');
        toggle.addEventListener('click', function() {
            var isOpen = section.classList.toggle('is-open');
            toggle.setAttribute('aria-expanded', isOpen ? 'true' : 'false');
        });
    }

    setupSectionToggle('dive-section');
    setupSectionToggle('hover-section');

    document.getElementById('dive-start-btn').addEventListener('click', function() {
        startDive();
    });
    document.getElementById('dive-reset-btn').addEventListener('click', function() {
        resetDive();
    });
    document.getElementById('dive-speed').addEventListener('input', function() {
        diveState.speed = parseFloat(this.value);
        var v = diveState.speed;
        document.getElementById('dive-speed-val').textContent =
            (v < 0.1 ? v.toFixed(2) : v.toFixed(1)) + '\u00d7';
    });
    document.getElementById('dive-cinematic').addEventListener('change', function() {
        diveState.cinematic = this.checked;
    });
    document.getElementById('dive-capture-btn').addEventListener('click', function() {
        toggleAnimationTimelineCapture('dive');
    });
    document.getElementById('dive-capture-smooth').addEventListener('change', function() {
        setAnimationTimelineCaptureCameraSmoothingEnabled(this.checked);
    });

    var diveTrack = document.getElementById('dive-horizon-track');
    var diveDragging = false;

    function handleDiveTrackSeek(e) {
        if (!diveState.active && !diveState.reachedSingularity) return;
        var rect = diveTrack.getBoundingClientRect();
        var clientX = e.touches ? e.touches[0].clientX : e.clientX;
        var x = Math.max(0, Math.min(clientX - rect.left, rect.width));
        var progress = x / rect.width;
        var startR = Math.max(diveState.prevDistance, 1);
        var targetR = startR * (1.0 - progress);
        seekDive(targetR);
    }

    diveTrack.addEventListener('mousedown', function(e) {
        diveDragging = true;
        handleDiveTrackSeek(e);
        e.preventDefault();
    });
    document.addEventListener('mousemove', function(e) {
        if (diveDragging) handleDiveTrackSeek(e);
    });
    document.addEventListener('mouseup', function() {
        diveDragging = false;
    });
    diveTrack.addEventListener('touchstart', function(e) {
        handleDiveTrackSeek(e);
        e.preventDefault();
    });
    diveTrack.addEventListener('touchmove', function(e) {
        handleDiveTrackSeek(e);
        e.preventDefault();
    });

    document.getElementById('hover-start-btn').addEventListener('click', function() {
        startHover();
    });
    document.getElementById('hover-reset-btn').addEventListener('click', function() {
        resetHover();
    });
    document.getElementById('hover-speed').addEventListener('input', function() {
        hoverState.speed = parseFloat(this.value);
        var v = hoverState.speed;
        document.getElementById('hover-speed-val').textContent =
            (v < 0.1 ? v.toFixed(2) : v.toFixed(1)) + '\u00d7';
    });
    document.getElementById('hover-capture-btn').addEventListener('click', function() {
        toggleAnimationTimelineCapture('hover');
    });
    document.getElementById('hover-capture-smooth').addEventListener('change', function() {
        setAnimationTimelineCaptureCameraSmoothingEnabled(this.checked);
    });

    var hoverTrack = document.getElementById('hover-horizon-track');
    var hoverDragging = false;

    function handleHoverTrackSeek(e) {
        if (!hoverState.active) return;
        var rect = hoverTrack.getBoundingClientRect();
        var clientX = e.touches ? e.touches[0].clientX : e.clientX;
        var x = Math.max(0, Math.min(clientX - rect.left, rect.width));
        var progress = x / rect.width;
        var startR = Math.max(hoverState.prevDistance, 1);
        var targetR = startR * (1.0 - progress);
        targetR = Math.max(hoverState.minR, targetR);
        seekHover(targetR);
    }

    hoverTrack.addEventListener('mousedown', function(e) {
        hoverDragging = true;
        handleHoverTrackSeek(e);
        e.preventDefault();
    });
    document.addEventListener('mousemove', function(e) {
        if (hoverDragging) handleHoverTrackSeek(e);
    });
    document.addEventListener('mouseup', function() {
        hoverDragging = false;
    });
    hoverTrack.addEventListener('touchstart', function(e) {
        handleHoverTrackSeek(e);
        e.preventDefault();
    });
    hoverTrack.addEventListener('touchmove', function(e) {
        handleHoverTrackSeek(e);
        e.preventDefault();
    });

    updateAnimationTimelineCaptureUi();

    buildTimelinePanel();

    var tlOpenBtn = document.createElement('button');
    tlOpenBtn.id = 'tl-open-btn';
    tlOpenBtn.type = 'button';
    tlOpenBtn.textContent = '\u25B2 TIMELINE';
    tlOpenBtn.title = 'Open the full timeline / dopesheet panel';
    document.body.appendChild(tlOpenBtn);

    tlOpenBtn.addEventListener('click', function() {
        var timelineBinding = getTimelinePanelBinding();
        if (timelineBinding && typeof timelineBinding.toggle === 'function') {
            timelineBinding.toggle();
        }
    });
}
