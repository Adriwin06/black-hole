export function setupControlsPanel(gui) {
    var STORED_WIDTH_KEY = 'black-hole.controls-panel.width';
    var DEFAULT_WIDTH = 370;
    var MIN_WIDTH = 300;
    var MAX_WIDTH = 600;

    var ctrlPanel = document.createElement('div');
    ctrlPanel.id = 'controls-panel';
    ctrlPanel.className = 'sp-panel sp-panel--right sp-panel--collapsed';
    ctrlPanel.innerHTML =
        '<div class="sp-resize sp-resize--left"></div>' +
        '<div class="sp-header">' +
            '<span class="sp-title">CONTROLS</span>' +
            '<button id="controls-close-btn" class="sp-close" type="button">&times;</button>' +
        '</div>' +
        '<div id="controls-content" class="sp-content"></div>';
    document.body.appendChild(ctrlPanel);

    var dgAc = gui.domElement.parentElement;
    document.getElementById('controls-content').appendChild(dgAc);

    var openBtn = document.createElement('button');
    openBtn.id = 'controls-open-btn';
    openBtn.className = 'sp-open-btn sp-open-btn--right';
    openBtn.type = 'button';
    openBtn.innerHTML = 'CONTROLS &#9654;';
    openBtn.title = 'Open Controls panel';
    document.body.appendChild(openBtn);

    var panelOpen = false;
    function setPanelOpen(open) {
        panelOpen = open;
        ctrlPanel.classList.toggle('sp-panel--collapsed', !open);
        openBtn.classList.toggle('sp-open-btn--hidden', open);
        document.body.classList.toggle('has-controls-panel', open);
    }

    openBtn.addEventListener('click', function() { setPanelOpen(true); });
    document.getElementById('controls-close-btn')
        .addEventListener('click', function() { setPanelOpen(false); });

    var resizer = ctrlPanel.querySelector('.sp-resize');
    var resizing = false;
    var startX = 0;
    var startW = DEFAULT_WIDTH;

    function readWidth() {
        try {
            var width = parseFloat(localStorage.getItem(STORED_WIDTH_KEY));
            return isFinite(width) ? width : null;
        } catch (err) {
            return null;
        }
    }

    function saveWidth(width) {
        try {
            localStorage.setItem(STORED_WIDTH_KEY, String(Math.round(width)));
        } catch (err) {}
    }

    function applyWidth(width, persist) {
        width = Math.round(Math.max(MIN_WIDTH, Math.min(MAX_WIDTH, width)));
        ctrlPanel.style.width = width + 'px';
        if (persist) saveWidth(width);
    }

    applyWidth(readWidth() || DEFAULT_WIDTH, false);

    resizer.addEventListener('pointerdown', function(e) {
        if (e.button !== 0) return;
        resizing = true;
        startX = e.clientX;
        startW = ctrlPanel.getBoundingClientRect().width || DEFAULT_WIDTH;
        document.body.classList.add('sp-resizing');
        if (resizer.setPointerCapture) resizer.setPointerCapture(e.pointerId);
        document.addEventListener('pointermove', onMove);
        document.addEventListener('pointerup', onEnd);
        e.preventDefault();
    });

    function onMove(e) {
        if (!resizing) return;
        applyWidth(startW + (startX - e.clientX), false);
        e.preventDefault();
    }

    function onEnd() {
        if (!resizing) return;
        resizing = false;
        document.body.classList.remove('sp-resizing');
        applyWidth(parseFloat(ctrlPanel.style.width) || DEFAULT_WIDTH, true);
        document.removeEventListener('pointermove', onMove);
        document.removeEventListener('pointerup', onEnd);
    }

    if (!(window.matchMedia && window.matchMedia('(max-width: 960px)').matches)) {
        setPanelOpen(true);
    }

    return {
        setOpen: setPanelOpen,
        isOpen: function() { return panelOpen; }
    };
}
