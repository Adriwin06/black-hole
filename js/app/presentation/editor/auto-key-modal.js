export function showAutoKeyChoiceModal(options) {
    var t = options.time;
    var pathCount = options.pathCount;
    var esc = options.esc;
    var onFull = options.onFull;
    var onDiff = options.onDiff;

    var overlay = document.createElement('div');
    overlay.className = 'tl-modal-overlay';
    var tLabel = t < 0.01 ? 't = 0' : 't = ' + t.toFixed(2) + 's';
    overlay.innerHTML =
        '<div class="tl-modal">' +
        '<div class="tl-modal-title">Auto Key &#8212; Initial State</div>' +
        '<p class="tl-modal-body">First Auto Key capture at <b>' + esc(tLabel) + '</b>.<br>' +
        'How should the initial state be recorded?</p>' +
        '<div class="tl-modal-choices">' +
        '<button class="tl-modal-btn tl-modal-btn--full" type="button">' +
        '<span class="tl-modal-btn-label">Full initial state</span>' +
        '<span class="tl-modal-btn-sub">Key every parameter now (' + pathCount + ' values) &mdash; guarantees a clean reset on play</span>' +
        '</button>' +
        '<button class="tl-modal-btn tl-modal-btn--diff" type="button">' +
        '<span class="tl-modal-btn-label">Changes only (diff mode)</span>' +
        '<span class="tl-modal-btn-sub">Only record what changes on the next press &mdash; lighter, manual approach</span>' +
        '</button>' +
        '</div>' +
        '<button class="tl-modal-btn tl-modal-btn--cancel" type="button">Cancel</button>' +
        '</div>';
    document.body.appendChild(overlay);

    function close() {
        if (overlay.parentNode) overlay.parentNode.removeChild(overlay);
    }

    overlay.querySelector('.tl-modal-btn--full').addEventListener('click', function() {
        close();
        onFull();
    });
    overlay.querySelector('.tl-modal-btn--diff').addEventListener('click', function() {
        close();
        onDiff();
    });
    overlay.querySelector('.tl-modal-btn--cancel').addEventListener('click', close);
    overlay.addEventListener('mousedown', function(e) {
        if (e.target === overlay) close();
    });
    overlay.addEventListener('keydown', function(e) {
        if (e.key === 'Escape') close();
    });
    overlay.setAttribute('tabindex', '-1');
    overlay.focus();
}
