// Role: Application entry point - fetches all GLSL shader shards in order and
//       concatenates them into a single Mustache template string, loads all
//       textures in parallel, then calls init() once everything is ready.

import { $, THREE } from './vendor.js';
import { SHADER_SHARDS } from './shaders/shader-shards.js';
import { init } from './core/renderer/renderer.js';
import { animate } from './core/renderer/render-loop.js';

var textures = {};
var glslSource = null;

function checkReady() {
    if (glslSource === null) return;
    for (var key in textures) {
        if (textures[key] === null) return;
    }
    $('#loader').hide();
    $('.initially-hidden').removeClass('initially-hidden');
    init(glslSource, textures);
    animate();
}

var texLoader = new THREE.TextureLoader();

function loadTexture(symbol, filename, interpolation) {
    textures[symbol] = null;
    texLoader.load(filename, function(tex) {
        tex.magFilter = interpolation;
        tex.minFilter = interpolation;
        textures[symbol] = tex;
        checkReady();
    });
}

loadTexture('galaxy', 'assets/img/milkyway.jpg', THREE.NearestFilter);
loadTexture('spectra', 'assets/img/spectra.png', THREE.LinearFilter);
loadTexture('moon', 'assets/img/beach-ball.png', THREE.LinearFilter);
loadTexture('stars', 'assets/img/stars.png', THREE.LinearFilter);

var shardRequests = SHADER_SHARDS.map(function(path) {
    return $.get(path);
});

$.when.apply($, shardRequests).done(function() {
    var args = Array.prototype.slice.call(arguments);
    var parts;
    if (shardRequests.length === 1) {
        parts = [args[0]];
    } else {
        parts = args.map(function(result) {
            return Array.isArray(result) ? result[0] : result;
        });
    }
    glslSource = parts.join('\n');
    checkReady();
}).fail(function(jqXHR, textStatus, errorThrown) {
    console.error('Failed to load GLSL shards:', textStatus, errorThrown);
});


