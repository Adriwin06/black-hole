// Role: Application entry point — fetches all GLSL shader shards in order and
//       concatenates them into a single Mustache template string, loads all
//       textures in parallel, then calls init() once everything is ready.
//       Replaces the old SHADER_LOADER single-file approach.

(function () {
    "use strict";

    // Ordered list of GLSL shards — concatenated top-to-bottom to form the
    // complete fragment shader before Mustache template expansion.
    var SHADER_SHARDS = [
        'shaders/raytracer/core/defines.glsl',
        'shaders/raytracer/core/math.glsl',
        'shaders/raytracer/physics/geodesics.glsl',
        'shaders/raytracer/physics/accretion.glsl',
        'shaders/raytracer/physics/jet.glsl',
        'shaders/raytracer/output/tonemapping.glsl',
        'shaders/raytracer/physics/planet.glsl',
        'shaders/raytracer/physics/background.glsl',
        'shaders/raytracer/output/trace_ray.glsl',
        'shaders/raytracer/output/main.glsl'
    ];

    var textures = {};
    var glslSource = null;
    var loadedCount = 0;
    var totalExpected = 0;

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

    // ── Texture loading ──────────────────────────────────────────────────────────
    var texLoader = new THREE.TextureLoader();

    function loadTexture(symbol, filename, interpolation) {
        textures[symbol] = null;
        totalExpected++;
        texLoader.load(filename, function(tex) {
            tex.magFilter = interpolation;
            tex.minFilter = interpolation;
            textures[symbol] = tex;
            checkReady();
        });
    }

    loadTexture('galaxy', 'img/milkyway.jpg',   THREE.NearestFilter);
    loadTexture('spectra', 'img/spectra.png',    THREE.LinearFilter);
    loadTexture('moon',   'img/beach-ball.png',  THREE.LinearFilter);
    loadTexture('stars',  'img/stars.png',       THREE.LinearFilter);

    // ── GLSL shard loading ───────────────────────────────────────────────────────
    // Fetch all shards in parallel then concatenate in declaration order.
    var shardRequests = SHADER_SHARDS.map(function(path) {
        return $.get(path);
    });

    $.when.apply($, shardRequests).done(function() {
        // $.when passes each result as a separate argument; each is [data, status, jqXHR]
        var args = Array.prototype.slice.call(arguments);
        // If only one shard, $.when passes the values directly (not as array)
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

}());

