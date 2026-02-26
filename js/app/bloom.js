// Role: Bloom post-processing pass — threshold → multi-level mip-chain Gaussian
//       blur → weighted composite. Simulates optical diffraction / lens glow from
//       bright emission. Called once during init() to build the pass; the returned
//       bloomPass object exposes render() and resize() for use in the render loop.

function setupBloom() {
    var ppVertexShader = [
        'varying vec2 vUv;',
        'void main() {',
        '    vUv = uv;',
        '    gl_Position = vec4(position, 1.0);',
        '}'
    ].join('\n');

    var thresholdFS = [
        'uniform sampler2D tDiffuse;',
        'uniform float threshold;',
        'uniform float softKnee;',
        'varying vec2 vUv;',
        'void main() {',
        '    vec4 c = texture2D(tDiffuse, vUv);',
        '    float br = max(c.r, max(c.g, c.b));',
        '    float knee = threshold * softKnee;',
        '    float soft = br - threshold + knee;',
        '    soft = clamp(soft, 0.0, 2.0 * knee);',
        '    soft = soft * soft / (4.0 * knee + 0.00001);',
        '    float contrib = max(soft, br - threshold) / max(br, 0.00001);',
        '    gl_FragColor = vec4(c.rgb * contrib, 1.0);',
        '}'
    ].join('\n');

    // 9-tap separable Gaussian blur (sigma ~ 1.77)
    // Weights: Pascal row 8 / 256 — sum = 1.0
    var blurFS = [
        'uniform sampler2D tDiffuse;',
        'uniform vec2 direction;',
        'varying vec2 vUv;',
        'void main() {',
        '    vec4 sum = vec4(0.0);',
        '    sum += texture2D(tDiffuse, vUv - 4.0 * direction) * 0.01621622;',
        '    sum += texture2D(tDiffuse, vUv - 3.0 * direction) * 0.05405405;',
        '    sum += texture2D(tDiffuse, vUv - 2.0 * direction) * 0.12162162;',
        '    sum += texture2D(tDiffuse, vUv - 1.0 * direction) * 0.19459459;',
        '    sum += texture2D(tDiffuse, vUv) * 0.22702703;',
        '    sum += texture2D(tDiffuse, vUv + 1.0 * direction) * 0.19459459;',
        '    sum += texture2D(tDiffuse, vUv + 2.0 * direction) * 0.12162162;',
        '    sum += texture2D(tDiffuse, vUv + 3.0 * direction) * 0.05405405;',
        '    sum += texture2D(tDiffuse, vUv + 4.0 * direction) * 0.01621622;',
        '    gl_FragColor = sum;',
        '}'
    ].join('\n');

    var copyFS = [
        'uniform sampler2D tDiffuse;',
        'varying vec2 vUv;',
        'void main() {',
        '    gl_FragColor = texture2D(tDiffuse, vUv);',
        '}'
    ].join('\n');

    // Composite: weighted sum of 5 bloom mip levels added to original
    // Wider mip levels are attenuated by bloomRadius^level for natural PSF falloff
    var compositeFS = [
        'uniform sampler2D tDiffuse;',
        'uniform sampler2D tBloom0;',
        'uniform sampler2D tBloom1;',
        'uniform sampler2D tBloom2;',
        'uniform sampler2D tBloom3;',
        'uniform sampler2D tBloom4;',
        'uniform float bloomStrength;',
        'uniform float bloomRadius;',
        'varying vec2 vUv;',
        'void main() {',
        '    vec4 orig = texture2D(tDiffuse, vUv);',
        '    float r = bloomRadius;',
        '    float w0 = 1.0, w1 = r, w2 = r*r, w3 = r*r*r, w4 = r*r*r*r;',
        '    float wsum = w0 + w1 + w2 + w3 + w4;',
        '    vec4 bloom = (texture2D(tBloom0, vUv) * w0',
        '        + texture2D(tBloom1, vUv) * w1',
        '        + texture2D(tBloom2, vUv) * w2',
        '        + texture2D(tBloom3, vUv) * w3',
        '        + texture2D(tBloom4, vUv) * w4) / wsum;',
        '    gl_FragColor = vec4(orig.rgb + bloom.rgb * bloomStrength, 1.0);',
        '}'
    ].join('\n');

    var BLOOM_LEVELS = 5;
    var rtParams = {
        minFilter: THREE.LinearFilter,
        magFilter: THREE.LinearFilter,
        format: THREE.RGBAFormat
    };

    var ppScene = new THREE.Scene();
    var ppCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
    var ppGeom = new THREE.PlaneBufferGeometry(2, 2);
    var ppMesh = new THREE.Mesh(ppGeom);
    ppScene.add(ppMesh);

    var thresholdMat = new THREE.ShaderMaterial({
        uniforms: {
            tDiffuse: { type: 't', value: null },
            threshold: { type: 'f', value: 0.65 },
            softKnee: { type: 'f', value: 0.5 }
        },
        vertexShader: ppVertexShader,
        fragmentShader: thresholdFS,
        depthWrite: false,
        depthTest: false
    });

    var blurMat = new THREE.ShaderMaterial({
        uniforms: {
            tDiffuse: { type: 't', value: null },
            direction: { type: 'v2', value: new THREE.Vector2() }
        },
        vertexShader: ppVertexShader,
        fragmentShader: blurFS,
        depthWrite: false,
        depthTest: false
    });

    var copyMat = new THREE.ShaderMaterial({
        uniforms: {
            tDiffuse: { type: 't', value: null }
        },
        vertexShader: ppVertexShader,
        fragmentShader: copyFS,
        depthWrite: false,
        depthTest: false
    });

    var compositeMat = new THREE.ShaderMaterial({
        uniforms: {
            tDiffuse: { type: 't', value: null },
            tBloom0: { type: 't', value: null },
            tBloom1: { type: 't', value: null },
            tBloom2: { type: 't', value: null },
            tBloom3: { type: 't', value: null },
            tBloom4: { type: 't', value: null },
            bloomStrength: { type: 'f', value: 0.35 },
            bloomRadius: { type: 'f', value: 0.85 }
        },
        vertexShader: ppVertexShader,
        fragmentShader: compositeFS,
        depthWrite: false,
        depthTest: false
    });

    function createTargets(w, h) {
        var mainRT = new THREE.WebGLRenderTarget(w, h, rtParams);
        var mips = [], temps = [];
        for (var i = 0; i < BLOOM_LEVELS; i++) {
            var mw = Math.max(1, Math.floor(w / Math.pow(2, i + 1)));
            var mh = Math.max(1, Math.floor(h / Math.pow(2, i + 1)));
            mips.push(new THREE.WebGLRenderTarget(mw, mh, rtParams));
            temps.push(new THREE.WebGLRenderTarget(mw, mh, rtParams));
        }
        return { mainRT: mainRT, mips: mips, temps: temps };
    }

    var targets = createTargets(1, 1);

    var bloomPass = {
        BLOOM_LEVELS: BLOOM_LEVELS,
        ppScene: ppScene,
        ppCamera: ppCamera,
        ppMesh: ppMesh,
        thresholdMat: thresholdMat,
        blurMat: blurMat,
        copyMat: copyMat,
        compositeMat: compositeMat,
        mainRT: targets.mainRT,
        bloomMips: targets.mips,
        bloomTemp: targets.temps,

        resize: function(w, h) {
            this.mainRT.dispose();
            for (var i = 0; i < BLOOM_LEVELS; i++) {
                this.bloomMips[i].dispose();
                this.bloomTemp[i].dispose();
            }
            var t = createTargets(w, h);
            this.mainRT = t.mainRT;
            this.bloomMips = t.mips;
            this.bloomTemp = t.temps;
        },

        render: function(rdr, mainScene, mainCamera, params) {
            var bp = this;

            // 1. Render main scene → full-res render target
            rdr.render(mainScene, mainCamera, bp.mainRT, true);

            // 2. Brightness threshold → first mip (half res)
            bp.thresholdMat.uniforms.tDiffuse.value = bp.mainRT;
            bp.thresholdMat.uniforms.threshold.value = params.threshold;
            bp.ppMesh.material = bp.thresholdMat;
            rdr.render(bp.ppScene, bp.ppCamera, bp.bloomMips[0], true);

            // 3. Progressive downsample + blur for each mip level
            for (var i = 0; i < bp.BLOOM_LEVELS; i++) {
                var mip = bp.bloomMips[i];
                var tmp = bp.bloomTemp[i];
                var mw = mip.width;
                var mh = mip.height;

                // Downsample from previous blurred level (bilinear)
                if (i > 0) {
                    bp.copyMat.uniforms.tDiffuse.value = bp.bloomMips[i - 1];
                    bp.ppMesh.material = bp.copyMat;
                    rdr.render(bp.ppScene, bp.ppCamera, mip, true);
                }

                // Horizontal blur: mip → temp
                bp.blurMat.uniforms.tDiffuse.value = mip;
                bp.blurMat.uniforms.direction.value.set(1.0 / mw, 0);
                bp.ppMesh.material = bp.blurMat;
                rdr.render(bp.ppScene, bp.ppCamera, tmp, true);

                // Vertical blur: temp → mip
                bp.blurMat.uniforms.tDiffuse.value = tmp;
                bp.blurMat.uniforms.direction.value.set(0, 1.0 / mh);
                rdr.render(bp.ppScene, bp.ppCamera, mip, true);
            }

            // 4. Composite: original + weighted bloom levels → screen
            bp.compositeMat.uniforms.tDiffuse.value = bp.mainRT;
            for (var j = 0; j < bp.BLOOM_LEVELS; j++) {
                bp.compositeMat.uniforms['tBloom' + j].value = bp.bloomMips[j];
            }
            bp.compositeMat.uniforms.bloomStrength.value = params.strength;
            bp.compositeMat.uniforms.bloomRadius.value = params.radius;
            bp.ppMesh.material = bp.compositeMat;
            rdr.render(bp.ppScene, bp.ppCamera);
        }
    };

    return bloomPass;
}
