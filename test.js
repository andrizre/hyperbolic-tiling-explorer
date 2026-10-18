// Headless verification of the tiling math extracted from main.js
const fs = require("fs");
const src = fs.readFileSync(__dirname + "/main.js", "utf8").replace('"use strict";', "");
const domShim = `
var document = {
  getElementById: (id) => id === "canvas" ? canvasShim : ({ addEventListener(){}, appendChild(){}, textContent: "", value: "", style: {} }),
  querySelectorAll: () => [],
  createElement: () => ({ style: {}, addEventListener(){}, click(){} }),
};
var canvasShim = {
  getContext: () => new Proxy({}, { get: () => () => {} }),
  getBoundingClientRect: () => ({ width: 800, height: 800 }),
  addEventListener(){}, style: {}, width: 800, height: 800,
};
var window = { addEventListener(){}, devicePixelRatio: 1 };
var requestAnimationFrame = () => {};
`;
const fn = new Function(domShim + src +
  "; return { computeFundamental, geodesicCircle, reflectMat, reflectMatLine, tCompose, tApply, IDENT, cx };");
const M = fn();

let fail = 0;
// Validate: every reflection fixes its edge endpoints, R∘R = identity, all vertices stay
// inside the unit disk, and BFS layer counts grow strictly monotonically.
for (const [p, q] of [[7, 3], [5, 4], [3, 7], [12, 3], [8, 3]]) {
  const fund = M.computeFundamental(p, q, 0);
  const r0 = Math.hypot(fund.verts[0].re, fund.verts[0].im);
  const edgeRefl = [];
  for (let i = 0; i < p; i++) {
    const v1 = fund.verts[i], v2 = fund.verts[(i + 1) % p];
    const g = M.geodesicCircle(v1, v2);
    edgeRefl.push(g
      ? { m: M.reflectMat(g), c: 1 }
      : { m: M.reflectMatLine(Math.atan2(v2.im, v2.re)), c: 1 });
  }
  for (let i = 0; i < p; i++) {
    for (const v of [fund.verts[i], fund.verts[(i + 1) % p]]) {
      const w = M.tApply(edgeRefl[i], v);
      if (Math.hypot(w.re - v.re, w.im - v.im) > 1e-9) {
        console.log(`FIXED-POINT FAIL at {${p},${q}} edge=${i}`);
        fail++;
      }
    }
    const rr = M.tCompose(edgeRefl[i], edgeRefl[i]);
    const probe = M.tApply(rr, M.cx(0.15, -0.25));
    if (Math.hypot(probe.re - 0.15, probe.im + 0.25) > 1e-9) {
      console.log(`INVOLUTION FAIL at {${p},${q}} edge=${i}`);
      fail++;
    }
  }
  const layerKey = (t) => fund.verts.map(v => {
    const w = M.tApply(t, v);
    return w.re.toFixed(9) + "," + w.im.toFixed(9);
  }).join(";");
  let layer = [{ m: M.IDENT, c: 0 }];
  const seen = new Set([layerKey(layer[0])]);
  const counts = [1];
  let maxAbs = 0;
  const DEPTH = 5;
  for (let d = 1; d <= DEPTH; d++) {
    const next = [];
    for (const t of layer) {
      for (let i = 0; i < p; i++) {
        const tn = M.tCompose(t, edgeRefl[i]);
        const k = layerKey(tn);
        if (seen.has(k)) continue;
        seen.add(k);
        next.push(tn);
        for (const v of fund.verts) {
          const w = M.tApply(tn, v);
          maxAbs = Math.max(maxAbs, Math.hypot(w.re, w.im));
        }
      }
    }
    layer = next;
    counts.push(layer.length);
  }
  const monotone = counts.every((c, i) => i === 0 || c > counts[i - 1]);
  const ok = maxAbs <= 1.0001 && r0 < 1 && monotone;
  if (!ok) fail++;
  console.log(`{${p},${q}}: layers=[${counts}] total=${seen.size}, maxVertexRadius=${maxAbs.toFixed(6)} ${ok ? "OK" : "FAIL"}`);
}
process.exit(fail ? 1 : 0);
