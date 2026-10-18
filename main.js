/* Hyperbolic Tiling Explorer
 * Renders regular {p,q} tilings of the hyperbolic plane on the Poincaré disk.
 * The fundamental polygon is reflected across its edges (geodesic circle
 * inversions); reflections composed recursively generate the whole tiling.
 */
"use strict";

const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");

const params = {
  p: 7, q: 3, depth: 5,
  palette: 0,
  colorFund: true,
  animate: false,
  showGeo: true,
};

const PALETTES = [
  ["#1b1035", "#ff9f43", "#ffd166", "#3a86ff"],
  ["#0e1a2b", "#4ecdc4", "#ff6b6b", "#ffe66d"],
  ["#101418", "#e63946", "#f1faee", "#a8dadc"],
  ["#14060f", "#ff8fab", "#fb6f92", "#ffc2d1"],
  ["#0a0f0a", "#95d5b2", "#70e000", "#ffffff"],
  ["#0b0e14", "#7f5af0", "#2cb67d", "#e0e0e0"],
];

/* ---------- complex helpers ---------- */
const cx = (re, im) => ({ re, im });
const cmul = (a, b) => ({ re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re });
const cdiv = (a, b) => {
  const d = b.re * b.re + b.im * b.im;
  return { re: (a.re * b.re + a.im * b.im) / d, im: (a.im * b.re - a.re * b.im) / d };
};
const cadd = (a, b) => ({ re: a.re + b.re, im: a.im + b.im });

/* ---------- Möbius transforms (a,b,c,d): z ↦ (az+b)/(cz+d) ---------- */
function mobApply(m, z) {
  const num = cadd(cmul(m[0], z), m[1]);
  const den = cadd(cmul(m[2], z), m[3]);
  return cdiv(num, den);
}
function mobMul(m, n) {
  return [
    cadd(cmul(m[0], n[0]), cmul(m[1], n[2])),
    cadd(cmul(m[0], n[1]), cmul(m[1], n[3])),
    cadd(cmul(m[2], n[0]), cmul(m[3], n[2])),
    cadd(cmul(m[2], n[1]), cmul(m[3], n[3])),
  ];
}
const IDENT = [cx(1, 0), cx(0, 0), cx(0, 0), cx(1, 0)];

/* Transform = (m, c): z ↦ m(conj^c z) — conjugation on the INPUT.
   A geodesic reflection R(z) = o + R²/conj(z − o) equals F*(conj z) with
   F(z) = (ō z − 1)/(z − o), so R = {m: F*, c: 1} — a single transform.
   Composition t2∘t1 (t1 first), derived & verified numerically:
     c2 = 0 → {m: m2∘m1,        c: c1}
     c2 = 1 → {m: m2∘conj(m1),  c: 1 − c1}   (conj(m1) = coeff-conjugated)
   Apply: conj the input when c=1; the output is always holomorphic. */
function conjCoeffs(m) {
  return m.map(w => cx(w.re, -w.im));
}
function tCompose(t2, t1) { // apply t1 first, then t2
  return t2.c === 0
    ? { m: mobMul(t2.m, t1.m), c: t1.c }
    : { m: mobMul(t2.m, conjCoeffs(t1.m)), c: 1 - t1.c };
}
function tApply(t, z) {
  return mobApply(t.m, t.c ? cx(z.re, -z.im) : z);
}

/* ---------- fundamental polygon ---------- */
function computeFundamental(p, q, rot) {
  const s = Math.sin(Math.PI / p) * Math.sin(Math.PI / q);
  const c = Math.cos(Math.PI / p) * Math.cos(Math.PI / q);
  const r = Math.sqrt((c - s) / (c + s)); // Euclidean circumradius in the disk
  const verts = [];
  for (let i = 0; i < p; i++) {
    const th = 2 * Math.PI * i / p + rot;
    verts.push(cx(r * Math.cos(th), r * Math.sin(th)));
  }
  return { p, q, verts };
}

/* Geodesic (arc orthogonal to the unit circle) through two disk points.
   Returns {ox, oy, R} or null if the geodesic is a straight diameter. */
function geodesicCircle(z1, z2) {
  const x1 = z1.re, y1 = z1.im, x2 = z2.re, y2 = z2.im;
  const dx = x2 - x1, dy = y2 - y1;
  const det = dx * y1 - dy * x1;
  if (Math.abs(det) < 1e-12) return null;
  // center o satisfies:  o·(z2 - z1) = (|z2|² - |z1|²)/2  and  2 o·z1 = |z1|² + 1
  const rhs = (x2 * x2 + y2 * y2 - x1 * x1 - y1 * y1) / 2;
  const r1 = (x1 * x1 + y1 * y1 + 1) / 2;
  const ox = (rhs * y1 - r1 * dy) / det;
  const oy = (r1 * dx - rhs * x1) / det;
  const R2 = ox * ox + oy * oy - 1;
  if (R2 <= 0) return null;
  return { ox, oy, R: Math.sqrt(R2) };
}

/* Reflection across a geodesic circle (|o|² = R² + 1):
   z ↦ o + R² / conj(z − o) = F*(conj z), where F(z) = (ō·z − 1)/(z − o)
   and F* is F with conjugated coefficients. Returns F* — use as
   transform {m: reflectMat(g), c: 1}. */
function reflectMat(g) {
  return [cx(g.ox, g.oy), cx(-1, 0), cx(1, 0), cx(-g.ox, g.oy)];
}
function reflectMatLine(ang) {
  // reflection across the diameter at angle ang: z ↦ e^{2i·ang} · conj(z)
  // = e^{2i·ang} applied to conj(z) → {m: [e,0,0,1], c: 1}
  const e = cx(Math.cos(2 * ang), Math.sin(2 * ang));
  return [e, cx(0, 0), cx(0, 0), cx(1, 0)];
}

/* ---------- rendering ---------- */
const MAX_TILES = 20000;
let tileCount = 0;

function drawTiling() {
  const { p, q } = params;
  const fund = computeFundamental(p, q, rotAngle);
  const W = canvas.width, H = canvas.height;
  const scale = Math.min(W, H) * 0.48;
  const ox = W / 2, oy = H / 2;
  const pal = PALETTES[params.palette];

  ctx.fillStyle = "#0b0e14";
  ctx.fillRect(0, 0, W, H);
  ctx.beginPath();
  ctx.arc(ox, oy, scale, 0, 2 * Math.PI);
  ctx.fillStyle = pal[0];
  ctx.fill();

  tileCount = 0;

  // Edge reflections of the central polygon
  const edgeRefl = [];
  for (let i = 0; i < p; i++) {
    const v1 = fund.verts[i], v2 = fund.verts[(i + 1) % p];
    const g = geodesicCircle(v1, v2);
    edgeRefl.push(g
      ? { m: reflectMat(g), c: 1 }
      : { m: reflectMatLine(Math.atan2(v2.im, v2.re)), c: 1 });
  }

  const drawn = new Set();

  const polyPts = (t) =>
    fund.verts.map(v => {
      const w = tApply(t, v);
      return [ox + w.re * scale, oy + w.im * scale];
    });

  // Exact-ish tile identity: all vertices rounded to 1e-9 (in disk units).
  // Centroid keys are unsafe here — distinct tiles can share a centroid.
  const tileKey = (t) => fund.verts.map(v => {
    const w = tApply(t, v);
    return w.re.toFixed(9) + "," + w.im.toFixed(9);
  }).join(";");

  const renderPoly = (t, pts, depth, key) => {
    drawn.add(key);
    tileCount++;

    ctx.beginPath();
    ctx.moveTo(pts[0][0], pts[0][1]);
    for (let i = 0; i < p; i++) {
      const a1 = fund.verts[i], a2 = fund.verts[(i + 1) % p];
      const p2 = pts[(i + 1) % p];
      if (!params.showGeo) {
        ctx.lineTo(p2[0], p2[1]);
      } else {
        // arc through transformed midpoint of the edge
        const mid = cx((a1.re + a2.re) / 2, (a1.im + a2.im) / 2);
        const wm = tApply(t, mid);
        ctx.quadraticCurveTo(ox + wm.re * scale, oy + wm.im * scale, p2[0], p2[1]);
      }
    }
    ctx.closePath();

    if (params.colorFund) {
      const tCol = Math.min(depth / Math.max(params.depth, 1), 1);
      ctx.fillStyle = depth === 0 ? pal[3] : lerpColor(pal[1], pal[2], tCol);
    } else {
      ctx.fillStyle = depth % 2 ? pal[1] : pal[2];
    }
    ctx.globalAlpha = 0.92;
    ctx.fill();
    ctx.globalAlpha = 1;
    if (params.showGeo) {
      ctx.strokeStyle = "rgba(10,12,18,0.55)";
      ctx.lineWidth = 1;
      ctx.stroke();
    }
  };

  const recurse = (t, depth) => {
    if (depth > params.depth || tileCount > MAX_TILES) return;
    const key = tileKey(t);
    if (drawn.has(key)) return;
    const pts = polyPts(t);
    renderPoly(t, pts, depth, key);
    if (depth === params.depth) return;
    for (let i = 0; i < p; i++) {
      recurse(tCompose(t, edgeRefl[i]), depth + 1);
    }
  };

  recurse({ m: IDENT, c: 0 }, 0);

  ctx.beginPath();
  ctx.arc(ox, oy, scale, 0, 2 * Math.PI);
  ctx.strokeStyle = "rgba(255,255,255,0.15)";
  ctx.lineWidth = 2;
  ctx.stroke();

  document.getElementById("info").textContent =
    `{${p},${q}} — ${tileCount.toLocaleString()} tiles` +
    (tileCount > MAX_TILES ? " (limit reached — reduce depth)" : "");
}

function lerpColor(a, b, t) {
  const pa = [1, 3, 5].map(i => parseInt(a.substr(i, 2), 16));
  const pb = [1, 3, 5].map(i => parseInt(b.substr(i, 2), 16));
  const c = pa.map((v, i) => Math.round(v + (pb[i] - v) * t));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}

/* ---------- interaction: drag to spin ---------- */
let rotAngle = 0;
let dragging = false, lastX = 0;

/* ---------- UI ---------- */
function fitCanvas() {
  const dpr = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  canvas.width = rect.width * dpr;
  canvas.height = rect.height * dpr;
  drawTiling();
}
window.addEventListener("resize", fitCanvas);

function bindRange(id, key, labelId) {
  const el = document.getElementById(id);
  el.addEventListener("input", () => {
    params[key] = +el.value;
    if (labelId) document.getElementById(labelId).textContent = el.value;
    if (key === "p" || key === "q") {
      // valid hyperbolic tiling requires 1/p + 1/q < 1/2
      let { p, q } = params;
      while (1 / p + 1 / q >= 0.5 && q > 3) q--;
      while (1 / p + 1 / q >= 0.5 && p > 3) p--;
      if (p !== params.p) { params.p = p; document.getElementById("p").value = p; }
      if (q !== params.q) { params.q = q; document.getElementById("q").value = q; }
      document.getElementById("pVal").textContent = params.p;
      document.getElementById("qVal").textContent = params.q;
    }
    drawTiling();
  });
}
bindRange("p", "p", "pVal");
bindRange("q", "q", "qVal");
bindRange("depth", "depth", "depthVal");

document.getElementById("colorFund").addEventListener("change", e => { params.colorFund = e.target.checked; drawTiling(); });
document.getElementById("animate").addEventListener("change", e => { params.animate = e.target.checked; });
document.getElementById("showGeo").addEventListener("change", e => { params.showGeo = e.target.checked; drawTiling(); });

document.getElementById("randomize").addEventListener("click", () => {
  const p = 3 + Math.floor(Math.random() * 10);
  const qMin = Math.ceil(2 * p / (p - 2)); // smallest q giving a hyperbolic tiling
  const q = Math.min(12, qMin + Math.floor(Math.random() * Math.max(1, 13 - qMin)));
  params.p = p; params.q = q;
  document.getElementById("p").value = p;
  document.getElementById("q").value = q;
  document.getElementById("pVal").textContent = p;
  document.getElementById("qVal").textContent = q;
  const d = p * q > 40 ? 4 : 5; // keep render cost sane
  params.depth = d;
  document.getElementById("depth").value = d;
  document.getElementById("depthVal").textContent = d;
  drawTiling();
});

document.getElementById("save").addEventListener("click", () => {
  const a = document.createElement("a");
  a.download = `hyperbolic-tiling-${params.p}-${params.q}.png`;
  a.href = canvas.toDataURL("image/png");
  a.click();
});

const palBox = document.getElementById("palettes");
PALETTES.forEach((pal, i) => {
  const sw = document.createElement("div");
  sw.className = "swatch" + (i === 0 ? " active" : "");
  sw.style.background = `linear-gradient(135deg, ${pal[1]}, ${pal[2]}, ${pal[3]})`;
  sw.addEventListener("click", () => {
    params.palette = i;
    document.querySelectorAll(".swatch").forEach(s => s.classList.remove("active"));
    sw.classList.add("active");
    drawTiling();
  });
  palBox.appendChild(sw);
});

canvas.addEventListener("pointerdown", e => {
  dragging = true; lastX = e.clientX; canvas.style.cursor = "grabbing";
});
window.addEventListener("pointerup", () => { dragging = false; canvas.style.cursor = "grab"; });
window.addEventListener("pointermove", e => {
  if (dragging) {
    rotAngle += (e.clientX - lastX) * 0.01;
    lastX = e.clientX;
    drawTiling();
  }
});

let animT = 0;
function loop() {
  if (params.animate && !dragging) {
    animT += 0.01;
    rotAngle = animT;
    drawTiling();
  }
  requestAnimationFrame(loop);
}

fitCanvas();
loop();
