'use strict';

/* ============================================================
 * Image Morphing — Beier–Neely 特征线形变
 * 纯原生 JS + Canvas 实现，无第三方依赖
 * ============================================================ */

// ---------- 参数 ----------
const MAX_WORK = 700;   // 工作分辨率最长边（像素）
const A = 0.01;         // Beier-Neely 参数 a：线与像素距离的衰减基准
const B = 2.0;          // Beier-Neely 参数 b：权重指数
const P = 0.5;          // Beier-Neely 参数 p：线长度影响
const HIT_RADIUS = 12;  // 鼠标命中半径（画布像素）
const LOOP_SECONDS = 4; // 自动播放单程时长（秒）
const MIN_LINE_LEN = 4; // 新建特征线的最小长度（图像像素）

// ---------- DOM ----------
const srcCanvas = document.getElementById('srcCanvas');
const tgtCanvas = document.getElementById('tgtCanvas');
const resCanvas = document.getElementById('resCanvas');
const srcCtx = srcCanvas.getContext('2d');
const tgtCtx = tgtCanvas.getContext('2d');
const resCtx = resCanvas.getContext('2d');
const tSlider = document.getElementById('tSlider');
const tLabel = document.getElementById('tLabel');
const playBtn = document.getElementById('playBtn');
const speedSel = document.getElementById('speedSel');
const lineCountEl = document.getElementById('lineCount');
const delBtn = document.getElementById('delBtn');
const clearBtn = document.getElementById('clearBtn');
const srcFile = document.getElementById('srcFile');
const tgtFile = document.getElementById('tgtFile');

// ---------- 状态 ----------
const state = {
  src: null,      // {canvas, w, h, data(ImageData)}
  tgt: null,
  lines: [],      // [{id, src:{x1,y1,x2,y2}, tgt:{...}}] 图像坐标
  t: 0.5,
  playing: false,
  playDir: 1,
  lastFrame: 0,
  selectedId: null,
  drag: null,     // {id, side, end:0|1}
  creating: null, // {id, side, from:{x,y}}
  nextId: 1,
};

// 复用的输出缓冲区
let scratch1 = null, scratch2 = null, scratchRes = null;

/* ============================================================
 * 图像加载
 * ============================================================ */

function makeWorkingCanvas(img) {
  const scale = Math.min(1, MAX_WORK / Math.max(img.naturalWidth, img.naturalHeight));
  const w = Math.max(1, Math.round(img.naturalWidth * scale));
  const h = Math.max(1, Math.round(img.naturalHeight * scale));
  const c = document.createElement('canvas');
  c.width = w;
  c.height = h;
  const ctx = c.getContext('2d');
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(img, 0, 0, w, h);
  return wrapCanvas(c);
}

function wrapCanvas(c) {
  return {
    canvas: c,
    w: c.width,
    h: c.height,
    data: c.getContext('2d').getImageData(0, 0, c.width, c.height),
  };
}

// 将 img 等比缩放入 w×h 的画布（留边透明）
function resizeTo(img, w, h) {
  const c = document.createElement('canvas');
  c.width = w;
  c.height = h;
  const ctx = c.getContext('2d');
  const scale = Math.min(w / img.w, h / img.h);
  const dw = img.w * scale, dh = img.h * scale;
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(img.canvas, (w - dw) / 2, (h - dh) / 2, dw, dh);
  return wrapCanvas(c);
}

// 保证源图与目标图工作分辨率一致（以源图尺寸为准）
function alignSizes() {
  if (!state.src || !state.tgt) return;
  if (state.tgt.w !== state.src.w || state.tgt.h !== state.src.h) {
    state.tgt = resizeTo(state.tgt, state.src.w, state.src.h);
  }
}

function loadImageFromSource(side) {
  const input = side === 'src' ? srcFile : tgtFile;
  const file = input.files && input.files[0];
  if (!file) return;
  const img = new Image();
  img.onload = () => {
    const im = makeWorkingCanvas(img);
    if (side === 'src') state.src = im;
    else state.tgt = im;
    alignSizes();
    state.lines = [];
    state.selectedId = null;
    state.drag = null;
    state.creating = null;
    updateLineUI();
    requestRender();
  };
  img.src = URL.createObjectURL(file);
}

/* ============================================================
 * Beier–Neely 核心算法
 * ============================================================ */

// 反向映射 warp：对输出图像每个像素 X，依据目标线(dstLines)与源线(srcLines)
// 求其在源图中的对应位置，再做双线性采样。
function warp(srcImg, srcLines, dstLines, outData) {
  const w = srcImg.w, h = srcImg.h;
  const sdata = srcImg.data.data;
  const odata = outData.data;
  const n = srcLines.length;

  // 无特征线时直接拷贝（退化为纯渐变）
  if (n === 0) {
    odata.set(sdata);
    return;
  }

  for (let y = 0; y < h; y++) {
    for (let x = 0; x < w; x++) {
      let sx = 0, sy = 0, sw = 0;

      for (let i = 0; i < n; i++) {
        const D = dstLines[i], S = srcLines[i];
        const dx = D.x2 - D.x1, dy = D.y2 - D.y1;
        const len2 = dx * dx + dy * dy;
        const len = Math.sqrt(len2) || 1e-6;
        const px = x - D.x1, py = y - D.y1;

        let u, v;
        if (len2 < 1e-9) {
          // 退化线（两端点重合）
          u = 0;
          v = Math.sqrt(px * px + py * py);
        } else {
          u = (px * dx + py * dy) / len2;
          v = (px * (-dy) + py * dx) / len;
        }

        // 在源线坐标系中的对应点 X'
        const dx2 = S.x2 - S.x1, dy2 = S.y2 - S.y1;
        const lenb = Math.sqrt(dx2 * dx2 + dy2 * dy2) || 1e-6;
        const sx_ = S.x1 + u * dx2 + v * (-dy2) / lenb;
        const sy_ = S.y1 + u * dy2 + v * (dx2) / lenb;

        // 点到线的距离：线段内用垂直距离，超出端点用端点距离
        let dist;
        if (u < 0) {
          dist = Math.sqrt(px * px + py * py);
        } else if (u > 1) {
          const qx = x - D.x2, qy = y - D.y2;
          dist = Math.sqrt(qx * qx + qy * qy);
        } else {
          dist = Math.abs(v);
        }

        // 影响权重：w = (len^p / (a + dist))^b
        const weight = Math.pow(Math.pow(len, P) / (A + dist), B);
        sx += sx_ * weight;
        sy += sy_ * weight;
        sw += weight;
      }

      if (sw > 1e-9) {
        sx /= sw;
        sy /= sw;
      }

      // 双线性采样（边缘 clamp）
      const x0 = Math.floor(sx), y0 = Math.floor(sy);
      const fx = sx - x0, fy = sy - y0;
      const xa = Math.max(0, Math.min(w - 1, x0));
      const xb = Math.max(0, Math.min(w - 1, x0 + 1));
      const ya = Math.max(0, Math.min(h - 1, y0));
      const yb = Math.max(0, Math.min(h - 1, y0 + 1));
      const ia = (ya * w + xa) * 4, ib = (ya * w + xb) * 4;
      const ic = (yb * w + xa) * 4, id2 = (yb * w + xb) * 4;
      const oi = (y * w + x) * 4;

      for (let c = 0; c < 3; c++) {
        const top = sdata[ia + c] + (sdata[ib + c] - sdata[ia + c]) * fx;
        const bot = sdata[ic + c] + (sdata[id2 + c] - sdata[ic + c]) * fx;
        odata[oi + c] = top + (bot - top) * fy;
      }
      odata[oi + 3] = 255;
    }
  }
}

function ensureScratch(w, h) {
  if (!scratch1 || scratch1.width !== w || scratch1.height !== h) {
    scratch1 = new ImageData(w, h);
    scratch2 = new ImageData(w, h);
    scratchRes = new ImageData(w, h);
  }
}

function renderMorph(t) {
  if (!state.src || !state.tgt) {
    drawPlaceholder(resCanvas, resCtx, '请加载图片');
    return;
  }

  const w = state.src.w, h = state.src.h;
  ensureScratch(w, h);

  // 中间帧的插值特征线
  const mid = state.lines.map(l => ({
    x1: l.src.x1 + (l.tgt.x1 - l.src.x1) * t,
    y1: l.src.y1 + (l.tgt.y1 - l.src.y1) * t,
    x2: l.src.x2 + (l.tgt.x2 - l.src.x2) * t,
    y2: l.src.y2 + (l.tgt.y2 - l.src.y2) * t,
  }));
  const srcLines = state.lines.map(l => l.src);
  const tgtLines = state.lines.map(l => l.tgt);

  // 分别 warp 源图与目标图到中间线
  warp(state.src, srcLines, mid, scratch1);
  warp(state.tgt, tgtLines, mid, scratch2);

  // 交叉融合：result = (1-t)*src + t*tgt
  const d1 = scratch1.data, d2 = scratch2.data, dr = scratchRes.data;
  const it = 1 - t;
  for (let i = 0; i < dr.length; i += 4) {
    dr[i]     = d1[i] * it     + d2[i] * t;
    dr[i + 1] = d1[i + 1] * it + d2[i + 1] * t;
    dr[i + 2] = d1[i + 2] * it + d2[i + 2] * t;
    dr[i + 3] = 255;
  }

  // 绘制到结果画布（等比缩放显示）
  fitCanvas(resCanvas, w, h);
  const disp = resCanvas._disp;
  resCtx.clearRect(0, 0, resCanvas.width, resCanvas.height);
  const tmp = document.createElement('canvas');
  tmp.width = w;
  tmp.height = h;
  tmp.getContext('2d').putImageData(scratchRes, 0, 0);
  resCtx.imageSmoothingEnabled = true;
  resCtx.drawImage(tmp, disp.ox, disp.oy, disp.dw, disp.dh);
}

/* ============================================================
 * 画布显示辅助
 * ============================================================ */

// 将画布尺寸设为父容器大小，并计算图像等比居中显示的变换
function fitCanvas(canvas, iw, ih) {
  const parent = canvas.parentElement;
  const cw = parent.clientWidth, ch = parent.clientHeight;
  if (cw < 2 || ch < 2) return;
  if (canvas.width !== cw || canvas.height !== ch) {
    canvas.width = cw;
    canvas.height = ch;
  }
  const scale = Math.min(cw / iw, ch / ih);
  canvas._disp = {
    ox: (cw - iw * scale) / 2,
    oy: (ch - ih * scale) / 2,
    dw: iw * scale,
    dh: ih * scale,
    scale: scale,
  };
}

function drawPlaceholder(canvas, ctx, text) {
  fitCanvas(canvas, 4, 3);
  if (!canvas._disp) return;
  ctx.fillStyle = '#fafafa';
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = '#888';
  ctx.font = '15px sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText(text, canvas.width / 2, canvas.height / 2);
  ctx.textAlign = 'start';
  ctx.textBaseline = 'alphabetic';
}

// 绘制一侧图像及其特征线覆盖层
function drawSide(canvas, ctx, img, side) {
  if (!img) {
    drawPlaceholder(canvas, ctx, side === 'src' ? '请加载源图' : '请加载目标图');
    return;
  }
  fitCanvas(canvas, img.w, img.h);
  if (!canvas._disp) return;
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(img.canvas, canvas._disp.ox, canvas._disp.oy, canvas._disp.dw, canvas._disp.dh);

  const disp = canvas._disp;
  const s = disp.dw / img.w;
  const sel = state.selectedId;

  for (const line of state.lines) {
    const p = line[side];
    const ax = disp.ox + p.x1 * s, ay = disp.oy + p.y1 * s;
    const bx = disp.ox + p.x2 * s, by = disp.oy + p.y2 * s;
    const isSel = line.id === sel;
    const baseColor = side === 'src' ? '#e74c3c' : '#2980b9';

    // 线体
    ctx.strokeStyle = isSel ? '#e6a700' : baseColor;
    ctx.lineWidth = isSel ? 3 : 2;
    ctx.beginPath();
    ctx.moveTo(ax, ay);
    ctx.lineTo(bx, by);
    ctx.stroke();

    // 端点
    for (const [ex, ey] of [[ax, ay], [bx, by]]) {
      ctx.fillStyle = isSel ? '#e6a700' : '#fff';
      ctx.strokeStyle = isSel ? '#e6a700' : baseColor;
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(ex, ey, isSel ? 6 : 5, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    }

    // 编号
    ctx.fillStyle = isSel ? '#b8860b' : '#666';
    ctx.font = '12px sans-serif';
    ctx.fillText(String(line.id), (ax + bx) / 2 + 7, (ay + by) / 2 - 7);
  }
}

// 整体刷新（重绘两侧 + 结果）
function refresh() {
  drawSide(srcCanvas, srcCtx, state.src, 'src');
  drawSide(tgtCanvas, tgtCtx, state.tgt, 'tgt');
  renderMorph(state.t);
}

let renderPending = false;
function requestRender() {
  if (renderPending) return;
  renderPending = true;
  requestAnimationFrame(() => {
    renderPending = false;
    refresh();
  });
}

/* ============================================================
 * 鼠标交互（联动绘制特征线）
 * ============================================================ */

function canvasPos(canvas, e) {
  const r = canvas.getBoundingClientRect();
  return {
    x: (e.clientX - r.left) * (canvas.width / r.width),
    y: (e.clientY - r.top) * (canvas.height / r.height),
  };
}

function canvasToImg(canvas, x, y) {
  const disp = canvas._disp;
  if (!disp) return null;
  return { x: (x - disp.ox) / disp.scale, y: (y - disp.oy) / disp.scale };
}

function clampImg(p, w, h) {
  return {
    x: Math.max(0, Math.min(w - 1, p.x)),
    y: Math.max(0, Math.min(h - 1, p.y)),
  };
}

function segDist(px, py, ax, ay, bx, by) {
  const dx = bx - ax, dy = by - ay;
  const len2 = dx * dx + dy * dy;
  let u;
  if (len2 < 1e-9) u = 0;
  else u = Math.max(0, Math.min(1, ((px - ax) * dx + (py - ay) * dy) / len2));
  const qx = ax + u * dx, qy = ay + u * dy;
  return Math.hypot(px - qx, py - qy);
}

function hitEndpoint(canvas, side, px, py) {
  const disp = canvas._disp;
  if (!disp) return null;
  const s = disp.dw / (state[side] ? state[side].w : 1);
  for (const line of state.lines) {
    const p = line[side];
    const ax = disp.ox + p.x1 * s, ay = disp.oy + p.y1 * s;
    const bx = disp.ox + p.x2 * s, by = disp.oy + p.y2 * s;
    if (Math.hypot(px - ax, py - ay) <= HIT_RADIUS) return { id: line.id, end: 0 };
    if (Math.hypot(px - bx, py - by) <= HIT_RADIUS) return { id: line.id, end: 1 };
  }
  return null;
}

function hitBody(canvas, side, px, py) {
  const disp = canvas._disp;
  if (!disp) return null;
  const s = disp.dw / (state[side] ? state[side].w : 1);
  let best = null, bestD = HIT_RADIUS;
  for (const line of state.lines) {
    const p = line[side];
    const ax = disp.ox + p.x1 * s, ay = disp.oy + p.y1 * s;
    const bx = disp.ox + p.x2 * s, by = disp.oy + p.y2 * s;
    const d = segDist(px, py, ax, ay, bx, by);
    if (d < bestD) { bestD = d; best = line.id; }
  }
  return best;
}

function onDown(side, e) {
  if (!state[side]) return;
  const canvas = side === 'src' ? srcCanvas : tgtCanvas;
  const { x, y } = canvasPos(canvas, e);
  const imgW = state[side].w, imgH = state[side].h;

  const hit = hitEndpoint(canvas, side, x, y);
  if (hit) {
    state.selectedId = hit.id;
    state.drag = { id: hit.id, side: side, end: hit.end };
    updateLineUI();
    requestRender();
    return;
  }

  const body = hitBody(canvas, side, x, y);
  if (body !== null) {
    state.selectedId = body;
    updateLineUI();
    requestRender();
    return;
  }

  // 新建成对特征线
  const p = canvasToImg(canvas, x, y);
  if (!p) return;
  const id = state.nextId++;
  const pt = { x1: p.x, y1: p.y, x2: p.x, y2: p.y };
  const line = { id: id, src: { ...pt }, tgt: { ...pt } };
  state.lines.push(line);
  state.selectedId = id;
  state.creating = { id: id, side: side, from: p };
  updateLineUI();
  requestRender();
}

function onMove(side, e) {
  const canvas = side === 'src' ? srcCanvas : tgtCanvas;
  const { x, y } = canvasPos(canvas, e);

  if (state.drag && state.drag.side === side) {
    const line = state.lines.find(l => l.id === state.drag.id);
    if (line && state[side]) {
      const p = clampImg(canvasToImg(canvas, x, y), state[side].w, state[side].h);
      const pt = line[side];
      if (state.drag.end === 0) { pt.x1 = p.x; pt.y1 = p.y; }
      else { pt.x2 = p.x; pt.y2 = p.y; }
      requestRender();
    }
  } else if (state.creating && state.creating.side === side) {
    const line = state.lines.find(l => l.id === state.creating.id);
    if (line && state[side]) {
      const p = clampImg(canvasToImg(canvas, x, y), state[side].w, state[side].h);
      const pt = line[side];
      pt.x2 = p.x;
      pt.y2 = p.y;
      requestRender();
    }
  }
}

function onUp() {
  if (state.creating) {
    const line = state.lines.find(l => l.id === state.creating.id);
    const p = line[state.creating.side];
    const len = Math.hypot(p.x2 - p.x1, p.y2 - p.y1);
    if (len < MIN_LINE_LEN) {
      state.lines = state.lines.filter(l => l.id !== state.creating.id);
      state.selectedId = null;
    }
    state.creating = null;
    updateLineUI();
  }
  state.drag = null;
  requestRender();
}

/* ============================================================
 * 控制交互
 * ============================================================ */

function updateLineUI() {
  lineCountEl.textContent = '特征线：' + state.lines.length;
  delBtn.disabled = state.selectedId == null;
  clearBtn.disabled = state.lines.length === 0;
}

function deleteSelected() {
  if (state.selectedId == null) return;
  state.lines = state.lines.filter(l => l.id !== state.selectedId);
  state.selectedId = null;
  state.drag = null;
  state.creating = null;
  updateLineUI();
  requestRender();
}

function clearLines() {
  state.lines = [];
  state.selectedId = null;
  state.drag = null;
  state.creating = null;
  updateLineUI();
  requestRender();
}

function togglePlay() {
  state.playing = !state.playing;
  playBtn.textContent = state.playing ? '暂停' : '播放';
  if (state.playing) {
    state.lastFrame = performance.now();
    requestAnimationFrame(playLoop);
  }
}

function playLoop(now) {
  if (!state.playing) return;
  const dt = Math.min(0.05, (now - state.lastFrame) / 1000);
  state.lastFrame = now;
  const speed = parseFloat(speedSel.value);
  state.t += (state.playDir * dt * speed) / LOOP_SECONDS;
  if (state.t >= 1) { state.t = 1; state.playDir = -1; }
  if (state.t <= 0) { state.t = 0; state.playDir = 1; }
  tSlider.value = Math.round(state.t * 100);
  tLabel.textContent = state.t.toFixed(2);
  requestRender();
  requestAnimationFrame(playLoop);
}

function onKey(e) {
  const tag = e.target && e.target.tagName;
  if (tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA') return;
  if (e.key === 'Delete' || e.key === 'Backspace') {
    e.preventDefault();
    deleteSelected();
  } else if (e.key === ' ') {
    e.preventDefault();
    togglePlay();
  }
}

/* ============================================================
 * 内置示例图片（Canvas 程序化生成）
 * ============================================================ */

function makeImage(w, h, draw) {
  const c = document.createElement('canvas');
  c.width = w;
  c.height = h;
  const ctx = c.getContext('2d');
  draw(ctx);
  return wrapCanvas(c);
}

function starPoints(cx, cy, outer, inner, n, rot) {
  const pts = [];
  for (let i = 0; i < n * 2; i++) {
    const r = i % 2 === 0 ? outer : inner;
    const a = rot + (i * Math.PI) / n;
    pts.push({ x: cx + Math.cos(a) * r, y: cy + Math.sin(a) * r });
  }
  return pts;
}

// 示例 1：圆形 → 五角星
function genSample1() {
  const N = 400, cx = 200, cy = 200;
  const fillShape = (ctx, pathFn) => {
    ctx.fillStyle = '#eef4ff';
    ctx.fillRect(0, 0, N, N);
    const g = ctx.createRadialGradient(cx - 30, cy - 30, 20, cx, cy, 160);
    g.addColorStop(0, '#ffd166');
    g.addColorStop(1, '#f4a261');
    ctx.fillStyle = g;
    ctx.beginPath();
    pathFn(ctx);
    ctx.closePath();
    ctx.fill();
    ctx.strokeStyle = '#d97706';
    ctx.lineWidth = 4;
    ctx.stroke();
  };

  const src = makeImage(N, N, ctx => {
    fillShape(ctx, c => c.arc(cx, cy, 150, 0, Math.PI * 2));
  });

  const star = starPoints(cx, cy, 150, 62, 5, -Math.PI / 2);
  const tgt = makeImage(N, N, ctx => {
    fillShape(ctx, c => {
      star.forEach((p, i) => (i === 0 ? c.moveTo(p.x, p.y) : c.lineTo(p.x, p.y)));
    });
  });

  const lines = [];
  for (let i = 0; i < 5; i++) {
    const tip = star[i * 2], val = star[i * 2 + 1];
    lines.push(
      { src: { x1: cx, y1: cy, x2: tip.x, y2: tip.y }, tgt: { x1: cx, y1: cy, x2: tip.x, y2: tip.y } },
      { src: { x1: cx, y1: cy, x2: val.x, y2: val.y }, tgt: { x1: cx, y1: cy, x2: val.x, y2: val.y } }
    );
  }
  return { src, tgt, lines };
}

// 示例 2：笑脸 → 严肃脸
function drawFace(ctx, N, smile) {
  ctx.fillStyle = '#eef2f7';
  ctx.fillRect(0, 0, N, N);
  // 头发
  ctx.fillStyle = '#6b4226';
  ctx.beginPath();
  ctx.arc(200, 160, 146, Math.PI * 1.05, Math.PI * 1.95);
  ctx.closePath();
  ctx.fill();
  // 脸
  ctx.fillStyle = '#f8c99b';
  ctx.beginPath();
  ctx.arc(200, 205, 140, 0, Math.PI * 2);
  ctx.fill();
  ctx.strokeStyle = 'rgba(0,0,0,0.15)';
  ctx.lineWidth = 3;
  ctx.stroke();
  // 眼睛
  ctx.fillStyle = '#2b2b2b';
  ctx.beginPath(); ctx.arc(150, 185, 15, 0, Math.PI * 2); ctx.fill();
  ctx.beginPath(); ctx.arc(250, 185, 15, 0, Math.PI * 2); ctx.fill();
  ctx.fillStyle = '#fff';
  ctx.beginPath(); ctx.arc(154, 180, 5, 0, Math.PI * 2); ctx.fill();
  ctx.beginPath(); ctx.arc(254, 180, 5, 0, Math.PI * 2); ctx.fill();
  // 嘴：微笑（二次曲线）或直线
  ctx.strokeStyle = '#8c4a3a';
  ctx.lineWidth = 6;
  ctx.lineCap = 'round';
  ctx.beginPath();
  if (smile) {
    ctx.moveTo(150, 250);
    ctx.quadraticCurveTo(200, 305, 250, 250);
  } else {
    ctx.moveTo(150, 250);
    ctx.lineTo(250, 250);
  }
  ctx.stroke();
}

function genSample2() {
  const N = 400;
  const src = makeImage(N, N, ctx => drawFace(ctx, N, true));
  const tgt = makeImage(N, N, ctx => drawFace(ctx, N, false));
  const mk = (x1, y1, x2, y2) => ({ src: { x1, y1, x2, y2 }, tgt: { x1, y1, x2, y2 } });
  const lines = [
    mk(135, 185, 165, 185), // 左眼
    mk(235, 185, 265, 185), // 右眼
    mk(150, 250, 250, 250), // 嘴
  ];
  return { src, tgt, lines };
}

// 示例 3：圆形 → 三角形
function genSample3() {
  const N = 400, cx = 200, cy = 200;
  const verts = [{ x: 200, y: 70 }, { x: 80, y: 320 }, { x: 320, y: 320 }];
  const mids = [
    { x: (verts[0].x + verts[1].x) / 2, y: (verts[0].y + verts[1].y) / 2 },
    { x: (verts[1].x + verts[2].x) / 2, y: (verts[1].y + verts[2].y) / 2 },
    { x: (verts[2].x + verts[0].x) / 2, y: (verts[2].y + verts[0].y) / 2 },
  ];
  const fill = (ctx, pathFn) => {
    ctx.fillStyle = '#eef4ff';
    ctx.fillRect(0, 0, N, N);
    const g = ctx.createRadialGradient(cx - 30, cy - 30, 20, cx, cy, 160);
    g.addColorStop(0, '#7ec8e3');
    g.addColorStop(1, '#2d6a9f');
    ctx.fillStyle = g;
    ctx.beginPath();
    pathFn(ctx);
    ctx.closePath();
    ctx.fill();
    ctx.strokeStyle = '#1d4e79';
    ctx.lineWidth = 4;
    ctx.stroke();
  };

  const src = makeImage(N, N, ctx => {
    fill(ctx, c => c.arc(cx, cy, 150, 0, Math.PI * 2));
  });

  const tgt = makeImage(N, N, ctx => {
    fill(ctx, c => {
      verts.forEach((p, i) => (i === 0 ? c.moveTo(p.x, p.y) : c.lineTo(p.x, p.y)));
    });
  });

  // 径向线：源图端点落在圆上同一角度，目标图端点为三角形顶点/边中点
  const lines = [...verts, ...mids].map(p => {
    const a = Math.atan2(p.y - cy, p.x - cx);
    const r = { x: cx + Math.cos(a) * 150, y: cy + Math.sin(a) * 150 };
    return {
      src: { x1: cx, y1: cy, x2: r.x, y2: r.y },
      tgt: { x1: cx, y1: cy, x2: p.x, y2: p.y },
    };
  });
  return { src, tgt, lines };
}

function assignIds(lines) {
  return lines.map(l => {
    const id = state.nextId++;
    return { id: id, src: { ...l.src }, tgt: { ...l.tgt } };
  });
}

function loadSample(n) {
  const s = genSample(n);
  state.src = s.src;
  state.tgt = s.tgt;
  state.lines = assignIds(s.lines);
  state.selectedId = null;
  state.drag = null;
  state.creating = null;
  state.t = 0.5;
  state.playing = false;
  state.playDir = 1;
  playBtn.textContent = '播放';
  tSlider.value = 50;
  tLabel.textContent = '0.50';
  updateLineUI();
  requestRender();
}

function genSample(n) {
  switch (n) {
    case 1: return genSample1();
    case 2: return genSample2();
    case 3: return genSample3();
    default: return genSample1();
  }
}

/* ============================================================
 * 初始化
 * ============================================================ */

function init() {
  tSlider.addEventListener('input', () => {
    state.t = tSlider.value / 100;
    tLabel.textContent = state.t.toFixed(2);
    requestRender();
  });
  playBtn.addEventListener('click', togglePlay);
  delBtn.addEventListener('click', deleteSelected);
  clearBtn.addEventListener('click', clearLines);
  srcFile.addEventListener('change', () => loadImageFromSource('src'));
  tgtFile.addEventListener('change', () => loadImageFromSource('tgt'));
  window.addEventListener('resize', requestRender);
  document.addEventListener('keydown', onKey);

  for (let i = 1; i <= 3; i++) {
    document.getElementById('sample' + i).addEventListener('click', () => loadSample(i));
  }

  srcCanvas.addEventListener('mousedown', e => onDown('src', e));
  tgtCanvas.addEventListener('mousedown', e => onDown('tgt', e));
  srcCanvas.addEventListener('mousemove', e => onMove('src', e));
  tgtCanvas.addEventListener('mousemove', e => onMove('tgt', e));
  window.addEventListener('mouseup', onUp);

  loadSample(1);
}

init();
