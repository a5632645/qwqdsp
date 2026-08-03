// ============================================================
//  piwarp5 复刻 — 纯 JS 实现
//  Kaiser 窗 + 帧内时间反转(×2×窗) + 重叠相加
//  与 C++ piwarp5 (M=32, hop=31) 的 DSP 保持一致
// ============================================================

'use strict';

const FS = 48000;   // 采样率
const M = 32;       // 子带数 (C++ filter_banks)
const HOP = M - 1;  // hop = 31

// ------------------------------------------------------------
//  状态
// ------------------------------------------------------------
const state = {
    // 第 0 列: null = 原始 (不处理); 之后为 { lenMult, beta }
    params: [
        null,
        { lenMult: 1, beta: 2 },
        { lenMult: 2, beta: 8 },
        { lenMult: 3, beta: 12 },
        { lenMult: 4, beta: 32 },
    ],
    whiteNoise: null,
    sweep: null,
};

// ------------------------------------------------------------
//  DSP 基础
// ------------------------------------------------------------

/** 复刻 C++ WhiteNoise (LCG, uint32 运算, [0,1] → [-1,1]) */
class WhiteNoise {
    constructor(seed) {
        this.reg = seed >>> 0;
    }
    next() {
        // reg_ = reg_ * 1103515245 + 12345 (uint32 回绕)
        this.reg = (Math.imul(this.reg, 1103515245) + 12345) >>> 0;
        return (this.reg / 4294967295) * 2 - 1; // [-1,1]
    }
}

/** 修正 Bessel I0 — 级数展开 */
function besselI0(x) {
    let sum = 1, term = 1;
    for (let k = 1; k < 100; k++) {
        term *= (x / (2 * k)) * (x / (2 * k));
        sum += term;
        if (term < 1e-12 * sum) break;
    }
    return sum;
}

/** Kaiser 窗 (对称, 与 C++ Kaiser::Window(win, beta, false) 一致) */
function kaiser(N, beta) {
    const w = new Float32Array(N);
    const i0b = besselI0(beta);
    for (let i = 0; i < N; i++) {
        const t = (2 * i) / (N - 1) - 1;
        w[i] = besselI0(beta * Math.sqrt(1 - t * t)) / i0b;
    }
    return w;
}

/** 指数扫频 f0 → f1 (对数频率线性变化) */
function expChirp(N, f0, f1) {
    const x = new Float32Array(N);
    const T = N / FS;
    const k = Math.log(f1 / f0);
    for (let i = 0; i < N; i++) {
        const t = i / FS;
        x[i] = Math.sin((2 * Math.PI * f0 * T * (Math.pow(f1 / f0, t / T) - 1)) / k);
    }
    return x;
}

/** 迭代基-2 FFT (原地, 复数 re/im) */
function fft(re, im) {
    const n = re.length;
    // 位反转重排
    for (let i = 1, j = 0; i < n; i++) {
        let bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            [re[i], re[j]] = [re[j], re[i]];
            [im[i], im[j]] = [im[j], im[i]];
        }
    }
    // 蝶形运算
    for (let len = 2; len <= n; len <<= 1) {
        const ang = (-2 * Math.PI) / len;
        const wRe = Math.cos(ang), wIm = Math.sin(ang);
        for (let i = 0; i < n; i += len) {
            let wr = 1, wi = 0;
            const half = len >> 1;
            for (let j = 0; j < half; j++) {
                const eRe = re[i + j], eIm = im[i + j];
                const oRe = re[i + j + half], oIm = im[i + j + half];
                const tr = wr * oRe - wi * oIm;
                const ti = wr * oIm + wi * oRe;
                re[i + j + half] = eRe - tr;
                im[i + j + half] = eIm - ti;
                re[i + j] = eRe + tr;
                im[i + j] = eIm + ti;
                const nwr = wr * wRe - wi * wIm;
                const nwi = wr * wIm + wi * wRe;
                wr = nwr;
                wi = nwi;
            }
        }
    }
}

/**
 * piwarp5: 帧内时间反转 ×2 × Kaiser 窗 + 重叠相加
 * out[m*hop + k] += 2 * win[k] * x[m*hop + (len-1-k)]
 * 归一化: / (len/hop), 再峰值归一化到 0.9 (与 C++ 端一致)
 */
function piwarp5(x, len, beta) {
    if (len > x.length) return x.slice(); // 保护
    const win = kaiser(len, beta);
    const N = x.length;
    const numFrames = Math.floor((N - len) / HOP) + 1;
    const outLen = (numFrames - 1) * HOP + len;
    const out = new Float32Array(outLen);
    for (let m = 0; m < numFrames; m++) {
        const start = m * HOP;
        for (let k = 0; k < len; k++) {
            out[start + k] += 2 * win[k] * x[start + (len - 1 - k)];
        }
    }
    // AnalyzeSynthsisOffline 归一化: / (len / hop)
    const g = HOP / len;
    for (let i = 0; i < outLen; i++) out[i] *= g;
    // 峰值归一化到 0.9, 避免 clip
    let peak = 0;
    for (let i = 0; i < outLen; i++) peak = Math.max(peak, Math.abs(out[i]));
    if (peak > 1e-9) {
        const gg = 0.9 / peak;
        for (let i = 0; i < outLen; i++) out[i] *= gg;
    }
    return out;
}

/** STFT 幅度谱 (dB) — 时频图数据 */
function spectrogram(x, nfft = 2048, hop = 512) {
    // Hann 窗 (对称)
    const win = new Float64Array(nfft);
    for (let i = 0; i < nfft; i++) {
        win[i] = 0.5 - 0.5 * Math.cos((2 * Math.PI * i) / (nfft - 1));
    }
    const re = new Float64Array(nfft);
    const im = new Float64Array(nfft);
    const bins = nfft / 2 + 1;
    const nframes = Math.floor((x.length - nfft) / hop) + 1;
    const spec = new Float32Array(nframes * bins);
    for (let m = 0; m < nframes; m++) {
        const off = m * hop;
        for (let i = 0; i < nfft; i++) {
            re[i] = x[off + i] * win[i];
            im[i] = 0;
        }
        fft(re, im);
        const base = m * bins;
        for (let k = 0; k < bins; k++) {
            spec[base + k] = 10 * Math.log10(re[k] * re[k] + im[k] * im[k] + 1e-12);
        }
    }
    return { spec, nframes, bins };
}

// ------------------------------------------------------------
//  渲染
// ------------------------------------------------------------

/** 波形 (min/max 包络) */
function drawWaveform(canvas, x) {
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height;
    ctx.clearRect(0, 0, W, H);
    // 中线
    ctx.strokeStyle = 'rgba(255,255,255,0.15)';
    ctx.beginPath();
    ctx.moveTo(0, H / 2);
    ctx.lineTo(W, H / 2);
    ctx.stroke();
    // min/max 包络
    const N = x.length;
    ctx.fillStyle = '#4fc3f7';
    for (let px = 0; px < W; px++) {
        const i0 = Math.floor((px * N) / W);
        const i1 = Math.max(i0 + 1, Math.floor(((px + 1) * N) / W));
        let mn = 1, mx = -1;
        for (let i = i0; i < i1; i++) {
            const v = x[i];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        const y0 = H / 2 - mx * (H / 2 - 2);
        const y1 = H / 2 - mn * (H / 2 - 2);
        ctx.fillRect(px, Math.min(y0, y1), 1, Math.max(1, Math.abs(y1 - y0)));
    }
}

/** 类 viridis 色图 */
const CMAP = [
    [68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37],
];
function colormap(t) {
    const x = Math.max(0, Math.min(1, t)) * (CMAP.length - 1);
    const i = Math.floor(x), f = x - i;
    const a = CMAP[i], b = CMAP[Math.min(i + 1, CMAP.length - 1)];
    return [
        Math.round(a[0] + (b[0] - a[0]) * f),
        Math.round(a[1] + (b[1] - a[1]) * f),
        Math.round(a[2] + (b[2] - a[2]) * f),
    ];
}

/** 时频图 (顶部=高频, 线性频率轴) */
function drawSpectrogram(canvas, spec, nframes, bins, minDb = -80, maxDb = 0) {
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height;
    const img = ctx.createImageData(W, H);
    const d = img.data;
    for (let y = 0; y < H; y++) {
        const bin = Math.round((bins - 1) * (1 - y / (H - 1)));
        const row = y * W * 4;
        for (let px = 0; px < W; px++) {
            const frame = Math.floor((px * nframes) / W);
            const v = spec[frame * bins + bin];
            const t = (v - minDb) / (maxDb - minDb);
            const c = colormap(t);
            const o = row + px * 4;
            d[o] = c[0];
            d[o + 1] = c[1];
            d[o + 2] = c[2];
            d[o + 3] = 255;
        }
    }
    ctx.putImageData(img, 0, 0);
}

// ------------------------------------------------------------
//  网格构建
// ------------------------------------------------------------
const paramsCol = document.getElementById('params-col');
const noiseCol = document.getElementById('noise-col');
const sweepCol = document.getElementById('sweep-col');
const seedInput = document.getElementById('seed-input');
const f0Input = document.getElementById('f0-input');
const f1Input = document.getElementById('f1-input');
const durInput = document.getElementById('dur-input');
const regenerateBtn = document.getElementById('regenerate-btn');
const statusEl = document.getElementById('status');

function buildGrid() {
    paramsCol.innerHTML = '';
    noiseCol.innerHTML = '';
    sweepCol.innerHTML = '';

    for (let i = 0; i < state.params.length; i++) {
        // --- 参数卡片 ---
        const pc = document.createElement('div');
        pc.className = 'cell param-cell';
        if (i === 0) {
            pc.innerHTML = '<div class="cell-head">#0 原始</div>'
                + '<div class="param-note">输入信号, 不处理</div>';
        }
        else {
            const p = state.params[i];
            pc.innerHTML = `<div class="cell-head">#${i} 组合</div>`
                + `<label>len_mult <input type="number" class="len-mult" value="${p.lenMult}" min="1" step="1"></label>`
                + `<label>beta <input type="number" class="beta" value="${p.beta}" min="0" step="0.5"></label>`
                + '<button class="del-btn" title="删除">×</button>';
            pc.querySelector('.len-mult').addEventListener('change', (e) => {
                p.lenMult = Math.max(1, parseInt(e.target.value, 10) || 1);
                renderColumn(i);
            });
            pc.querySelector('.beta').addEventListener('change', (e) => {
                p.beta = Math.max(0, parseFloat(e.target.value) || 0);
                renderColumn(i);
            });
            pc.querySelector('.del-btn').addEventListener('click', () => {
                state.params.splice(i, 1);
                rebuild();
            });
        }
        paramsCol.appendChild(pc);

        // --- 白噪声输出单元格 ---
        const nc = document.createElement('div');
        nc.className = 'cell out-cell';
        nc.innerHTML = `<div class="cell-head">白噪声 #${i}</div>`
            + '<canvas class="wave"></canvas>'
            + '<canvas class="spec"></canvas>'
            + '<button class="play-btn">▶ 播放</button>';
        noiseCol.appendChild(nc);

        // --- 扫频输出单元格 ---
        const sc = document.createElement('div');
        sc.className = 'cell out-cell';
        sc.innerHTML = `<div class="cell-head">扫频 #${i}</div>`
            + '<canvas class="wave"></canvas>'
            + '<canvas class="spec"></canvas>'
            + '<button class="play-btn">▶ 播放</button>';
        sweepCol.appendChild(sc);
    }

    // --- 添加按钮 (参数列底部) ---
    const addBtn = document.createElement('button');
    addBtn.className = 'add-btn';
    addBtn.textContent = '+ 添加组合';
    addBtn.addEventListener('click', () => {
        state.params.push({ lenMult: 2, beta: 8 });
        rebuild();
    });
    paramsCol.appendChild(addBtn);
}

// ------------------------------------------------------------
//  渲染单列 / 全部
// ------------------------------------------------------------
/** 画布像素尺寸与 CSS 显示尺寸同步 (保证清晰) */
function fitCanvas(c) {
    const w = c.clientWidth, h = c.clientHeight;
    if (w > 0 && h > 0 && (c.width !== w || c.height !== h)) {
        c.width = w;
        c.height = h;
    }
}

/** 渲染单列 (一个参数索引) */
async function renderColumn(i) {
    const p = state.params[i];
    const len = p ? M * p.lenMult : 0;
    const beta = p ? p.beta : 0;

    // 输出: 第 0 列为原始输入
    const noiseOut = i === 0 ? state.whiteNoise : piwarp5(state.whiteNoise, len, beta);
    const sweepOut = i === 0 ? state.sweep : piwarp5(state.sweep, len, beta);

    const nc = noiseCol.children[i];
    const sc = sweepCol.children[i];

    // 波形
    fitCanvas(nc.querySelector('.wave'));
    fitCanvas(sc.querySelector('.wave'));
    drawWaveform(nc.querySelector('.wave'), noiseOut);
    drawWaveform(sc.querySelector('.wave'), sweepOut);

    // 时频图 (计算量大, 逐索引让出主线程)
    fitCanvas(nc.querySelector('.spec'));
    const specN = spectrogram(noiseOut);
    drawSpectrogram(nc.querySelector('.spec'), specN.spec, specN.nframes, specN.bins);
    fitCanvas(sc.querySelector('.spec'));
    const specS = spectrogram(sweepOut);
    drawSpectrogram(sc.querySelector('.spec'), specS.spec, specS.nframes, specS.bins);

    // 播放
    nc.querySelector('.play-btn').onclick = () => playBuffer(noiseOut);
    sc.querySelector('.play-btn').onclick = () => playBuffer(sweepOut);

    await new Promise((r) => setTimeout(r, 0));
}

async function renderAll() {
    statusEl.textContent = '计算中…';
    for (let i = 0; i < state.params.length; i++) {
        await renderColumn(i);
    }
    statusEl.textContent = '';
}

function rebuild() {
    buildGrid();
    renderAll();
}

// ------------------------------------------------------------
//  播放 (Web Audio)
// ------------------------------------------------------------
let audioCtx = null;
let currentSrc = null;

function playBuffer(buffer) {
    if (!audioCtx) {
        audioCtx = new (window.AudioContext || window.webkitAudioContext)();
    }
    audioCtx.resume();
    if (currentSrc) {
        try { currentSrc.stop(); } catch (e) { /* ignore */ }
        currentSrc = null;
    }
    const src = audioCtx.createBufferSource();
    const ab = audioCtx.createBuffer(1, buffer.length, FS);
    ab.getChannelData(0).set(buffer);
    src.buffer = ab;
    src.connect(audioCtx.destination);
    src.start();
    currentSrc = src;
}

// ------------------------------------------------------------
//  输入生成 & 启动
// ------------------------------------------------------------
function regenerate() {
    const seed = (parseInt(seedInput.value, 10) || 0) >>> 0;
    const f0 = parseFloat(f0Input.value) || 20;
    const f1 = parseFloat(f1Input.value) || 20000;
    const dur = parseFloat(durInput.value) || 10;
    const N = Math.round(dur * FS);

    // 白噪声 (复刻 C++ WhiteNoise LCG)
    const wn = new WhiteNoise(seed);
    state.whiteNoise = new Float32Array(N);
    for (let i = 0; i < N; i++) state.whiteNoise[i] = wn.next();

    // 指数扫频
    state.sweep = expChirp(N, f0, f1);

    renderAll();
}

regenerateBtn.addEventListener('click', regenerate);
seedInput.addEventListener('change', regenerate);
f0Input.addEventListener('change', regenerate);
f1Input.addEventListener('change', regenerate);
durInput.addEventListener('change', regenerate);

buildGrid();
regenerate();
