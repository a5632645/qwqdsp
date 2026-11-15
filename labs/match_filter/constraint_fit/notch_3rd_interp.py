"""三阶 notch：结构约束 + BRASIL 区间等化。

模拟 notch 的功率（``r = w/wc``）

    Pa(w) = (1 - r^2)^2 / ((1 - r^2)^2 + (r/Q)^2)

硬约束：

    P(0) = 1,   P(x_fc) = 0

三阶数字滤波器的功率在 ``x = sin^2(w/2)`` 下是 3/3 有理函数：

    N(x) = (1 - x/xc)^2 (1 + c x)
    D(x) = 1 + d1 x + d2 x^2 + d3 x^3

缺口二重零点固定；``c`` 是额外分子因子，``c >= -1`` 保证分子非负。
共 4 个自由参数（c,d1,d2,d3）→ 4 个自由节点 / 5 段，BRASIL 等化
5 段 max|ΔdB|。缺口处的 0/0 取二阶曲率比极限。

本脚本只设计功率响应；``D > 0``、``N >= 0`` 保证存在稳定的三阶
谱分解。

结果（fs = 48 kHz，max|ΔdB|）：fc <= 4 kHz 为 0.0000–0.0136，
23.9 kHz 为 0.430–0.596；全表 0–0.60。相比二阶 notch（fc <= 4 kHz
0.103、23.9 kHz 1.47–1.70），高频档改善约 3 倍。

运行：``python qwqdsp/labs/match_filter/constraint_fit/notch_3rd_interp.py``
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "best_pole_fit"))
import brasil as BRASIL  # noqa: E402
import power_match as PM  # noqa: E402

PI = PM.PI
FS = PM.FS
FC_LIST = PM.FC_LIST
Q_LIST = (0.5, 1.0 / np.sqrt(2.0), 1.0, 2.0)


def analog_power(w, wc, Q):
    r = np.asarray(w, dtype=float) / wc
    a = (1.0 - r * r) ** 2
    return a / (a + (r / Q) ** 2)


class Notch3:
    """三阶 notch 的功率响应；不保存 b/a 相位分解。"""

    def __init__(self, fc, Q, c, d1, d2, d3, fs=FS):
        self.fc = fc
        self.Q = Q
        self.fs = fs
        self.c = c
        self.d = (d1, d2, d3)
        self.wc = 2.0 * PI * fc / fs
        self.xc = np.sin(self.wc / 2.0) ** 2

    def power(self, w):
        x = np.sin(np.asarray(w, dtype=float) / 2.0) ** 2
        d1, d2, d3 = self.d
        N = (1.0 - x / self.xc) ** 2 * (1.0 + self.c * x)
        D = 1.0 + d1 * x + d2 * x * x + d3 * x ** 3
        return N / D


def _match(fc, Q, nodes, fs=FS):
    wc = 2.0 * PI * fc / fs
    xc = np.sin(wc / 2.0) ** 2
    xs = np.sort(np.asarray(nodes, float))
    if (xs.size != 4 or np.min(np.diff(xs)) <= 1e-12
            or np.any(np.abs(xs - xc) < 1e-12)):
        return None, "nodes"
    w = 2.0 * np.arcsin(np.sqrt(xs))
    vals = analog_power(w, wc, Q)
    tmp = 1.0 - xs / xc
    N0 = tmp * tmp
    Nc = xs * N0
    M = np.stack([Nc, -vals * xs, -vals * xs * xs, -vals * xs ** 3], axis=1)
    scale = np.linalg.norm(M, axis=0)
    try:
        c, d1, d2, d3 = np.linalg.solve(M / scale, vals - N0) / scale
    except np.linalg.LinAlgError:
        return None, "singular"
    if c < -1.0:
        return None, "N < 0"
    xx = np.linspace(0.0, 1.0, 2001)
    if np.min(1.0 + d1 * xx + d2 * xx ** 2 + d3 * xx ** 3) <= 0.0:
        return None, "D <= 0"
    return Notch3(fc, Q, c, d1, d2, d3, fs), "ok"


def error_db(match, x, fc, Q, fs=FS):
    x = np.asarray(x, dtype=float)
    wc = 2.0 * PI * fc / fs
    xc = np.sin(wc / 2.0) ** 2
    d1, d2, d3 = match.d
    D = 1.0 + d1 * x + d2 * x * x + d3 * x ** 3
    w = 2.0 * np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0)))
    r = w / wc
    A = (1.0 - r * r) ** 2
    B = A + (r / Q) ** 2
    dx = x - xc
    tol = 1e-13 * max(xc, 1e-300)
    dwdx = np.empty_like(x)
    near = np.abs(dx) <= tol
    dwdx[near] = 1.0 / np.sqrt(xc * (1.0 - xc))
    dwdx[~near] = (w[~near] - wc) / dx[~near]
    a_over = ((dwdx * (1.0 + r)) / wc) ** 2
    n_over = (1.0 + match.c * x) / xc ** 2
    ratio = (n_over / a_over) * (B / D)
    return 10.0 * np.log10(np.maximum(ratio, 1e-300))


def _measure_grid(fc, fs=FS):
    xc = np.sin(PI * fc / fs) ** 2
    xs = np.linspace(0.0, 1.0, 4001)
    width = max(4.0 * xc * (1.0 - xc), 1e-12)
    loc = xc + np.linspace(-1.0, 1.0, 201) * min(width, 0.25)
    return np.unique(np.clip(np.concatenate((xs, loc, [xc])), 0.0, 1.0))


def design(fc, Q, fs=FS):
    xc = np.sin(PI * fc / fs) ** 2
    xgrid = _measure_grid(fc, fs)

    def measure(nodes):
        if np.any(np.abs(nodes - xc) < 1e-12):
            return None
        match, _ = _match(fc, Q, nodes, fs)
        if match is None:
            return None
        return BRASIL.segment_maxima(
            xgrid, np.abs(error_db(match, xgrid, fc, Q, fs)), nodes)

    near = max(2.0 * xc, 1e-8)
    sets = ((0.02, 0.2, 0.6, 0.95), (0.01, 0.1, 0.5, 0.9),
            (0.05, 0.3, 0.7, 0.99), (0.3, 0.6, 0.9, 0.999),
            (0.001, 0.05, 0.4, 0.9), (0.1, 0.5, 0.85, 0.99),
            (max(0.1 * xc, 1e-8), max(0.5 * xc, 1e-8), near, 0.7),
            (0.01, 0.3, 0.95, 0.999),
            (0.9, 0.999, 0.9999, 0.99999),
            (0.99, 0.9999, 0.99999, 0.999999))
    best = (np.inf, None)
    for t in sets:
        nodes0 = np.sort(np.clip(np.asarray(t, float), 1e-9, 1.0 - 1e-9))
        if (np.min(np.diff(nodes0)) <= 1e-12
                or np.any(np.abs(nodes0 - xc) < 1e-12)):
            continue
        nodes, _, err = BRASIL.equalize(nodes0, measure)
        if nodes is not None and err < best[0]:
            best = (err, nodes)
    if best[1] is None:
        return None, None, np.inf
    match, _ = _match(fc, Q, best[1], fs)
    if match is None:
        return None, None, np.inf
    return match, best[1], float(np.max(np.abs(
        error_db(match, xgrid, fc, Q, fs))))


def main():
    print(f"fs = {FS:g} Hz；三阶 notch：P(0)=1, P(fc)=0 + BRASIL")
    print(f"{'':>12}" + "".join(f"{fc:>12.1f}" for fc in FC_LIST))
    rows = {}
    for Q in Q_LIST:
        cells = []
        for fc in FC_LIST:
            match, _, err = design(fc, Q)
            rows[(Q, fc)] = match
            cells.append(f"{err:12.4f}")
        print(f"{f'Q={Q:g}':>12}" + "".join(cells))
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    _figure(rows, out)
    print(f"图：{out / 'notch_3rd_interp.png'}")


def _figure(rows, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(4, 2, figsize=(12.5, 14.0),
                             gridspec_kw={"wspace": 0.22, "hspace": 0.35})
    f = np.geomspace(10.0, FS / 2.0, 3001)
    w = 2.0 * PI * f / FS
    xs = np.sin(np.linspace(0.0, PI, 20001) / 2.0) ** 2
    ww = 2.0 * np.arcsin(np.sqrt(xs))
    for r, Q in enumerate(Q_LIST):
        ax_r, ax_e = axes[r]
        for j, fc in enumerate(FC_LIST):
            match = rows[(Q, fc)]
            if match is None:
                continue
            color = f"C{j}"
            ax_r.semilogx(
                f, 10.0 * np.log10(np.maximum(
                    analog_power(w, 2.0 * PI * fc / FS, Q), 1e-300)),
                color=color, lw=3.0, ls="--", alpha=0.45)
            ax_r.semilogx(
                f, 10.0 * np.log10(np.maximum(match.power(w), 1e-300)),
                color=color, lw=1.2, label=f"{fc:g} Hz")
            err = error_db(match, xs, fc, Q)
            ax_e.plot(ww / PI * FS / 2.0, np.maximum(np.abs(err), 1e-12),
                      color=color, lw=1.1)
        ax_r.set_xscale("log")
        ax_r.set_xlim(10.0, FS / 2.0)
        ax_r.set_ylim(-110.0, 25.0)
        ax_r.set_ylabel(f"Q = {Q:g}\n幅度 / dB")
        ax_r.grid(True, which="both", alpha=0.25)
        ax_e.set_xlim(0.0, FS / 2.0)
        ax_e.set_yscale("log")
        ax_e.set_ylim(1e-8, 1e1)
        ax_e.set_ylabel("|ΔdB|")
        ax_e.grid(True, which="both", alpha=0.25)
    axes[0][0].legend(fontsize=7, title="fc", ncol=2)
    axes[-1][0].set_xlabel("频率 / Hz")
    axes[-1][1].set_xlabel("频率 / Hz")
    fig.suptitle("三阶 notch: P(0)=1, P(fc)=0+BRASIL 4点有理插值", y=0.995)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.05)
    fig.savefig(out / "notch_3rd_interp.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
