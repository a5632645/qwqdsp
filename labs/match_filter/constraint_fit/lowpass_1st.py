"""一阶低通：硬约束法的负结果。

数字二阶节 5 个系数（分子 3、分母 2）。固定三条约束

    P(0) = 1,    P(fc) = 1/2,    P'(fc) = Pa'(fc),

再给分母两条额外约束，就得到 5 元非线性方程组。本文件枚举候选约束对。

**严格残差检查下，5 个 fc 上没有任何一对约束可解**：``root`` 返回的
``success=True`` 解残差在 1e-6 … 1e-2 量级，属于最小二乘折中而非精确解。
这说明“锁定分子 + 五条硬约束”本身不相容；正确做法是把匹配视为有理插值，
见 ``lowpass_1st_interp.py``。

运行：``python qwqdsp/labs/match_filter/constraint_fit/lowpass_1st.py``
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
from scipy.optimize import root

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "best_pole_fit"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "ivantsovy"))
import power_match as PM  # noqa: E402
import check_two_point as C  # noqa: E402
from two_point import lowpass as iv_lowpass  # noqa: E402

PI = PM.PI
"""圆周率。"""

FS = PM.FS
"""实验采样率。"""

FC_LIST = PM.FC_LIST
"""与其余实验一致的 5 个截止频率。"""

H = 1e-6
"""数值微分步长。"""


# --------------------------------------------------------------------- 模拟量
def analog_P(w, wc):
    r = np.asarray(w, dtype=float) / wc
    return 1.0 / (1.0 + r * r)


def analog_dP(w, wc):
    r = np.asarray(w, dtype=float) / wc
    return -2.0 * r / (wc * (1.0 + r * r) ** 2)


def analog_d2P(w, wc):
    r = np.asarray(w, dtype=float) / wc
    return 2.0 * (3.0 * r * r - 1.0) / (wc * wc * (1.0 + r * r) ** 3)


def analog_phase(w, wc):
    return -np.arctan(np.asarray(w, dtype=float) / wc)


def analog_delay(w, wc):
    r = np.asarray(w, dtype=float) / wc
    return 1.0 / (wc * (1.0 + r * r))


# --------------------------------------------------------------------- 数字量
def digital_P(B, A, w):
    ph = PM.phi(np.asarray(w, dtype=float))
    return (ph @ B) / (ph @ A)


def digital_dP(B, A, w):
    return (digital_P(B, A, w + H) - digital_P(B, A, w - H)) / (2.0 * H)


def digital_d2P(B, A, w):
    return (digital_P(B, A, w + H) - 2.0 * digital_P(B, A, w)
            + digital_P(B, A, w - H)) / (H * H)


def digital_H(b, a, w):
    z = np.exp(-1j * np.asarray(w, dtype=float))
    return (b[0] + b[1] * z + b[2] * z * z) / (1.0 + a[0] * z + a[1] * z * z)


def digital_delay(b, a, w):
    ph = np.unwrap(np.angle(digital_H(b, a, np.array([w - H, w, w + H]))))
    return -(ph[2] - ph[0]) / (2.0 * H)


# --------------------------------------------------------------------- 约束
_KEYS = ("P_pi", "dP_pi", "d2P_0", "d2P_pi", "d2P_wc", "phase_wc", "delay_0")


def digital_cond(key, B, A, b, a, wc):
    if key == "P_pi":
        return float(digital_P(B, A, np.pi))
    if key == "dP_pi":
        return float(digital_dP(B, A, np.pi))
    if key == "d2P_0":
        return float(digital_d2P(B, A, 0.0))
    if key == "d2P_pi":
        return float(digital_d2P(B, A, np.pi))
    if key == "d2P_wc":
        return float(digital_d2P(B, A, wc))
    if key == "phase_wc":
        return float(np.angle(digital_H(b, a, np.array([wc]))[0]))
    if key == "delay_0":
        return float(digital_delay(b, a, 0.0))
    raise KeyError(key)


def analog_cond(key, wc):
    if key == "P_pi":
        return float(analog_P(np.pi, wc))
    if key == "dP_pi":
        return float(analog_dP(np.pi, wc))
    if key == "d2P_0":
        return float(analog_d2P(0.0, wc))
    if key == "d2P_pi":
        return float(analog_d2P(np.pi, wc))
    if key == "d2P_wc":
        return float(analog_d2P(wc, wc))
    if key == "phase_wc":
        return float(analog_phase(wc, wc))
    if key == "delay_0":
        return float(analog_delay(0.0, wc))
    raise KeyError(key)


PAIRS = (
    ("P_pi", "dP_pi"),
    ("P_pi", "d2P_0"),
    ("P_pi", "d2P_pi"),
    ("dP_pi", "d2P_0"),
    ("d2P_0", "d2P_pi"),
    ("P_pi", "d2P_wc"),
    ("P_pi", "phase_wc"),
    ("P_pi", "delay_0"),
)
"""分母的两条约束（从 ``_KEYS`` 里选）。"""


# --------------------------------------------------------------------- 求解
def _a_of(p, q):
    return np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])


def _numerator3(wc, p, q):
    """0/fc/slope 三条对 B 线性，给定分母直接解出分子。"""
    A = PM.alpha_of(p, q)
    ph0 = PM.phi(np.array([0.0]))[0]
    phc = PM.phi(np.array([wc]))[0]
    sc = np.sin(wc / 2.0) ** 2
    dsc = np.sin(wc) / 2.0
    dphc = np.array([-1.0, 1.0, 4.0 * (1.0 - 2.0 * sc)]) * dsc
    ac = float(phc @ A)
    row = dphc * ac - phc * float(dphc @ A)
    rows = np.stack([ph0, phc, row])
    rhs = np.array([1.0 * A[0], 0.5 * ac, analog_dP(wc, wc) * ac * ac])
    return np.linalg.solve(rows, rhs)


def _make_match(B, p, q):
    A = PM.alpha_of(p, q)
    if not PM.feasible(B):
        return None
    return PM.Match(PM.solve_b(B), _a_of(p, q), float(p), float(q), B, A)


def error_db(match, fc, fs, n=20001):
    s = np.linspace(0.0, 1.0, n)
    ph = np.stack([1.0 - s, s, 4.0 * s * (1.0 - s)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    w = 2.0 * np.arcsin(np.sqrt(s))
    ha2 = analog_P(w, 2.0 * PI * fc / fs)
    return 10.0 * np.log10(np.maximum(hd2, 1e-300) / np.maximum(ha2, 1e-300)), s, w


def design(pair, fc, fs=FS):
    """固定 0/fc/slope + 额外两条约束，解 5 元非线性方程组。"""
    wc = 2.0 * PI * fc / fs
    if pair[0] == "ivantsov":
        _, a = C.expand(iv_lowpass(fc, fs))
        a1 = float(a[1])
        a2 = float(a[2]) if len(a) > 2 else 0.0
        p = 1.0 + a1 + a2
        q = 1.0 - a1 + a2
        B = _numerator3(wc, p, q)
        match = _make_match(B, p, q)
        if match is None:
            return None, np.inf
        err, _, _ = error_db(match, fc, fs, 400001)
        return match, float(np.max(np.abs(err)))
    r1 = PM.bilinear_pole(fc, fs)
    p0 = 1.0 - r1
    q0 = 1.0 + r1
    B0 = _numerator3(wc, p0, q0)
    x0 = np.concatenate([B0, [p0, q0]])

    def residual(x):
        B = x[:3]
        p, q = x[3], x[4]
        if not (p > 0.0 and q > 0.0 and p + q < 4.0):
            return np.full(5, 1e3)
        A = PM.alpha_of(p, q)
        b, a = PM.solve_b(B), _a_of(p, q)
        out = np.array([
            float(P_digital(B, A, 0.0) - 1.0),
            float(P_digital(B, A, wc) - 0.5),
            float(dP_digital(B, A, wc) - analog_dP(wc, wc)),
            float(digital_cond(pair[0], B, A, b, a, wc)
                  - analog_cond(pair[0], wc)),
            float(digital_cond(pair[1], B, A, b, a, wc)
                  - analog_cond(pair[1], wc)),
        ])
        return out

    res = root(residual, x0, method="hybr", options={"maxfev": 4000})
    if not res.success or np.max(np.abs(res.fun)) > 1e-8:
        return None, np.inf
    B, p, q = res.x[:3], float(res.x[3]), float(res.x[4])
    match = _make_match(B, p, q)
    if match is None:
        return None, np.inf
    err, _, _ = error_db(match, fc, fs, 400001)
    return match, float(np.max(np.abs(err)))


def P_digital(B, A, w):
    return float(digital_P(B, A, np.array([w]))[0])


def dP_digital(B, A, w):
    return float(digital_dP(B, A, np.array([w]))[0])


# --------------------------------------------------------------------- 实验
def _format_pole(z):
    if abs(z.imag) < 1e-9:
        return f"{z.real:.6f}"
    return f"{z.real:.6f}{z.imag:+.6f}j"


def main():
    sets = PAIRS + (("ivantsov", "den"),)
    print(f"fs = {FS:g} Hz；分子锁定 0/fc/slope，分母由两条约束锁定")
    print(f"{'fc':>8}" + "".join(f"{a + '+' + b:>22}" for a, b in sets))
    results = {}
    for fc in FC_LIST:
        row = []
        for pair in sets:
            match, err = design(pair, fc)
            results[(fc, pair)] = (match, err)
            row.append(f"{err:10.4f}")
        print(f"{fc:8.1f}" + "".join(f"{v:>22}" for v in row))
    print()

    for pair in sets:
        errs = [results[(fc, pair)][1] for fc in FC_LIST]
        print(f"{pair[0] + '+' + pair[1]:>22}: "
              f"max {np.max(errs):8.4f}  全部 " +
              " ".join(f"{e:.4f}" for e in errs))
    print()

    feasible_pairs = [pr for pr in PAIRS
                      if any(np.isfinite(results[(fc, pr)][1]) for fc in FC_LIST)]
    if not feasible_pairs:
        print("严格残差检查下，所有硬约束对在 5 个 fc 均无解。")
        print("结论：5 条硬约束不相容；改用插值视角（lowpass_1st_interp.py）。")
        return

    best_pair = min(feasible_pairs, key=lambda pr: _rank(pr, results))
    print(f"推荐约束对（按可解档数和有限误差）：{best_pair[0]} + {best_pair[1]}")

    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    _figure(results, out, best_pair)
    print(f"图：{out / 'lowpass_1st_constraints.png'}")


def _rank(pair, results):
    errs = [results[(fc, pair)][1] for fc in FC_LIST]
    finite = [e for e in errs if np.isfinite(e)]
    return (-len(finite), max(finite) if finite else np.inf)


def _figure(results, out, best_pair):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    fig, (ax_r, ax_e) = plt.subplots(
        2, 1, figsize=(10.5, 7.4), sharex=False,
        gridspec_kw={"height_ratios": (2.0, 1.2), "hspace": 0.28})
    f = np.geomspace(10.0, FS / 2.0, 4001)
    w = 2.0 * PI * f / FS
    skipped = []
    for j, fc in enumerate(FC_LIST):
        match, _ = results[(fc, best_pair)]
        if match is None:
            skipped.append(fc)
            continue
        color = f"C{j}"
        hd2 = match.power(w)
        ha2 = analog_P(w, 2.0 * PI * fc / FS)
        ax_r.semilogx(f, 10.0 * np.log10(np.maximum(ha2, 1e-300)),
                      color=color, lw=3.4, ls="--", alpha=0.30)
        ax_r.semilogx(f, 10.0 * np.log10(np.maximum(hd2, 1e-300)),
                      color=color, lw=1.3, label=f"{fc:g} Hz")
        err, _, ww = error_db(match, fc, FS, 200001)
        ax_e.plot(ww / PI * FS / 2.0, np.maximum(np.abs(err), 1e-12),
                  color=color, lw=1.2)

    ax_r.set_xscale("log")
    ax_r.set_xlim(10.0, FS / 2.0)
    ax_r.set_ylabel("幅度 / dB")
    ax_r.set_ylim(-90.0, 5.0)
    ax_r.grid(True, which="both", alpha=0.25)
    ax_r.legend(fontsize=8, title="fc", ncol=2)

    ax_e.set_ylabel("|ΔdB|")
    ax_e.set_yscale("log")
    ax_e.set_ylim(1e-6, 1e2)
    ax_e.set_xlim(0.0, FS / 2.0)
    ax_e.set_xlabel("频率 / Hz")
    ax_e.grid(True, which="both", alpha=0.25)

    note = "" if not skipped else "；无解：" + " / ".join(f"{x:g}" for x in skipped) + " Hz"
    ax_r.set_title(f"一阶低通，约束 {best_pair[0]} + {best_pair[1]}" + note)
    fig.tight_layout()
    fig.savefig(out / "lowpass_1st_constraints.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
