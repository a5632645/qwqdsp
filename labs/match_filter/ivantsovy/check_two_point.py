"""自检：一阶高通 / 低通 / shelf 的双点匹配（中点与 Nyquist）。

验证四件事：

1. 闭式解与直接解方程一致（用各自的自由参数做独立求根，二者互相印证）；
2. 匹配点在所有尺度上都精确到机器精度；
3. 稳定性、``sigma <= 2`` 的上界、``fc >= fs/2`` 的拒绝；
4. 各类型的结构性事实：高通的低频斜率恰好缩到 ``omega/u`` 倍，低通的
   DC 增益恰好为 1（因此低通实际精确匹配 DC、截止频率、Nyquist 三点），
   shelf 的 DC 增益自动精确且 boost 与 cut 严格互逆。

用法: python qwqdsp/labs/match_filter/ivantsovy/check_two_point.py
"""
from __future__ import annotations

import numpy as np
from scipy import optimize

import ivantsovy as I
import two_point as TP

TEST_CUTOFFS = tuple(float(v) for v in np.geomspace(20.0, 23900.0, 5))
"""测试截止频率：20 Hz … 23.9 kHz，对数 5 点（含端点），fs 统一 48 kHz。"""

CASES = tuple((fc, 48000.0) for fc in TEST_CUTOFFS)

POLES = np.array([-1.0])
KINDS = (
    ("高通", np.array([0.0]), TP.highpass),
    ("低通", np.array([]), TP.lowpass),
)

SHELF_KINDS = (
    ("高 shelf", "high", TP.high_shelf),
    ("低 shelf", "low", TP.low_shelf),
)
SHELF_CASES = CASES
SHELF_BOOSTS = (12.0, -12.0, 6.0, -6.0)

ZETAS = (0.5, 1.0 / np.sqrt(2.0), 0.9, 1.2)
ORDER2_CASES = CASES
ORDER2_KINDS = (("二阶高通", "hp", TP.highpass2), ("二阶低通", "lp", TP.lowpass2))
BAND_KINDS = (("二阶带通", "bp", TP.bandpass), ("二阶带阻", "bs", TP.bandstop))
SHELF2_QS = (np.sqrt(2.0) / 2.0, 6.0)
"""二阶 shelf 的两组阻尼：Q = 0.707（Butterworth）与 Q = 6（谐振）。"""

SHELF2_KINDS = (("二阶高 shelf", "high", TP.high_shelf2),
                ("二阶低 shelf", "low", TP.low_shelf2))


def require(cond: bool, message: str) -> None:
    if not cond:
        raise AssertionError(message)


def analog_mag(zeros, f, fc: float) -> np.ndarray:
    return I.analog_magnitude(zeros, POLES, f, fc)


def check_exact_two_points() -> None:
    print("双点匹配（截止频率处 1/sqrt(2)，Nyquist 处与模拟原型相等）")
    worst = 0.0
    for label, zeros, build in KINDS:
        print(f"  [{label}]")
        for fc, fs in CASES:
            design = build(fc, fs)
            omega, sigma = TP.omega_sigma(fc, fs)
            h_cut = float(I.digital_magnitude(design, np.array([fc]), fs)[0])
            h_nyq = float(I.digital_magnitude(design, np.array([fs / 2.0]), fs)[0])
            a_nyq = float(analog_mag(zeros, np.array([fs / 2.0]), fc)[0])
            rel_cut = abs(h_cut * np.sqrt(2.0) - 1.0)
            rel_nyq = abs(h_nyq / a_nyq - 1.0)
            worst = max(worst, rel_cut, rel_nyq)
            require(rel_cut < 1e-11 and rel_nyq < 1e-11,
                    f"{label} fc={fc}, fs={fs}: 匹配点误差 {rel_cut:.2e} / {rel_nyq:.2e}")
            print(f"    fc={fc:8.0f}  fs={fs:6.0f}  u={fs / fc:9.4f}  "
                  f"omega={omega:12.6f}  sigma={sigma:9.6f}  "
                  f"截止 {rel_cut:.1e}  Nyquist {rel_nyq:.1e}")
    print(f"PASS  两种类型的截止频率与 Nyquist 同时精确（最大相对误差 {worst:.1e}）\n")


def free_highpass(fc: float, fs: float) -> tuple[float, float]:
    """高通 ``g (1-z^-1)/(1+beta z^-1)`` 的独立求根，返回 ``(omega, sigma)``。"""
    u = fs / fc
    target = float(analog_mag(KINDS[0][1], np.array([fs / 2.0]), fc)[0])
    theta = 2.0 * np.pi / u

    def residual(beta: float) -> float:
        gain = target * (1.0 - beta) / 2.0
        return gain * 2.0 * np.sin(theta / 2.0) / abs(1.0 + beta * np.exp(-1j * theta)) \
            - 1.0 / np.sqrt(2.0)

    beta = optimize.brentq(residual, -1.0 + 1e-12, 1.0 - 1e-12, xtol=1e-15, rtol=1e-15)
    omega = 2.0 * np.pi * target * (1.0 - beta) / 2.0 / (1.0 + beta)
    sigma = np.sqrt(max(np.pi ** 2 * ((1.0 - beta) / (1.0 + beta)) ** 2 - omega ** 2, 0.0))
    return float(omega), float(sigma)


def lowpass_root_quality(fc: float, fs: float) -> tuple[float, float]:
    """低通两方程的根质量：闭式解处的残差与雅可比行列式（非零即局部唯一）。

    贴近 Nyquist 时（``u -> 2``）方程极平，多初值的迭代可能停在伪极小，
    因此这里直接检查根条件，而不是从远处迭代。
    """
    u = fs / fc
    target = float(analog_mag(KINDS[1][1], np.array([fs / 2.0]), fc)[0])
    theta = 2.0 * np.pi / u

    def eqs(x):
        alpha, beta = x
        scale = (1.0 + beta) / (1.0 + alpha)
        cut = scale * abs(1.0 + alpha * np.exp(-1j * theta)) \
            / abs(1.0 + beta * np.exp(-1j * theta)) - 1.0 / np.sqrt(2.0)
        nyq = scale * (1.0 - alpha) / (1.0 - beta) - target
        return np.array([cut, nyq])

    x0 = np.array([I.phi(0.0, TP.omega_sigma(fc, fs)[1]),
                   I.phi(TP.omega_sigma(fc, fs)[0], TP.omega_sigma(fc, fs)[1])])
    step = 1e-7
    jac = np.zeros((2, 2))
    for j in range(2):
        delta = np.zeros(2)
        delta[j] = step
        jac[:, j] = (eqs(x0 + delta) - eqs(x0 - delta)) / (2.0 * step)
    return float(np.max(np.abs(eqs(x0)))), float(abs(np.linalg.det(jac)))


def check_uniqueness() -> None:
    worst = 0.0
    for fc, fs in CASES:
        omega, sigma = TP.omega_sigma(fc, fs)
        om, sg = free_highpass(fc, fs)
        err = max(abs(om / omega - 1.0), abs(sg / sigma - 1.0))
        worst = max(worst, err)
        require(err < 1e-8,
                f"高通 fc={fc}, fs={fs}: 独立求根与闭式解不一致 {err:.2e}")

        res, det = lowpass_root_quality(fc, fs)
        require(res < 1e-12 and det > 1e-6,
                f"低通 fc={fc}: 残差 {res:.1e} / det {det:.1e}")
    print(f"PASS  闭式解唯一：高通用 brentq 独立求根（偏差 {worst:.1e}），"
          f"低通为两方程的孤立根（残差 < 1e-12 且雅可比非奇异）\n")


def check_parameters() -> None:
    for fc, fs in CASES:
        omega, sigma = TP.omega_sigma(fc, fs)
        require(sigma <= 2.0 + 1e-12, f"fc={fc}: sigma={sigma} 超过 2")
        require(omega < fs / fc, f"fc={fc}: omega 应小于 fs/fc")
    _, sigma = TP.omega_sigma(1e-9, 48000.0)
    require(abs(sigma - 2.0) < 1e-9, f"sigma 的低频极限应为 2，得到 {sigma}")
    omega, _ = TP.omega_sigma(1e-9, 48000.0)
    require(abs(omega / (48000.0 / 1e-9) - 1.0) < 1e-9, "omega 的低频极限应为 fs/fc")
    print("PASS  sigma <= 2 恒成立；u -> infinity 时 (omega, sigma) -> (fs/fc, 2)\n")


def expand(design: I.Design) -> tuple[np.ndarray, np.ndarray]:
    """把因子串联展开成 ``(b, a)``，并归一化到 ``a[0] = 1``。"""
    b = np.array([design.gain])
    for fac in design.num:
        b = np.convolve(b, fac.poly * (fac.gain / fac.norm))
    a = np.array([1.0])
    for fac in design.den:
        a = np.convolve(a, fac.poly * (fac.gain / fac.norm))
    return b / a[0], a / a[0]


def check_paper_coefficients() -> None:
    """把双点解代入 (2.8) / Table 2.1 的显式系数，逐项核对。"""
    two_pi = 2.0 * np.pi

    for fc, fs in CASES:
        omega, sigma = TP.omega_sigma(fc, fs)
        alpha, beta = I.phi(0.0, sigma), I.phi(omega, sigma)

        b, a = expand(TP.lowpass(fc, fs))
        g = 1.0 / (1.0 + alpha)
        b_ref = np.array([g * (1.0 + beta), g * (1.0 + beta) * alpha])
        require(np.max(np.abs(a - [1.0, beta])) < 1e-12
                and np.max(np.abs(b - b_ref)) < 1e-12,
                f"低通 fc={fc}: 系数与 Table 2.1 不一致 b={b}, a={a}")

        b, a = expand(TP.highpass(fc, fs))
        gain = omega / two_pi * (1.0 + beta)
        require(np.max(np.abs(a - [1.0, beta])) < 1e-12
                and np.max(np.abs(b - [gain, -gain])) < 1e-12,
                f"高通 fc={fc}: 系数与 Table 2.1 不一致 b={b}, a={a}")
    print("PASS  双点解逐系数等于论文 (2.8) / Table 2.1（alpha、beta、G 用双点 omega/sigma）\n")


def check_structural_points() -> None:
    """高通：低频斜率缩到 omega/u；低通：DC 增益恒为 1。"""
    for fc, fs in CASES:
        omega, _ = TP.omega_sigma(fc, fs)
        f = np.array([fc * 1e-4])
        hd = float(I.digital_magnitude(TP.highpass(fc, fs), f, fs)[0])
        ha = float(analog_mag(KINDS[0][1], f, fc)[0])
        err = abs(hd / ha / (omega / (fs / fc)) - 1.0)
        require(err < 1e-6, f"高通 fc={fc}: 低频比值偏离 omega/u（{err:.2e}）")

        dc = float(I.digital_magnitude(TP.lowpass(fc, fs), np.array([0.0]), fs)[0])
        require(abs(dc - 1.0) < 1e-15, f"低通 fc={fc}: DC 增益 {dc} 应为 1")
    print("PASS  高通低频增益比 = omega/u；低通 DC 增益 = 1（机器精度）\n")


def check_stability() -> None:
    worst = 0.0
    count = 0
    for fc in np.geomspace(1.0, 23990.0, 400):
        built = [build(float(fc), 48000.0) for _, _, build in KINDS]
        for _, _, build in SHELF_KINDS:
            for boost in (12.0, -12.0):
                built.append(build(float(fc), 48000.0, boost))
        for _, _, build in ORDER2_KINDS + BAND_KINDS:
            for zeta in (0.5, 1.0 / np.sqrt(2.0), 1.2):
                built.append(build(float(fc), 48000.0, zeta))
        for _, _, build in SHELF2_KINDS:
            for Q in SHELF2_QS:
                built.append(build(float(fc), 48000.0, 12.0, 1.0 / (2.0 * Q)))
        for design in built:
            _, poles = I.digital_zeros_poles(design)
            radius = float(np.max(np.abs(poles)))
            worst = max(worst, radius)
            require(radius < 1.0, f"fc={fc}: 极点在单位圆外（{radius}）")
            count += 1
    print(f"PASS  全部 {count} 组稳定（最大极点模 {worst:.9f}）\n")


def check_infeasible() -> None:
    for fc in (24000.0, 30000.0):
        try:
            TP.omega_sigma(fc, 48000.0)
        except ValueError:
            continue
        raise AssertionError(f"fc={fc} 在 fs=48 kHz 下本应无解")
    try:
        I.design(np.array([0.0]), POLES, 48000.0, 20000.0, 1.2, prewarp=True, omega=2.0)
    except ValueError:
        pass
    else:
        raise AssertionError("omega 与 prewarp 同时给出本应报错")
    print("PASS  fc >= fs/2 无解；omega 与 prewarp 互斥\n")


def shelf_analog(kind: str, f, fc: float, boost_db: float) -> np.ndarray:
    """模拟 shelf 幅度（中点增益 ``sqrt(g)`` 在 ``fc``），按 Table 2.1 的定义式。"""
    g = 10.0 ** (boost_db / 20.0)
    root = np.sqrt(g)
    x = np.atleast_1d(np.asarray(f, dtype=float)) / fc
    if kind == "high":
        return np.abs(1.0 + 1j * root * x) / np.abs(1.0 + 1j * x / root)
    return g * np.abs(1.0 + 1j * x / root) / np.abs(1.0 + 1j * root * x)


def free_shelf(kind: str, fc: float, fs: float, boost_db: float,
               seed: tuple[float, float]) -> tuple[float, float]:
    """shelf 的独立求根（自由 ``alpha, beta``，DC 归一化固定），返回 ``(omega, sigma)``。"""
    g = 10.0 ** (boost_db / 20.0)
    mid = np.sqrt(g)
    target = float(shelf_analog(kind, fs / 2.0, fc, boost_db)[0])
    theta = 2.0 * np.pi * (fc / fs)
    dc = g if kind == "low" else 1.0

    def eqs(x):
        alpha, beta = x
        if not (-1.0 < alpha < 1.0 and -1.0 < beta < 1.0):
            return [1e3, 1e3]
        scale = dc * (1.0 + beta) / (1.0 + alpha)
        cut = scale * abs(1.0 + alpha * np.exp(-1j * theta)) \
            / abs(1.0 + beta * np.exp(-1j * theta)) - mid
        nyq = scale * (1.0 - alpha) / (1.0 - beta) - target
        return [cut, nyq]

    out = optimize.least_squares(
        eqs, seed,
        bounds=([-0.999999999, -0.999999999], [0.999999999, 0.999999999]),
        xtol=1e-15, ftol=1e-15, gtol=1e-15)
    alpha, beta = out.x
    require(float(np.max(np.abs(eqs(out.x)))) < 1e-12,
            f"{kind} fc={fc}: shelf 独立求根未收敛")
    p2 = (np.pi * (1.0 - alpha) / (1.0 + alpha)) ** 2
    q2 = (np.pi * (1.0 - beta) / (1.0 + beta)) ** 2
    a2, b2 = (g, 1.0 / g) if kind == "high" else (1.0 / g, g)
    omega2 = (p2 - q2) / (a2 - b2)
    sigma2 = p2 - a2 * omega2
    require(omega2 > 0.0 and sigma2 > 0.0, f"{kind} fc={fc}: 反推 (omega, sigma) 非正")
    return float(np.sqrt(omega2)), float(np.sqrt(sigma2))


def check_shelves() -> None:
    print("shelf 双点匹配（中点 sqrt(g) 与 Nyquist 精确；DC 自动精确）")
    worst = 0.0
    for label, kind, build in SHELF_KINDS:
        for boost in SHELF_BOOSTS:
            g = 10.0 ** (boost / 20.0)
            root = np.sqrt(g)
            errs = [0.0, 0.0, 0.0, 0.0]
            for fc, fs in SHELF_CASES:
                design = build(fc, fs, boost)
                omega, sigma = TP.omega_sigma(fc, fs)
                h_mid = float(I.digital_magnitude(design, np.array([fc]), fs)[0])
                h_nyq = float(I.digital_magnitude(design, np.array([fs / 2.0]), fs)[0])
                h_dc = float(I.digital_magnitude(design, np.array([0.0]), fs)[0])
                e_mid = abs(h_mid / np.sqrt(g) - 1.0)
                e_nyq = abs(h_nyq / float(shelf_analog(kind, fs / 2.0, fc, boost)[0]) - 1.0)
                e_dc = abs(h_dc / (g if kind == "low" else 1.0) - 1.0)

                alpha = I.phi(omega * (root if kind == "high" else 1.0 / root), sigma)
                beta = I.phi(omega * (1.0 / root if kind == "high" else root), sigma)
                gain = 1.0 / (1.0 + alpha) if kind == "high" else g / (1.0 + alpha)
                b_ref = np.array([gain * (1.0 + beta), gain * (1.0 + beta) * alpha])
                b, a = expand(design)
                e_coef = max(float(np.max(np.abs(b - b_ref))),
                             float(np.max(np.abs(a - [1.0, beta]))))

                require(e_mid < 1e-11 and e_nyq < 1e-11 and e_dc < 1e-13
                        and e_coef < 1e-12,
                        f"{label} boost={boost}, fc={fc}: "
                        f"中点 {e_mid:.1e} / Nyquist {e_nyq:.1e} / "
                        f"DC {e_dc:.1e} / 系数 {e_coef:.1e}")
                errs = [max(errs[0], e_mid), max(errs[1], e_nyq),
                        max(errs[2], e_dc), max(errs[3], e_coef)]
            worst = max(worst, errs[0], errs[1])
            print(f"  [{label}] boost={boost:+6.1f} dB: 中点 {errs[0]:.1e}  "
                  f"Nyquist {errs[1]:.1e}  DC {errs[2]:.1e}  系数 {errs[3]:.1e}")
    print(f"PASS  shelf 中点与 Nyquist 精确（最大相对误差 {worst:.1e}）\n")


def check_shelf_inversion() -> None:
    worst = 0.0
    for _, kind, build in SHELF_KINDS:
        for boost in (12.0, 6.0):
            for fc, fs in SHELF_CASES:
                f = np.geomspace(1.0, fs / 2.0, 500)
                up = I.digital_magnitude(build(fc, fs, boost), f, fs)
                down = I.digital_magnitude(build(fc, fs, -boost), f, fs)
                err = float(np.max(np.abs(up * down - 1.0)))
                worst = max(worst, err)
                require(err < 1e-12,
                        f"{kind} fc={fc}: boost 与 cut 不互逆（{err:.1e}）")
    print(f"PASS  shelf 的 boost 与 cut 严格互逆（|H+b|·|H-b| = 1，偏差 {worst:.1e}）\n")


def check_shelf_uniqueness() -> None:
    worst = 0.0
    for _, kind, _ in SHELF_KINDS:
        for boost in (12.0, -6.0):
            for fc, fs in ((2000.0, 48000.0), (20000.0, 48000.0)):
                omega, sigma = TP.omega_sigma(fc, fs)
                for seed in ((0.2, -0.8), (0.8, 0.5), (0.5, -0.2)):
                    om, sg = free_shelf(kind, fc, fs, boost, seed)
                    err = max(abs(om / omega - 1.0), abs(sg / sigma - 1.0))
                    worst = max(worst, err)
                    require(err < 1e-8,
                            f"{kind} boost={boost}, fc={fc}: 独立求根与闭式解不一致 "
                            f"{err:.2e}")
    print(f"PASS  shelf 闭式解唯一：自由 (alpha, beta) 求根一致（最大相对偏差 "
          f"{worst:.1e}）\n")


def report_shelves() -> None:
    print("shelf 幅度对照（+12 dB；−12 dB 由 boost/cut 互逆性给出相同数值）")
    print(f"{'kind':>7}{'fc':>8}{'u':>10}{'max|ΔdB|':>10}{'最差处':>11}{'模拟电平':>10}")
    for label, kind, build in SHELF_KINDS:
        for fc, fs in SHELF_CASES:
            design = build(fc, fs, 12.0)
            f = np.linspace(0.0, fs / 2.0, 60001)
            ha = shelf_analog(kind, f, fc, 12.0)
            hd = I.digital_magnitude(design, f, fs)
            dev = 20.0 * np.log10(np.maximum(hd, 1e-300) / ha)
            j = int(np.argmax(np.abs(dev)))
            print(f"{label:>7}{fc:8.0f}{fs / fc:10.3f}{np.max(np.abs(dev)):10.4f}"
                  f"{f[j]:10.0f}Hz{20.0 * np.log10(ha[j]):10.2f}")
    print()


def analog2(kind: str, f, fc: float, zeta: float) -> np.ndarray:
    """二阶原型幅度（(3.1) 的等价形式）。"""
    x = np.atleast_1d(np.asarray(f, dtype=float)) / fc
    b = 2.0 * zeta ** 2 - 1.0
    den = np.sqrt(x ** 4 + 2.0 * x ** 2 * b + 1.0)
    return x ** 2 / den if kind == "hp" else 1.0 / den


def design2(kind: str, omega: float, sigma: float, fc: float, fs: float,
            zeta: float) -> I.Design:
    zeros = np.array([0.0, 0.0]) if kind == "hp" else np.array([])
    return I.design(zeros, TP.second_order_poles(zeta), fs, fc, sigma, omega=omega)


def shelf2_analog(kind: str, f, fc: float, boost_db: float, zeta: float) -> np.ndarray:
    zeros, poles, gain = TP.shelf_prototype2(kind, boost_db, zeta)
    return I.analog_magnitude(zeros, poles, np.atleast_1d(np.asarray(f, dtype=float)),
                              fc, gain)


def check_shelves2() -> None:
    print("二阶 shelf（Q = 0.707 / 6）：中点 sqrt(g)、Nyquist、DC 精确")
    for Q in SHELF2_QS:
        zeta = 1.0 / (2.0 * Q)
        errs = [0.0, 0.0, 0.0, 0.0]
        for label, kind, build in SHELF2_KINDS:
            for boost in (12.0, -12.0, 6.0, -6.0):
                g = 10.0 ** (boost / 20.0)
                root = g ** 0.25
                for fc, fs in CASES:
                    design = build(fc, fs, boost, zeta)
                    omega, sigma = TP.omega_sigma(fc, fs)
                    h_mid = float(I.digital_magnitude(design, np.array([fc]), fs)[0])
                    h_nyq = float(I.digital_magnitude(design, np.array([fs / 2.0]), fs)[0])
                    h_dc = float(I.digital_magnitude(design, np.array([0.0]), fs)[0])
                    e_mid = abs(h_mid / np.sqrt(g) - 1.0)
                    e_nyq = abs(h_nyq / float(shelf2_analog(kind, fs / 2.0, fc, boost, zeta)[0])
                                - 1.0)
                    e_dc = abs(h_dc / (g if kind == "low" else 1.0) - 1.0)

                    alpha1, alpha2 = I.phi12(omega * (root if kind == "high" else 1.0 / root),
                                             zeta, sigma)
                    beta1, beta2 = I.phi12(omega * (1.0 / root if kind == "high" else root),
                                           zeta, sigma)
                    gain = (1.0 / (1.0 + alpha1 + alpha2) if kind == "high"
                            else g / (1.0 + alpha1 + alpha2))
                    b_ref = gain * (1.0 + beta1 + beta2) * np.array([1.0, alpha1, alpha2])
                    b, a = expand(design)
                    e_coef = max(float(np.max(np.abs(b - b_ref))),
                                 float(np.max(np.abs(a - [1.0, beta1, beta2]))))

                    require(e_mid < 1e-9 and e_nyq < 1e-9 and e_dc < 1e-12
                            and e_coef < 1e-10,
                            f"{label} Q={Q}, boost={boost}, fc={fc}: "
                            f"{e_mid:.1e}/{e_nyq:.1e}/{e_dc:.1e}/{e_coef:.1e}")
                    errs = [max(errs[0], e_mid), max(errs[1], e_nyq),
                            max(errs[2], e_dc), max(errs[3], e_coef)]
        print(f"  Q={Q:.4f} (zeta={zeta:.4f}): 中点 {errs[0]:.1e}  Nyquist {errs[1]:.1e}  "
              f"DC {errs[2]:.1e}  系数 {errs[3]:.1e}")
    print("PASS  二阶 shelf 与一阶共用同一组 (omega, sigma)\n")


def check_shelves2_reciprocity() -> None:
    worst = 0.0
    for Q in SHELF2_QS:
        zeta = 1.0 / (2.0 * Q)
        for _, kind, build in SHELF2_KINDS:
            for fc, fs in ((4065.0, 48000.0), (23900.0, 48000.0)):
                f = np.geomspace(20.0, fs / 2.0, 500)
                up = I.digital_magnitude(build(fc, fs, 12.0, zeta), f, fs)
                down = I.digital_magnitude(build(fc, fs, -12.0, zeta), f, fs)
                err = float(np.max(np.abs(up * down - 1.0)))
                worst = max(worst, err)
                require(err < 1e-12, f"{kind} Q={Q} fc={fc}: 互逆性 {err:.1e}")
    print(f"PASS  二阶 shelf 的 boost / cut 严格互逆（偏差 {worst:.1e}）\n")


def report_shelves2() -> None:
    print("二阶 shelf 幅度对照（+12 dB；|Ha| >= 0.1 内 max|ΔdB|）")
    print(f"{'Q':>7}{'kind':>13}{'fc':>8}{'u':>10}{'max|ΔdB|':>10}{'最差处':>11}")
    for Q in SHELF2_QS:
        zeta = 1.0 / (2.0 * Q)
        for label, kind, build in SHELF2_KINDS:
            for fc, fs in CASES:
                design = build(fc, fs, 12.0, zeta)
                f = np.linspace(0.0, fs / 2.0, 120001)
                ha = shelf2_analog(kind, f, fc, 12.0, zeta)
                hd = I.digital_magnitude(design, f, fs)
                dev = 20.0 * np.log10(np.maximum(hd, 1e-300) / np.maximum(ha, 1e-300))
                mask = ha >= 0.1
                j = int(np.argmax(np.abs(dev[mask])))
                print(f"{Q:7.3f}{label:>13}{fc:8.0f}{fs / fc:10.3f}"
                      f"{np.max(np.abs(dev[mask])):10.4f}{f[mask][j]:10.0f}Hz")
    print()


def check_second_order() -> None:
    print("二阶：|H(fc)| = |Ha(fc)|、|H(fs/2)| = |Ha(fs/2)|（与 zeta、u 无关）")
    for zeta in ZETAS:
        errs = [0.0, 0.0, 0.0]
        for fc, fs in ORDER2_CASES:
            omega, sigma = TP.omega_sigma(fc, fs)
            for label, kind, build in ORDER2_KINDS:
                design = build(fc, fs, zeta)
                h_cut = float(I.digital_magnitude(design, np.array([fc]), fs)[0])
                h_nyq = float(I.digital_magnitude(design, np.array([fs / 2.0]), fs)[0])
                a_cut = float(analog2(kind, fc, fc, zeta)[0])
                a_nyq = float(analog2(kind, fs / 2.0, fc, zeta)[0])
                e_cut = abs(h_cut / a_cut - 1.0)
                e_nyq = abs(h_nyq / a_nyq - 1.0)
                # u = 2400（fc = 20 Hz）时 biquad 系数相差 3 个数量级，
                # 双精度求值本身只剩 ~1e-11 的相对精度，故阈值取 1e-9。
                require(e_cut < 1e-9 and e_nyq < 1e-9,
                        f"{label} zeta={zeta}, fc={fc}: 误差 {e_cut:.1e} / {e_nyq:.1e}")
                errs[0] = max(errs[0], e_cut)
                errs[1] = max(errs[1], e_nyq)

                if zeta < 1.0:
                    phi0 = I.phi(0.0, sigma)
                    beta1, beta2 = I.phi12(omega, zeta, sigma)
                    if kind == "hp":
                        alpha1, alpha2 = -2.0, 1.0
                        gain = omega ** 2 / (4.0 * np.pi ** 2)
                    else:
                        alpha1, alpha2 = 2.0 * phi0, phi0 ** 2
                        gain = 1.0 / (1.0 + alpha1 + alpha2)
                    b_ref = gain * (1.0 + beta1 + beta2) * np.array([1.0, alpha1, alpha2])
                    b, a = expand(design)
                    e_coef = max(float(np.max(np.abs(b - b_ref))),
                                 float(np.max(np.abs(a - [1.0, beta1, beta2]))))
                    require(e_coef < 1e-12,
                            f"{label} zeta={zeta}, fc={fc}: 系数与 Table 3.1 不一致 "
                            f"{e_coef:.1e}")
                    errs[2] = max(errs[2], e_coef)

                if zeta == 1.0 / np.sqrt(2.0):
                    require(abs(h_cut * np.sqrt(2.0) - 1.0) < 1e-9,
                            f"{label} fc={fc}: zeta=1/sqrt(2) 时 |H(fc)| 应为 1/sqrt(2)")
        print(f"  zeta={zeta:.3f}: |H(fc)| {errs[0]:.1e}  Nyquist {errs[1]:.1e}  "
              f"Table 3.1 系数 {errs[2]:.1e}")
    print("PASS  二阶梯形的 fc 与 Nyquist 同时精确（含实极点 zeta > 1）\n")


def check_second_order_uniqueness() -> None:
    worst = 0.0
    for zeta in (0.5, 1.0 / np.sqrt(2.0)):
        for fc, fs in ((2000.0, 48000.0), (20000.0, 48000.0)):
            omega, sigma = TP.omega_sigma(fc, fs)
            for label, kind, _ in ORDER2_KINDS:
                target_cut = float(analog2(kind, fc, fc, zeta)[0])
                target_nyq = float(analog2(kind, fs / 2.0, fc, zeta)[0])

                def eqs(x):
                    w, s = x
                    if w <= 0.0 or s <= 0.0:
                        return [1e3, 1e3]
                    design = design2(kind, w, s, fc, fs, zeta)
                    return [
                        float(I.digital_magnitude(design, np.array([fc]), fs)[0])
                        - target_cut,
                        float(I.digital_magnitude(design, np.array([fs / 2.0]), fs)[0])
                        - target_nyq,
                    ]

                sol, info, _, _ = optimize.fsolve(
                    eqs, [omega * 1.05, sigma * 0.95], full_output=True, xtol=1e-13)
                res = float(np.max(np.abs(info["fvec"])))
                require(res < 1e-10,
                        f"{label} zeta={zeta}, fc={fc}: 二阶独立求根残差 {res:.1e}")
                err = max(abs(sol[0] / omega - 1.0), abs(sol[1] / sigma - 1.0))
                worst = max(worst, err)
                require(err < 1e-4,
                        f"{label} zeta={zeta}, fc={fc}: 独立求根与闭式解不一致 {err:.1e}")
    print(f"PASS  二阶闭式解唯一：直接解 (omega, sigma) 两方程一致（最大相对偏差 "
          f"{worst:.1e}）\n")


def analog_band(kind: str, f, fc: float, zeta: float) -> np.ndarray:
    """二阶带通 / 带阻的原型幅度（峰值或陷波在 ``fc``）。"""
    s = 1j * np.atleast_1d(np.asarray(f, dtype=float)) / fc
    den = np.abs(s ** 2 + 2.0 * zeta * s + 1.0)
    return np.abs(2.0 * zeta * s) / den if kind == "bp" else np.abs(s ** 2 + 1.0) / den


def check_bands() -> None:
    print("二阶带通 / 带阻：中心点与 Nyquist 精确")
    for zeta in ZETAS:
        errs = [0.0, 0.0, 0.0]
        for fc, fs in ORDER2_CASES:
            omega, sigma = TP.omega_sigma(fc, fs)
            theta = 2.0 * np.pi * fc / fs
            for label, kind, build in BAND_KINDS:
                design = build(fc, fs, zeta)
                h_cut = float(I.digital_magnitude(design, np.array([fc]), fs)[0])
                h_nyq = float(I.digital_magnitude(design, np.array([fs / 2.0]), fs)[0])
                a_nyq = float(analog_band(kind, fs / 2.0, fc, zeta)[0])
                e_nyq = abs(h_nyq / a_nyq - 1.0)

                if kind == "bp":
                    e_center = abs(h_cut - 1.0)
                else:
                    e_center = h_cut
                    zeros, _ = I.digital_zeros_poles(design)
                    on = zeros[np.abs(zeros.imag) > 1e-12]
                    require(len(on) == 2, f"{label} fc={fc}: 应有共轭零点对")
                    require(float(np.max(np.abs(np.abs(on) - 1.0))) < 1e-12,
                            f"{label} fc={fc}: 零点不在单位圆上")
                    require(float(np.max(np.abs(np.abs(np.angle(on)) - theta))) < 1e-12,
                            f"{label} fc={fc}: 陷波频率偏移")
                require(e_center < 1e-9 and e_nyq < 1e-9,
                        f"{label} zeta={zeta}, fc={fc}: {e_center:.1e} / {e_nyq:.1e}")
                errs[0] = max(errs[0], e_center)
                errs[1] = max(errs[1], e_nyq)

                if zeta < 1.0:
                    phi0 = I.phi(0.0, sigma)
                    beta1, beta2 = I.phi12(omega, zeta, sigma)
                    if kind == "bp":
                        alpha1, alpha2 = phi0 - 1.0, -phi0
                        gain = (omega / (2.0 * np.pi)) * 2.0 * zeta / (2.0 + alpha1)
                    else:
                        alpha1, alpha2 = I.phi12(omega, 0.0, sigma)
                        gain = 1.0 / (1.0 + alpha1 + alpha2)
                    b_ref = gain * (1.0 + beta1 + beta2) * np.array([1.0, alpha1, alpha2])
                    b, a = expand(design)
                    e_coef = max(float(np.max(np.abs(b - b_ref))),
                                 float(np.max(np.abs(a - [1.0, beta1, beta2]))))
                    require(e_coef < 1e-12,
                            f"{label} zeta={zeta}, fc={fc}: 系数与 Table 3.1 不一致 "
                            f"{e_coef:.1e}")
                    errs[2] = max(errs[2], e_coef)
        print(f"  zeta={zeta:.3f}: 中心点 {errs[0]:.1e}  Nyquist {errs[1]:.1e}  "
              f"Table 3.1 系数 {errs[2]:.1e}")
    print("PASS  带通峰值 1 与带阻零点 0 都精确落在 fc，Nyquist 同时精确\n")


def report_bands() -> None:
    zeta = 1.0 / np.sqrt(2.0)
    print("二阶带通 / 带阻幅度对照（zeta = 1/sqrt(2)；|Ha| >= 0.5 内取最大偏差）")
    print(f"{'类型':>9}{'fc':>8}{'u':>10}{'|H(fc)|':>11}{'max|ΔdB|':>10}"
          f"{'最差处':>11}{'模拟电平':>10}")
    for label, kind, build in BAND_KINDS:
        for fc, fs in ORDER2_CASES:
            design = build(fc, fs, zeta)
            f = np.linspace(0.0, fs / 2.0, 120001)
            ha = analog_band(kind, f, fc, zeta)
            hd = I.digital_magnitude(design, f, fs)
            dev = 20.0 * np.log10(np.maximum(hd, 1e-300) / np.maximum(ha, 1e-300))
            mask = ha >= 0.5
            j = int(np.argmax(np.abs(dev[mask])))
            h_cut = float(I.digital_magnitude(design, np.array([fc]), fs)[0])
            print(f"{label:>9}{fc:8.0f}{fs / fc:10.3f}{h_cut:11.2e}"
                  f"{np.max(np.abs(dev[mask])):10.3f}{f[mask][j]:10.0f}Hz"
                  f"{20.0 * np.log10(ha[mask][j]):10.2f}")
    print()


def report_second_order() -> None:
    zeta = 1.0 / np.sqrt(2.0)
    print(f"二阶幅度对照（zeta = 1/sqrt(2)；ΔdB = 20 log10(Hd/Ha)）")
    print(f"{'类型':>9}{'fc':>8}{'u':>10}{'max|ΔdB|':>10}{'>=-20dB':>9}"
          f"{'最差处':>11}{'模拟电平':>10}")
    for label, kind, build in ORDER2_KINDS:
        for fc, fs in ORDER2_CASES:
            design = build(fc, fs, zeta)
            f = np.linspace(0.0, fs / 2.0, 60001)
            ha = analog2(kind, f, fc, zeta)
            hd = I.digital_magnitude(design, f, fs)
            dev = 20.0 * np.log10(np.maximum(hd, 1e-300) / np.maximum(ha, 1e-300))
            j = int(np.argmax(np.abs(dev)))
            mask = ha >= 0.1
            print(f"{label:>9}{fc:8.0f}{fs / fc:10.3f}{np.max(np.abs(dev)):10.4f}"
                  f"{np.max(np.abs(dev[mask])):9.4f}{f[j]:10.0f}Hz"
                  f"{20.0 * np.log10(ha[j]):10.1f}")
    print()


def report_highpass() -> None:
    print("高通幅度对照（ΔdB = 20 log10(Hd/Ha)；低频偏移 = 20 log10(omega/u)）")
    print(f"{'fc':>8}{'fs':>8}{'u':>9}{'sigma':>9}{'max|ΔdB|':>10}"
          f"{'>=fc/10':>10}{'>=0.9fc':>10}{'低频偏移':>10}")
    for fc, fs in CASES:
        design = TP.highpass(fc, fs)
        omega, sigma = TP.omega_sigma(fc, fs)
        full = np.concatenate([np.geomspace(fc * 1e-3, fs / 2.0, 40000), [fs / 2.0]])
        bands = {"full": full, "db20": full[full >= fc / 10.0],
                 "pass": full[full >= 0.9 * fc]}
        stats = {}
        for key, band in bands.items():
            ha = analog_mag(KINDS[0][1], band, fc)
            hd = I.digital_magnitude(design, band, fs)
            stats[key] = float(np.max(np.abs(20.0 * np.log10(hd / ha))))
        offset = 20.0 * np.log10(omega / (fs / fc))
        print(f"{fc:8.0f}{fs:8.0f}{fs / fc:9.3f}{sigma:9.6f}{stats['full']:10.4f}"
              f"{stats['db20']:10.4f}{stats['pass']:10.4f}{offset:10.3f}")
    print()


def report_lowpass() -> None:
    print("低通幅度对照（DC、截止频率、Nyquist 三点精确；残差集中在高频端）")
    print(f"{'fc':>8}{'fs':>8}{'u':>9}{'sigma':>9}{'max|ΔdB|':>10}"
          f"{'<=10fc':>10}{'<=0.7fc':>10}{'最差处':>12}{'模拟电平':>10}")
    for fc, fs in CASES:
        design = TP.lowpass(fc, fs)
        omega, sigma = TP.omega_sigma(fc, fs)
        full = np.concatenate([np.linspace(0.0, fs / 2.0, 60001)])
        ha = analog_mag(KINDS[1][1], full, fc)
        hd = I.digital_magnitude(design, full, fs)
        dev = 20.0 * np.log10(np.maximum(hd, 1e-300) / np.maximum(ha, 1e-300))
        worst = int(np.argmax(np.abs(dev)))
        stats = {"full": float(np.max(np.abs(dev)))}
        for key, band in (("db20", full <= 10.0 * fc), ("pass", full <= 0.7 * fc)):
            stats[key] = float(np.max(np.abs(dev[band])))
        print(f"{fc:8.0f}{fs:8.0f}{fs / fc:9.3f}{sigma:9.6f}{stats['full']:10.4f}"
              f"{stats['db20']:10.4f}{stats['pass']:10.4f}"
              f"{full[worst]:11.0f}Hz{20.0 * np.log10(ha[worst]):10.1f}")
    print()


def report_summary() -> None:
    """统一口径：模拟原型 >= -20 dB 的频段内 max|ΔdB|，按 5 个测试截止频率排。"""
    zeta = 1.0 / np.sqrt(2.0)

    columns = (
        ("一阶HP", lambda fc, fs: TP.highpass(fc, fs),
         lambda f, fc: I.analog_magnitude(KINDS[0][1], POLES, f, fc)),
        ("一阶LP", lambda fc, fs: TP.lowpass(fc, fs),
         lambda f, fc: I.analog_magnitude(KINDS[1][1], POLES, f, fc)),
        ("高shelf", lambda fc, fs: TP.high_shelf(fc, fs, 12.0),
         lambda f, fc: shelf_analog("high", f, fc, 12.0)),
        ("低shelf", lambda fc, fs: TP.low_shelf(fc, fs, 12.0),
         lambda f, fc: shelf_analog("low", f, fc, 12.0)),
        ("二阶HP", lambda fc, fs: TP.highpass2(fc, fs, zeta),
         lambda f, fc: analog2("hp", f, fc, zeta)),
        ("二阶LP", lambda fc, fs: TP.lowpass2(fc, fs, zeta),
         lambda f, fc: analog2("lp", f, fc, zeta)),
        ("带通", lambda fc, fs: TP.bandpass(fc, fs, zeta),
         lambda f, fc: analog_band("bp", f, fc, zeta)),
        ("带阻", lambda fc, fs: TP.bandstop(fc, fs, zeta),
         lambda f, fc: analog_band("bs", f, fc, zeta)),
    )

    print("汇总（模拟原型 >= -20 dB 的频段内 max|ΔdB|；shelf 为 +12 dB，二阶 ζ=1/√2）")
    print(f"{'fc':>8}{'u':>10}" + "".join(f"{name:>9}" for name, _, _ in columns))
    for fc, fs in CASES:
        f = np.linspace(0.0, fs / 2.0, 120001)
        cells = []
        for _, build, analog in columns:
            ha = analog(f, fc)
            hd = I.digital_magnitude(build(fc, fs), f, fs)
            dev = 20.0 * np.log10(np.maximum(hd, 1e-300) / np.maximum(ha, 1e-300))
            mask = ha >= 0.1
            cells.append(float(np.max(np.abs(dev[mask]))))
        print(f"{fc:8.0f}{fs / fc:10.3f}" + "".join(f"{v:9.4f}" for v in cells))
    print()


def main() -> int:
    check_exact_two_points()
    check_uniqueness()
    check_parameters()
    check_paper_coefficients()
    check_structural_points()
    check_shelves()
    check_shelf_inversion()
    check_shelf_uniqueness()
    check_second_order()
    check_second_order_uniqueness()
    check_bands()
    check_shelves2()
    check_shelves2_reciprocity()
    check_stability()
    check_infeasible()
    report_highpass()
    report_lowpass()
    report_shelves()
    report_second_order()
    report_bands()
    report_shelves2()
    report_summary()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
