"""NonIntegerSRC (jatinchowdhury18) 的逐行 Python 复刻 + 行为测量。

来源:
    https://github.com/jatinchowdhury18/NonIntegerSRC
    src/HPResampler.h            —— 活动分支 (`#else`) 与禁用分支 (`#if 0`)
    src/src_utils/HPFilters.h    —— FilterSpec（论文 Table 1 系数）+ 两个 filter bank
    src/src_utils/FastMath.h     —— fast sin/cos、vSum
论文: Holters & Parker, DAFx-18, "A Combined Model for a Bucket Brigade Device and its Input and
      Output Filters" — https://www.dafx.de/paper-archive/2018/papers/DAFx2018_paper_12.pdf

这里**逐行照抄** C++ 的数值路径（SSE 被展平成 4 个复标量），不做任何"修正"，
以便判定该实现的行为。用法:

    python noninteger_src_repro.py
"""

import numpy as np

# ------------------------------------------------------------
# FilterSpec：论文 Table 1 的系数（2 个共轭极点对；实极点 p1 / 实残差 r1 被丢弃）
# 输入滤波器的残差被整体缩放 1/12.6271（极点未缩放）——见 notes 的核对表
# ------------------------------------------------------------
iRoot = np.array([-10329.2715 - 329.848j, -10329.2715 + 329.848j,
                  366.990557 - 1811.4318j, 366.990557 + 1811.4318j])
iPole = np.array([-55482.0 - 25082.0j, -55482.0 + 25082.0j,
                  -26292.0 - 59437.0j, -26292.0 + 59437.0j])
oRoot = np.array([-11256.0 - 99566.0j, -11256.0 + 99566.0j,
                  -13802.0 - 24606.0j, -13802.0 + 24606.0j])
oPole = np.array([-51468.0 - 21437.0j, -51468.0 + 21437.0j,
                  -26276.0 - 59699.0j, -26276.0 + 59699.0j])


class OutputFilterBank:
    """`OutputFilterBank`（HPFilters.h）的照抄。Ts 为 filter 内部采样周期。"""

    def __init__(self, sample_time):
        self.Ts = sample_time
        self.g_coef = oRoot / oPole
        self.x = np.zeros(4, complex)

    def calc_h0(self):
        """源码里定义了但从未被调用。"""
        return -self.g_coef.real.sum()

    def set_freq(self, freq):
        freq_factor = freq / 9500.0                       # originalCutoff
        self.pole_corr = np.exp(oPole * freq_factor * self.Ts)
        self.pole_corr_angle = np.angle(self.pole_corr)
        self.amult = self.g_coef * self.pole_corr

    def set_time(self, tn):
        self.Gcalc = self.amult * self.pole_corr ** (1.0 - tn)

    def set_delta(self, delta):
        self.Aplus = np.exp(1j * self.pole_corr_angle * delta)

    def calc_g(self):
        self.Gcalc = self.Gcalc * self.Aplus

    def process(self, u):
        self.x = self.pole_corr * self.x + u


class InputFilterBank:
    """`InputFilterBank`（HPFilters.h）的照抄，只被 `#if 0` 分支使用。"""

    def __init__(self, sample_time):
        self.Ts = sample_time
        self.x = np.zeros(4, complex)

    def set_freq(self, freq):
        freq_factor = freq / 9900.0                       # originalCutoff
        self.root_corr = iRoot * freq_factor
        self.pole_corr = np.exp(iPole * freq_factor * self.Ts)
        self.pole_corr_angle = np.angle(self.pole_corr)
        self.g_coef = self.root_corr * self.Ts

    def set_time(self, tn):
        self.Gcalc = self.g_coef * self.pole_corr ** tn

    def set_delta(self, delta):
        self.Aplus = np.exp(1j * self.pole_corr_angle * delta)

    def calc_g(self):
        self.Gcalc = self.Gcalc * self.Aplus

    def process(self, u):
        self.x = self.pole_corr * self.x + u


def hp_resampler(signal, sample_rate, ratio, variant="as-coded"):
    """`HPResampler::prepare` + 活动分支 `process` 的照抄。

    ratio = 输出采样率 / 输入采样率（src_test.cpp: 48k->96k 用 2.0，96k->48k 用 0.5）
    variant: "as-coded" 源码原样 | "state" 输出读 filter 状态 | "h0" 同时用 calcH0()
    """
    Ts = 1.0 / sample_rate
    Ts_in = 1.0 / (sample_rate * ratio)      # == 输出采样周期
    Ts_out = 1.0 / (sample_rate / ratio)     # == ratio^2 * 输出采样周期

    flt = OutputFilterBank(Ts_in)
    flt.set_freq(sample_rate * 0.5)
    flt.set_delta(Ts_out)

    tn, y_old, i = 0.0, 0.0, 0
    flt.set_time(tn)
    out = []
    while i < len(signal):
        acc = np.zeros(4, complex)
        while tn < Ts and i < len(signal):
            y = signal[i]; i += 1
            delta = y - y_old
            y_old = y
            flt.calc_g()
            acc += flt.Gcalc * delta
            tn += Ts_out
        tn -= Ts

        flt.process(acc)                                  # 状态 x 被更新……
        if variant == "as-coded":
            out.append(y_old + acc.real.sum())            # ……但这里用的是驱动项累加器
        elif variant == "state":
            out.append(y_old + flt.x.real.sum())
        elif variant == "h0":
            out.append(flt.calc_h0() * y_old + flt.x.real.sum())
        else:
            raise ValueError(variant)
    return np.array(out)


def hp_resampler_input_filter(signal, sample_rate, ratio):
    """`#if 0 // use input filter (Dirac)` 分支的照抄（读状态，是"对"的那一半）。"""
    Ts = 1.0 / sample_rate
    Ts_in = 1.0 / (sample_rate * ratio)

    flt = InputFilterBank(Ts)
    flt.set_freq(sample_rate * 0.5)
    flt.set_delta(Ts_in)

    tn, i = 0.0, 0
    flt.set_time(tn)
    out = []
    while i < len(signal):
        while tn < Ts:
            flt.calc_g()
            out.append(-0.95097 * (flt.Gcalc * flt.x).real.sum())
            tn += Ts_in
        tn -= Ts
        flt.process(signal[i]); i += 1
    return np.array(out)


def zoh_naive(signal, ratio):
    n = int(len(signal) * ratio)
    return signal[np.minimum((np.arange(n) / ratio).astype(int), len(signal) - 1)]


def bin_mag(sig, fs, f):
    n = len(sig); w = np.hanning(n)
    sp = np.abs(np.fft.rfft(sig * w)) * 2 / w.sum()
    return sp[np.argmin(np.abs(np.fft.rfftfreq(n, 1 / fs) - f))]


def ls_amp_resid(sig, fo, f, skip=400):
    """最小二乘拟合同频正弦，返回 (幅度, 残差 rms)。"""
    k = np.arange(len(sig)) * 1.0
    A = np.stack([np.sin(2 * np.pi * f * k / fo)[skip:],
                  np.cos(2 * np.pi * f * k / fo)[skip:]], 1)
    coef, *_ = np.linalg.lstsq(A, sig[skip:], rcond=None)
    return np.hypot(*coef), np.sqrt(((sig[skip:] - A @ coef) ** 2).mean())


def main():
    N = 8192
    print("=" * 78)
    print("A) 降采样 96k->48k，输入 30 kHz（必须被抗混叠滤掉；不滤会折到 18 kHz）")
    print("=" * 78)
    t = np.arange(N)
    x = np.sin(2 * np.pi * 30000 * t / 96000)
    cases = [("outputFilter as-coded", hp_resampler(x, 96000, 0.5, "as-coded")),
             ("outputFilter state   ", hp_resampler(x, 96000, 0.5, "state")),
             ("inputFilter (#if 0)  ", hp_resampler_input_filter(x, 96000, 0.5)),
             ("naive ZOH decimation ", zoh_naive(x, 0.5))]
    for name, y in cases:
        print(f"   {name}: 18 kHz 处的幅度 = {bin_mag(y, 48000, 18000):.4f}   (满幅输入 = 1.0)")

    print()
    print("=" * 78)
    print("B) 升采样 48k->96k，单音保真度（最小二乘拟合幅度 / 残差 rms）")
    print("=" * 78)
    for f in (1000.0, 5000.0, 10000.0):
        x = np.sin(2 * np.pi * f * t / 48000)
        rows = [("as-coded", hp_resampler(x, 48000, 2.0, "as-coded")),
                ("state", hp_resampler(x, 48000, 2.0, "state")),
                ("ZOH", zoh_naive(x, 2.0))]
        txt = "  ".join(f"{n}: amp={ls_amp_resid(y, 96000, f)[0]:.3f} resid={ls_amp_resid(y, 96000, f)[1]:.4f}"
                        for n, y in rows)
        print(f"   f={f:6.0f} Hz  {txt}")

    print()
    print("=" * 78)
    print("C) 代码事实")
    print("=" * 78)
    flt = OutputFilterBank(1.0 / (48000 * 2))
    flt.set_freq(24000.0)
    print(f"   calcH0() = -sum Re(r/p) = {flt.calc_h0():.4f}   （定义于源码，但从未被调用）")
    flt.set_delta(1.0 / (48000 / 2))
    print(f"   |Aplus| = {np.abs(flt.Aplus)}   （恒等于 1：Gcalc 只有相位、无衰减记忆）")


if __name__ == "__main__":
    main()
