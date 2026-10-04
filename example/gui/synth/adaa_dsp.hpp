#pragma once

// ============================================================
// ADAA（反导数抗混叠）振荡器 / waveshaper —— GUI 例子自带的 DSP
// ============================================================
// 算法见论文：
//   [1] L. Gabrielli, S. D'Angelo, P. P. La Pastina, S. Squartini,
//       "Antiderivative Antialiasing for Arbitrary Waveform Generation",
//       IEEE/ACM TASLP 30:2743-2753, 2022.（振荡器：DPW ≡ AA-FIR、AA-IIR）
//   [2] P. P. La Pastina, S. D'Angelo, L. Gabrielli,
//       "Arbitrary-Order IIR Antiderivative Antialiasing", DAFx20in21.（静态非线性：AA-IIR）
// 本文件是 `qwqdsp/labs/adaa_iir`（Python 参考实现，21 项交叉验证全过）的 C++ 移植：
// 数学逐式对应，见各函数注释与 `temp/probe_adaa_dsp.cpp`（与 Python 逐点比对的探针）。
//
// 统一记号：非线性函数的输入记 x，逐样本差分 Δ = x_{n+1} - x_n；抗混叠核是若干复共轭
// 极点对 (B, β)（每个极点对的连续时域贡献为 2Re(B·e^{βt})）。递推：
//
//     ŷ_{n+1} = e^β·ŷ_n + 2B·I_n,      y = Re(ŷ)
//
// 其中加权平均积分（论文 [1] 式 (26) 的等价形式，但改用区间参数 u，避免 γ = β/Δ 在
// |Δ|→0 时的病态）：
//
//     I_n = ∫_0^1 f(x_n + uΔ)·e^{β(1-u)} du
//
// 分段线性的 f 时，把区间按断点切成「槽位」，每个槽位有闭式：
//
//     ∫_0^h (γ₀ + αz)e^{βz}dz = γ₀(e^{βh}-1)/β + α((βh-1)e^{βh}+1)/β²
//
// 代价 ∝ 每个采样区间跨过的段数：硬削波 / 波折叠 / 锯齿 / 三角只有 1~5 个折点，
// 因此每采样通常只有 1 个槽位（此时只需预计算好的 exp(β)，连 exp 都不用算）。

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <span>

#include "qwqdsp/filter/iir_design.hpp"

namespace adaa {

inline constexpr float kEps = 1e-9f;
inline constexpr double kTwoPi = 6.283185307179586476925286766559;
inline constexpr int kMaxPairs = 5;

// ------------------------------------------------------------
// 抗混叠滤波器：极点 / 留数
// ------------------------------------------------------------
// 极点本身就是「每采样」衰减率（与两篇论文的参考代码一致，见 labs/adaa_iir）。
//   AA-IIR-1 = butter(2, 2π·0.45)                      —— 固化
//   AA-IIR-2 = 10 阶椭圆 rs=80dB（库里的 IIRDesign，偶数阶修正）—— 现场设计
// ⚠ AA-IIR-2 的直接项 A0 = 1e-3（= 阻带增益）被忽略 —— 论文参考代码同样丢掉了它
//   （`[r,p,k] = residue(b,a)` 不用 k）。代价：通带增益 1-1e-3，即 -0.0087 dB。
struct PolePair {
    std::complex<float> residue;
    std::complex<float> pole;
};

struct Filter {
    int n_pairs{};
    std::array<PolePair, kMaxPairs> pairs{};
};

/// AA-IIR-1：2 阶 Butterworth，Fc = 0.45·Fs。
[[nodiscard]] inline Filter MakeButter2() {
    Filter f;
    f.n_pairs = 1;
    f.pairs[0] = {{0.0f, -1.9992973222f}, {-1.9992973222f, 1.9992973222f}};
    return f;
}

/// AA-IIR-2（旧）：10 阶 Chebyshev II，rs = 60 dB，阻带边沿 0.61·Fs。
///
/// ⚠ 保留只为对照：它的阻带边沿在 **Nyquist 之上**（0.61·fs），于是在
/// 0.5–0.61·fs 之间几乎没有衰减 —— 而折叠回基带的正是高于 Nyquist 的那部分，
/// 所以它的抗混叠在这些频率上只有十几 dB。新的 `MakeElliptic10` 就是修这个。
[[nodiscard]] inline Filter MakeCheby2() {
    Filter f;
    f.n_pairs = 5;
    f.pairs[0] = {{2.6281951982e-01f, -4.1698263224e-01f}, {-2.9931827180e-01f, 2.9476407233f}};
    f.pairs[1] = {{2.3416143280f, 1.0489688371f}, {-9.7441783728e-01f, 2.9828621117f}};
    f.pairs[2] = {{-2.6343139669f, 7.5174544707f}, {-1.8900320239f, 2.9479695675f}};
    f.pairs[3] = {{-1.6677182592e+01f, -6.7624788765f}, {-3.1558039886f, 2.5080127194f}};
    f.pairs[4] = {{1.6696365333e+01f, -2.3511355358e+01f}, {-4.3778067068f, 1.0814910500f}};
    return f;
}

/// AA-IIR-2（新）：**10 阶椭圆低通**，rs = 80 dB、偶数阶修正，用库里的
/// `qwqdsp_filter::IIRDesign::Elliptic` 现场设计，再自己算部分分式留数。
///
/// 两点关键：
/// 1. **偶数阶修正**（`even_order_modify = true`）。偶数阶椭圆在高频端的增益不为零
///    （|H(j∞)| = 阻带电平 d ≠ 0），对应时域里一个 δ（直接项），而 AA-IIR 只处理
///    「极点对」；库里的修正在 w² 上做 Möbius 变换把最低反射零点映到 0、最高传输零点
///    映到无穷远，于是严格真（|H(∞)| = 0）、DC 增益仍为 1、两端等波纹深度不变
///    （代价是阻带边沿略微外移）。scipy 没有这个修正。
/// 2. **阻带边沿必须落在 Nyquist(0.5·fs) 以下**：折叠回基带的是高于 Nyquist 的分量，
///    所以 `stop_normalized` 取 0.49 —— 整个折叠频带都被压到 -80 dB 以下。原型里先扫出
///    「首次到 -rs dB」的频率，再把原型整体缩放到那里（通带边沿随之落在 ~0.415·fs）。
///
/// 实测（labs/adaa_iir 的削波 SNR，88 键）：平均 68.0 → **83.5 dB**，最差音 33.1 → 49.3 dB；
/// 0.5·fs 处 -21 → -86 dB。代价与旧档相同（都是 5 对极点）。
[[nodiscard]] inline Filter MakeElliptic10(double stop_normalized = 0.49,
                                           double rs_db = 80.0, double rp_db = 0.1) {
    using qwqdsp_filter::IIRDesign;
    using Zpk = IIRDesign::ZPK;
    using C = std::complex<double>;
    constexpr size_t kPairs = 5;      // 10 阶 = 5 对共轭极点

    std::array<Zpk, kPairs> zpk{};
    if (!IIRDesign::Elliptic(std::span{zpk}, kPairs, rp_db, rs_db, /*even_order_modify=*/true)) {
        return MakeCheby2();          // 规格不合法（rs <= rp 之类）时退回旧档
    }
    const auto section = [](const Zpk& s, C x) {
        const C num = s.z.has_value() ? (x - *s.z) * (x - std::conj(*s.z)) : C{1.0, 0.0};
        return s.k * num / ((x - s.p) * (x - std::conj(s.p)));
    };
    // 原型（通带边沿 = 1 rad/s）里首次到 -rs dB 的频率
    double w_stop = 0.0;
    for (double w = 1.0; w < 4.0; w += 1e-4) {
        C h{1.0, 0.0};
        for (const auto& s : zpk) { h *= section(s, C{0.0, w}); }
        if (20.0 * std::log10(std::abs(h)) <= -rs_db) { w_stop = w; break; }
    }
    if (w_stop <= 0.0) { return MakeCheby2(); }
    const double scale = kTwoPi * stop_normalized / w_stop;

    Filter out;
    out.n_pairs = static_cast<int>(kPairs);
    for (size_t j = 0; j < kPairs; ++j) {
        const C p = zpk[j].p;
        // 留数 = k_j·N_j(p)/(p - conj(p)) · Π_{i≠j} H_i(p)
        C r = zpk[j].k
            * (zpk[j].z.has_value() ? (p - *zpk[j].z) * (p + *zpk[j].z) : C{1.0, 0.0})
            / (p - std::conj(p));
        for (size_t i = 0; i < kPairs; ++i) {
            if (i != j) { r *= section(zpk[i], p); }
        }
        out.pairs[j].residue = static_cast<std::complex<float>>(scale * r);
        out.pairs[j].pole = static_cast<std::complex<float>>(scale * p);
    }
    return out;
}

// ------------------------------------------------------------
// 单槽位的解析积分
// ------------------------------------------------------------

/**
 * 槽位 [a,b] 上的 ∫ f(ξ)·e^{β(x1-ξ)/Δ} dξ。
 *
 * @param a,b      槽位在 ξ 轴上的左右端（a < b，落在 [min(x0,x1), max(x0,x1)] 内）
 * @param x0,x1    采样区间端点（x1 = x0 + Δ）
 * @param m,q,off  槽位内的直线 f(ξ) = m·(ξ - off) + q
 * @param beta     极点
 * @param exp_beta 预计算的 exp(β)（整区间单槽时用它，省掉两次 exp）
 */
[[nodiscard]] inline std::complex<float> SlotIntegral(
    float a, float b, float x0, float x1,
    float m, float q, float off,
    std::complex<float> beta, std::complex<float> exp_beta) noexcept
{
    const float delta = x1 - x0;
    const float inv_delta = 1.0f / delta;
    const float ua = (a - x0) * inv_delta;
    const float ub = (b - x0) * inv_delta;
    const float u_hi = std::max(ua, ub);
    const float u_lo = std::min(ua, ub);
    const float w_a = 1.0f - u_hi;              ///< 槽位起点到区间末端的「采样距离」
    const float h = u_hi - u_lo;                ///< 槽位在采样尺度上的长度 ∈ [0,1]
    const float xi_start = x0 + u_hi * delta;   ///< 槽位起点的 ξ（沿路径）
    const float gamma0 = m * (xi_start - off) + q;
    const float alpha = -m * delta;

    // 整区间单槽（区间内没有断点）：h=1、w_a=0 → 用预计算值，省两次 exp（最常见的情形）
    const bool single = (h > 1.0f - 1e-6f) && (w_a < 1e-6f);
    const std::complex<float> e_h = single ? exp_beta : std::exp(beta * h);
    const std::complex<float> term = gamma0 * (e_h - 1.0f) / beta
                                   + alpha * ((beta * h - 1.0f) * e_h + 1.0f) / (beta * beta);
    return single ? term : term * std::exp(beta * w_a);
}

// ------------------------------------------------------------
// 非周期分段线性函数（waveshaper 用）
// ------------------------------------------------------------

/**
 * 分段线性非线性：段 j 覆盖 [x[j], x[j+1]]，f = m[j]·x + q[j]；域外由两条延长线覆盖
 * （钳位 / 折回 / 线性外推）。
 */
struct Pwl {
    static constexpr int kMaxSeg = 64;
    int n_seg{};
    std::array<float, kMaxSeg + 1> x{};
    std::array<float, kMaxSeg> m{};
    std::array<float, kMaxSeg> q{};
    float ml{}, ql{};   ///< x < x[0] 处的延长线
    float mr{}, qr{};   ///< x > x[n_seg] 处的延长线

    void AddSeg(float x0, float x1, float slope, float intercept) noexcept {
        x[n_seg] = x0;
        m[n_seg] = slope;
        q[n_seg] = intercept;
        ++n_seg;
        x[n_seg] = x1;
    }

    /// 预计算原函数用的断点累计量（构造完必须调用）。
    void Finish() noexcept {
        accum1_[0] = 0.0f;
        accum2_[0] = 0.0f;
        float a1 = 0.0f;
        float a2 = 0.0f;
        for (int i = 0; i < n_seg; ++i) {
            const float a = x[i];
            const float b = x[i + 1];
            const float d = b - a;
            a2 += a1 * d + m[i] * ((b * b * b - a * a * a) / 3.0f - a * a * d) * 0.5f
                + q[i] * ((b * b - a * a) * 0.5f - a * d);
            a1 += m[i] * (b * b - a * a) * 0.5f + q[i] * d;
            accum1_[i + 1] = a1;
            accum2_[i + 1] = a2;
        }
    }

    /// 取值（域外走延长线）。
    [[nodiscard]] float Value(float v) const noexcept {
        const int i = SegIndex(v);
        if (i == 0) { return ml * v + ql; }
        if (i == n_seg + 1) { return mr * v + qr; }
        return m[i - 1] * v + q[i - 1];
    }

    /// 扩展段索引：0 = 左延长线，1..n_seg = 内部段，n_seg+1 = 右延长线。
    [[nodiscard]] int SegIndex(float v) const noexcept {
        if (v < x[0]) { return 0; }
        if (v >= x[n_seg]) { return n_seg + 1; }
        int lo = 0;
        int hi = n_seg;
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            if (v < x[mid]) { hi = mid; } else { lo = mid; }
        }
        return lo + 1;
    }

    /// 一阶原函数（未归一化；AA-FIR 只用差分，常数项无关，但必须与 Raw2 配套）。
    [[nodiscard]] float Raw1(float v) const noexcept {
        const int i = SegIndex(v);
        if (i == 0) {
            const float lo = x[0];
            return -(ml * (lo * lo - v * v) * 0.5f + ql * (lo - v));
        }
        if (i == n_seg + 1) {
            const float hi = x[n_seg];
            return accum1_[n_seg] + mr * (v * v - hi * hi) * 0.5f + qr * (v - hi);
        }
        const float a = x[i - 1];
        return accum1_[i - 1] + m[i - 1] * (v * v - a * a) * 0.5f + q[i - 1] * (v - a);
    }

    /// 二阶原函数（满足 dRaw2/dv = Raw1 —— 这一条不成立时 AA-FIR-2 会全错）。
    [[nodiscard]] float Raw2(float v) const noexcept {
        const int i = SegIndex(v);
        if (i == 0) {
            const float lo = x[0];
            const float d = lo - v;
            return ml * 0.5f * (lo * lo * d - (lo * lo * lo - v * v * v) / 3.0f) + ql * d * d * 0.5f;
        }
        if (i == n_seg + 1) {
            const float hi = x[n_seg];
            const float d = v - hi;
            return accum2_[n_seg] + accum1_[n_seg] * d
                 + mr * ((v * v * v - hi * hi * hi) / 3.0f - hi * hi * d) * 0.5f
                 + qr * ((v * v - hi * hi) * 0.5f - hi * d);
        }
        const float a = x[i - 1];
        const float d = v - a;
        return accum2_[i - 1] + accum1_[i - 1] * d
             + m[i - 1] * ((v * v * v - a * a * a) / 3.0f - a * a * d) * 0.5f
             + q[i - 1] * ((v * v - a * a) * 0.5f - a * d);
    }

private:
    std::array<float, kMaxSeg + 1> accum1_{};
    std::array<float, kMaxSeg + 1> accum2_{};
};

/// 硬削波 f(x) = clip(x, ±threshold)（2 个折点）。
[[nodiscard]] inline Pwl MakeHardClip(float threshold = 1.0f) {
    Pwl f;
    f.AddSeg(-threshold, threshold, 1.0f, 0.0f);
    f.ml = 0.0f;
    f.ql = -threshold;
    f.mr = 0.0f;
    f.qr = threshold;
    f.Finish();
    return f;
}

/// 非对称波折叠（[DAFx25] 式 (15)）：|x| ≤ τ 直通，之外折回（2 个折点）。
[[nodiscard]] inline Pwl MakeWavefold(float tau = 0.7f) {
    Pwl f;
    f.AddSeg(-tau, tau, 1.0f, 0.0f);
    f.ml = -1.0f;
    f.ql = -2.0f * tau;
    f.mr = -1.0f;
    f.qr = 2.0f * tau;
    f.Finish();
    return f;
}

/// tanh(x/beta)·alpha 的分段线性拟合（域外按端值钳位）。
[[nodiscard]] inline Pwl MakeTanhFit(float alpha = 1.0f, float beta = 0.3f,
                                     float range = 8.0f, int segs = 48) {
    const auto fn = [&](float v) { return alpha * std::tanh(v / beta); };
    const float step = 2.0f * range / static_cast<float>(segs);
    Pwl f;
    for (int i = 0; i < segs; ++i) {
        const float a = -range + step * static_cast<float>(i);
        const float b = a + step;
        const float slope = (fn(b) - fn(a)) / step;
        f.AddSeg(a, b, slope, fn(a) - slope * a);
    }
    f.ml = 0.0f;
    f.ql = fn(-range);
    f.mr = 0.0f;
    f.qr = fn(range);
    f.Finish();
    return f;
}

// ------------------------------------------------------------
// 三次多项式非线性（「多项式」shape）
// ------------------------------------------------------------

/**
 * 三次多项式 f(x) = b1·x + b2·x² + b3·x³。
 *
 * GUI 的「多项式」shape 用它。两种用法：
 *   * 实数链路（trivial / AA-FIR / AA-IIR）：直接对实数输入作用，系数取 (1, d, d)（d = 驱动量）；
 *   * 解析链路（两个新增的过采样方法）：把 H(z) = z + z² + z³ 作用在复数解析信号上并取实部，
 *     归一化 `F = Re(H(d·z))/d` 是 Vicanek, *Complex Waveshapers* (2025) 的做法 —— 基波增益
 *     恒为 1、n 次谐波 ∝ d^(n−1)，d=1 时逐点等于 labs/adaa_iir 参考实现的 Re(z+z²+z³)。
 */
struct Cubic {
    float b1{1.0f};
    float b2{};
    float b3{};

    [[nodiscard]] float Value(float x) const noexcept {
        return b1 * x + b2 * x * x + b3 * x * x * x;
    }

    /// 一阶原函数（F(0) = 0）。
    [[nodiscard]] float Raw1(float x) const noexcept {
        const float x2 = x * x;
        return b1 * 0.5f * x2 + b2 * (1.0f / 3.0f) * x2 * x + b3 * 0.25f * x2 * x2;
    }

    /// 二阶原函数（F(0) = F'(0) = 0）。
    [[nodiscard]] float Raw2(float x) const noexcept {
        const float x2 = x * x;
        return b1 * (1.0f / 6.0f) * x2 * x + b2 * (1.0f / 12.0f) * x2 * x2
             + b3 * (1.0f / 20.0f) * x2 * x2 * x;
    }

    /**
     * AA-IIR 的加权积分 I = ∫_0^1 f(x0+uΔ)·e^{β(1-u)} du。
     *
     * 换元 v = 1−u 后 I = ∫_0^1 f(x1−vΔ)·e^{βv} dv；被积函数是 v 的三次多项式，
     * 用 J_n = ∫_0^1 v^n e^{βv} dv = (e^β − n·J_{n−1})/β（J_0 = (e^β−1)/β）闭式求和。
     */
    [[nodiscard]] std::complex<float> Forcing(
        float x0, float x1, std::complex<float> beta, std::complex<float> exp_beta) const noexcept
    {
        const float m = x0 - x1;                     // = −Δ
        const float c0 = b1 * x1 + b2 * x1 * x1 + b3 * x1 * x1 * x1;
        const float c1 = m * (b1 + 2.0f * b2 * x1 + 3.0f * b3 * x1 * x1);
        const float c2 = m * m * (b2 + 3.0f * b3 * x1);
        const float c3 = m * m * m * b3;
        const std::complex<float> j0 = (exp_beta - 1.0f) / beta;
        const std::complex<float> j1 = (exp_beta - j0) / beta;
        const std::complex<float> j2 = (exp_beta - 2.0f * j1) / beta;
        const std::complex<float> j3 = (exp_beta - 3.0f * j2) / beta;
        return c0 * j0 + c1 * j1 + c2 * j2 + c3 * j3;
    }
};

// ------------------------------------------------------------
// 周期分段线性波形（振荡器用）
// ------------------------------------------------------------

/**
 * 周期波形：段 j 覆盖 [x[j], x[j+1]]（周期 T = x[n_seg]），周期内坐标下 f = m[j]·x + q[j]。
 * 相位由调用方保持**单调递增**，跨周期通过「段号 + 周期号」还原（与论文把时间/相位作为
 * 非线性输入的做法一致：不是先回绕再插值）。
 */
struct PeriodicWave {
    static constexpr int kMaxSeg = 8;
    int n_seg{};
    float T{1.0f};
    std::array<float, kMaxSeg + 1> x{};
    std::array<float, kMaxSeg> m{};
    std::array<float, kMaxSeg> q{};
    std::array<float, kMaxSeg + 1> accum1{};
    std::array<float, kMaxSeg + 1> accum2{};

    void AddSeg(float x0, float x1, float slope, float intercept) noexcept {
        x[n_seg] = x0;
        m[n_seg] = slope;
        q[n_seg] = intercept;
        ++n_seg;
        x[n_seg] = x1;
        T = x1;
    }

    /// 预计算断点累计量与每周期漂移（构造完必须调用）。
    void Finish() noexcept {
        accum1[0] = 0.0f;
        accum2[0] = 0.0f;
        float a1 = 0.0f;
        float a2 = 0.0f;
        for (int i = 0; i < n_seg; ++i) {
            const float a = x[i];
            const float b = x[i + 1];
            const float d = b - a;
            a2 += a1 * d + m[i] * ((b * b * b - a * a * a) / 3.0f - a * a * d) * 0.5f
                + q[i] * ((b * b - a * a) * 0.5f - a * d);
            a1 += m[i] * (b * b - a * a) * 0.5f + q[i] * d;
            accum1[i + 1] = a1;
            accum2[i + 1] = a2;
        }
        drift1 = a1;
        drift2 = a2;
    }

    /// 取值（相位可任意，自动按周期还原）。
    /// ⚠ 表内直线存的是**周期内坐标**：f(r) = m[j]·r + q[j]（r = 相位 - 周期号·T）。
    [[nodiscard]] float Value(float phi) const noexcept {
        const float p = std::floor(phi / T);
        const float r = phi - p * T;
        const int j = LocalSeg(r);
        return m[j] * r + q[j];
    }

    /// 一阶原函数（含跨周期漂移）。``p_ref`` 是参考周期号：返回值相差一个常数，
    /// 差分不受影响，但可以把数值量级压住（见 Antideriv2 的说明）。
    [[nodiscard]] float Antideriv1(float phi, int p_ref = 0) const noexcept {
        const int p = static_cast<int>(std::floor(phi / T));
        const float r = phi - static_cast<float>(p) * T;
        const int j = LocalSeg(r);
        return static_cast<float>(p - p_ref) * drift1 + accum1[j]
             + m[j] * (r * r - x[j] * x[j]) * 0.5f + q[j] * (r - x[j]);
    }

    /// 二阶原函数（含跨周期漂移；满足 dAntideriv2/dφ = Antideriv1）。
    ///
    /// ⚠ ``p_ref`` 必须传：F₂ 里含 ``p·drift2`` 项，相位不回绕时 p 会一直涨，
    /// float32 下差分精度会随时间流失（长跑会出现杂音）。传参考周期号后量级被压住，
    /// 而 AA-FIR-2 用的是差分，常数项自动抵消。
    [[nodiscard]] float Antideriv2(float phi, int p_ref = 0) const noexcept {
        const int p = static_cast<int>(std::floor(phi / T));
        const float r = phi - static_cast<float>(p) * T;
        const int j = LocalSeg(r);
        const float a = x[j];
        const float d = r - a;
        const float loc2 = accum2[j] + accum1[j] * d
                         + m[j] * ((r * r * r - a * a * a) / 3.0f - a * a * d) * 0.5f
                         + q[j] * ((r * r - a * a) * 0.5f - a * d);
        const float dp = static_cast<float>(p - p_ref);
        // p(p-1) - p_ref(p_ref-1) = (p-p_ref)(p+p_ref-1)，用这个形式避免大数相减
        const float drift_quad = (drift1 * T) * 0.5f * dp
            * static_cast<float>(p + p_ref - 1);
        return dp * drift2 + drift_quad + loc2 + static_cast<float>(p) * drift1 * r;
    }

    /// 周期内段号（按真实断点查找，支持非均匀分段）。
    [[nodiscard]] int LocalSeg(float r) const noexcept {
        int lo = 0;
        int hi = n_seg;
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            if (r < x[mid]) { hi = mid; } else { lo = mid; }
        }
        return lo;
    }

    float drift1{};
    float drift2{};
};

/// 锯齿波 f(φ) = 2φ - 1（1 段）。
[[nodiscard]] inline PeriodicWave MakeSaw() {
    PeriodicWave w;
    w.AddSeg(0.0f, 1.0f, 2.0f, -1.0f);
    w.Finish();
    return w;
}

/// 三角波（2 段）。
[[nodiscard]] inline PeriodicWave MakeTriangle() {
    PeriodicWave w;
    w.AddSeg(0.0f, 0.5f, 4.0f, -1.0f);
    w.AddSeg(0.5f, 1.0f, -4.0f, 3.0f);
    w.Finish();
    return w;
}

/// 梯形「方波」：上升/下降沿横跨 1/edge_div 个周期，边沿**跨在周期边界上**
/// （f(0)=0、f(1)=0 连续）。edge_div 越大越接近理想方波。
[[nodiscard]] inline PeriodicWave MakeSquare(float edge_div = 128.0f) {
    const float e = 1.0f / edge_div;
    const float s = 2.0f / e;
    PeriodicWave w;
    w.AddSeg(0.0f, 0.5f * e, s, 0.0f);                       // 0 → 1
    w.AddSeg(0.5f * e, 0.5f - 0.5f * e, 0.0f, 1.0f);         // 平顶
    w.AddSeg(0.5f - 0.5f * e, 0.5f + 0.5f * e, -s, 1.0f / e); // 1 → -1
    w.AddSeg(0.5f + 0.5f * e, 1.0f - 0.5f * e, 0.0f, -1.0f); // 平底
    w.AddSeg(1.0f - 0.5f * e, 1.0f, s, -2.0f / e);           // -1 → 0
    w.Finish();
    return w;
}

// ------------------------------------------------------------
// 加权平均积分 I（两条路线：非周期 Pwl / 周期波形）
// ------------------------------------------------------------

/// 非周期分段线性 f：I = ∫_0^1 f(x0+uΔ)e^{β(1-u)}du。
[[nodiscard]] inline std::complex<float> ForcingPwl(
    float x0, float x1, std::complex<float> beta, std::complex<float> exp_beta,
    const Pwl& f) noexcept
{
    const float delta = x1 - x0;
    if (std::abs(delta) < kEps) {
        // Δ=0（输入停滞）：f 在区间上是常数 → 取极限 I = f(x)·(e^β-1)/β
        // （与论文参考代码的 tol 分支同值）
        return f.Value(x0) * (exp_beta - 1.0f) / beta;
    }
    const float lo = std::min(x0, x1);
    const float hi = std::max(x0, x1);
    const int i_lo = f.SegIndex(lo);
    const int i_hi = f.SegIndex(hi);
    std::complex<float> acc{};
    for (int i = i_lo; i <= i_hi; ++i) {
        const float left = (i == 0) ? lo : std::max(lo, f.x[i - 1]);
        const float right = (i == f.n_seg + 1) ? hi : std::min(hi, f.x[i]);
        if (right <= left) { continue; }
        float m = 0.0f;
        float q = 0.0f;
        // ⚠ Pwl 的直线是**绝对坐标**下的 f(ξ) = m·ξ + q，所以槽位参考点恒为 0
        //   （周期波形才需要按周期号偏移，见 ForcingPeriodic）。
        if (i == 0) { m = f.ml; q = f.ql; }
        else if (i == f.n_seg + 1) { m = f.mr; q = f.qr; }
        else { m = f.m[i - 1]; q = f.q[i - 1]; }
        acc += SlotIntegral(left, right, x0, x1, m, q, 0.0f, beta, exp_beta);
    }
    return acc;
}

/// 周期分段线性波形：同上，但沿**真实断点**（支持非均匀分段）逐段推进，跨周期由周期号
/// 还原局部坐标。
[[nodiscard]] inline std::complex<float> ForcingPeriodic(
    float x0, float x1, std::complex<float> beta, std::complex<float> exp_beta,
    const PeriodicWave& w) noexcept
{
    const float delta = x1 - x0;
    if (std::abs(delta) < kEps) {
        return w.Value(x0) * (exp_beta - 1.0f) / beta;
    }
    float cur = x0;
    int j = w.LocalSeg(x0 - std::floor(x0 / w.T) * w.T);
    int p = static_cast<int>(std::floor(x0 / w.T));
    std::complex<float> acc{};
    // 每次推进一个断点；Δ < 一个周期时最多 n_seg + 1 次
    for (int guard = 0; guard <= w.n_seg + 1; ++guard) {
        const float seg_right = w.x[j + 1] + static_cast<float>(p) * w.T;
        const float right = std::min(x1, seg_right);
        if (right > cur) {
            acc += SlotIntegral(cur, right, x0, x1, w.m[j], w.q[j],
                                static_cast<float>(p) * w.T, beta, exp_beta);
        }
        if (seg_right >= x1) { break; }
        cur = seg_right;
        ++j;
        if (j >= w.n_seg) { j = 0; ++p; }
    }
    return acc;
}

// ------------------------------------------------------------
// 方法枚举
// ------------------------------------------------------------

enum class Method {
    Trivial = 0,
    AaFir1,     ///< 一阶 ADAA（矩形核）＝ DPW-2
    AaFir2,     ///< 二阶 ADAA（三角核）＝ DPW-3
    AaIIR1,     ///< AA-IIR，2 阶 Butterworth
    AaIIR2,     ///< AA-IIR，10 阶 Chebyshev II
    NumMethods
};

inline constexpr const char* kMethodNames[] = {
    "trivial",
    "AA-FIR-1 / DPW-2",
    "AA-FIR-2 / DPW-3",
    "AA-IIR-1 (butter2)",
    "AA-IIR-2 (cheby2-10)",
};

// ------------------------------------------------------------
// AA-IIR 状态（并行结构：每个极点对一条一阶复递推）
// ------------------------------------------------------------

class AaIIRState {
public:
    void SetFilter(const Filter& f) noexcept {
        filter_ = f;
        for (int i = 0; i < f.n_pairs; ++i) {
            exp_[i] = std::exp(f.pairs[i].pole);
        }
        Reset();
    }

    void Reset() noexcept {
        y_hat_.fill({});
    }

    [[nodiscard]] const Filter& GetFilter() const noexcept { return filter_; }
    [[nodiscard]] std::complex<float> Exp(int i) const noexcept { return exp_[i]; }

    /// ŷ ← e^β·ŷ + forcing，返回新的 ŷ。
    [[nodiscard]] std::complex<float> Step(int i, std::complex<float> forcing) noexcept {
        y_hat_[i] = exp_[i] * y_hat_[i] + forcing;
        return y_hat_[i];
    }

private:
    Filter filter_{};
    std::array<std::complex<float>, kMaxPairs> y_hat_{};
    std::array<std::complex<float>, kMaxPairs> exp_{};
};

// ------------------------------------------------------------
// waveshaper：正弦（或任意信号）→ 分段线性非线性
// ------------------------------------------------------------

/**
 * 抗混叠 waveshaper。`Process(x)` 推入一个输入样本，返回一个输出样本。
 * 状态与论文参考代码一致：首样本按「前一个输入 = 0」处理。
 */
class Shaper {
public:
    void Reset() noexcept {
        iir_.Reset();
        prev_ = 0.0f;
        last_s_ = 0.0f;
        has_s_ = false;
        have_prev_ = false;
    }

    /// 分段线性 shape（默认）。
    void SetShape(const Pwl* f) noexcept { shape_ = f; cubic_ = nullptr; }

    /// 三次多项式 shape（「多项式」，见 Cubic）。
    void SetCubic(const Cubic* c) noexcept { cubic_ = c; shape_ = nullptr; }

    void SetMethod(Method m, const Filter& iir_filter) noexcept {
        method_ = m;
        iir_.SetFilter(iir_filter);
    }

    [[nodiscard]] float Process(float x) noexcept {
        float y = 0.0f;
        switch (method_) {
            case Method::Trivial:
                y = Value(x);
                break;
            case Method::AaFir1:
                y = Fir1(x);
                break;
            case Method::AaFir2:
                y = Fir2(x);
                break;
            case Method::AaIIR1:
            case Method::AaIIR2: {
                std::complex<float> acc{};
                for (int i = 0; i < iir_.GetFilter().n_pairs; ++i) {
                    const auto beta = iir_.GetFilter().pairs[i].pole;
                    const auto I = Forcing(prev_, x, beta, iir_.Exp(i));
                    acc += iir_.Step(i, 2.0f * iir_.GetFilter().pairs[i].residue * I);
                }
                y = acc.real();
                break;
            }
            default:
                break;
        }
        prev2_ = prev_;     ///< AA-FIR-2 的分母要 x_n - x_{n-2}（首样本按 0 处理，与参考实现一致）
        prev_ = x;
        have_prev_ = true;
        return y;
    }

private:
    // ----- shape 取值 / 原函数 / 加权积分（Pwl 与 Cubic 二选一）-----

    [[nodiscard]] float Value(float x) const noexcept {
        return (cubic_ != nullptr) ? cubic_->Value(x) : shape_->Value(x);
    }

    [[nodiscard]] float Raw1(float x) const noexcept {
        return (cubic_ != nullptr) ? cubic_->Raw1(x) : shape_->Raw1(x);
    }

    [[nodiscard]] float Raw2(float x) const noexcept {
        return (cubic_ != nullptr) ? cubic_->Raw2(x) : shape_->Raw2(x);
    }

    [[nodiscard]] std::complex<float> Forcing(float x0, float x1, std::complex<float> beta,
                                              std::complex<float> exp_beta) const noexcept {
        return (cubic_ != nullptr)
            ? cubic_->Forcing(x0, x1, beta, exp_beta)
            : ForcingPwl(x0, x1, beta, exp_beta, *shape_);
    }

    /// 一阶 ADAA：y = (F₁(x_n)-F₁(x_{n-1}))/(x_n-x_{n-1})，退化时取中点函数值。
    [[nodiscard]] float Fir1(float x) const noexcept {
        if (!have_prev_) { return Value(x); }
        const float d = x - prev_;
        if (std::abs(d) < kEps) { return Value(0.5f * (x + prev_)); }
        return (Raw1(x) - Raw1(prev_)) / d;
    }

    /// 二阶 ADAA（三角核）：
    ///   s_i = (F₂(x_i)-F₂(x_{i-1}))/(x_i-x_{i-1}),  y = 2·(s_n - s_{n-1})/(x_n-x_{n-2})
    ///
    /// ⚠ 两处都随 float32 精度做了退化保护（`labs` 的 Python 是 float64，看不出这个坑）：
    ///   * |Δ| 很小时 F₂ 的差商会抵消爆炸 → 改用极限值 s ≈ F₁(中点)（中值定理）；
    ///   * 分母 x_n - x_{n-2} 很小时退回首一阶结果。
    /// 阈值取相对量（信号幅度的 1e-3），否则正弦在削波饱和区（相邻样本几乎相等）会冒出
    /// 幅值上百的尖峰。
    [[nodiscard]] float Fir2(float x) noexcept {
        if (!have_prev_) {
            last_s_ = 0.0f;
            has_s_ = false;
            return Value(x);
        }
        const float scale = std::max({1.0f, std::abs(x), std::abs(prev_)});
        const float thr = 1e-3f * scale;
        const float d_now = x - prev_;
        const float s_now = (std::abs(d_now) < thr)
            ? Raw1(0.5f * (x + prev_))                      // Δ→0 的极限
            : (Raw2(x) - Raw2(prev_)) / d_now;
        float y = 0.0f;
        if (!has_s_) {
            y = Fir1(x);                        // 还没有 s_{n-1}
        }
        else {
            const float d_total = x - prev2_;
            y = (std::abs(d_total) < thr) ? Fir1(x)
                                          : 2.0f * (s_now - last_s_) / d_total;
        }
        last_s_ = s_now;
        has_s_ = true;
        return y;
    }

    const Pwl* shape_{};
    const Cubic* cubic_{};
    Method method_{Method::Trivial};
    AaIIRState iir_{};
    float prev_{};
    float prev2_{};
    float last_s_{};    ///< 上一时刻的 s（AA-FIR-2 用）
    bool has_s_{};
    bool have_prev_{};
};

// ------------------------------------------------------------
// 振荡器：相位（单调递增）→ 周期分段线性波形
// ------------------------------------------------------------

/**
 * 抗混叠振荡器。相位保持**单调递增**（Δ = f0/fs 恒定），跨周期由波形类还原 ——
 * 与论文「把时间/相位本身作为非线性输入」的做法一致。polyBLEP 档不在这里（见例子主文件，
 * 用库里的 `qwqdsp_oscillator::PolyBlep`）。
 */
class Oscillator {
public:
    void Reset() noexcept {
        iir_.Reset();
        phase_ = 0.0f;
        last_s_ = 0.0f;
        has_s_ = false;
    }

    void SetWave(const PeriodicWave* w) noexcept { wave_ = w; }

    void SetMethod(Method m, const Filter& iir_filter) noexcept {
        method_ = m;
        iir_.SetFilter(iir_filter);
    }

    void SetDelta(float delta) noexcept { delta_ = delta; }

    /// 当前（未回绕的）相位，供显示用。
    [[nodiscard]] double Phase() const noexcept { return phase_; }

    /// 产生下一个样本。
    ///
    /// 相位用 double 累加：float32 累加在长跑时会漂（每秒 ~1e-3 个周期量级），
    /// 漂移虽听不出来，但会让「跳变落在哪个采样区间」越来越偏，影响抗混叠修正的位置。
    [[nodiscard]] float Process() noexcept {
        // 相位按**整周期**重基准：相位不回绕时 float32 的分辨率随相位线性变差
        // （相位 1e3 时约 6e-5 个周期），跳变位置会抖、抗混叠残留被抬到 -50 dB 量级。
        // 减整周期是精确操作，且 forcing / AA-FIR 只用「相对周期差」，不影响结果。
        const double period = wave_->T;
        if (phase_ > 4.0 * period) {
            phase_ -= std::floor(phase_ / period) * period;
        }  // REBASE
        const auto x_prev = static_cast<float>(phase_);
        phase_ += static_cast<double>(delta_);
        const float x = static_cast<float>(phase_);
        const int p_ref = static_cast<int>(std::floor(x_prev / wave_->T));
        float y = 0.0f;
        switch (method_) {
            case Method::Trivial:
                y = wave_->Value(x);
                break;
            case Method::AaFir1: {
                const float d = x - x_prev;
                y = (std::abs(d) < kEps)
                    ? wave_->Value(x)
                    : (wave_->Antideriv1(x, p_ref) - wave_->Antideriv1(x_prev, p_ref)) / d;
                break;
            }
            case Method::AaFir2: {
                const float d_now = x - x_prev;
                const float s_now = (std::abs(d_now) < kEps)
                    ? wave_->Value(x)
                    : (wave_->Antideriv2(x, p_ref) - wave_->Antideriv2(x_prev, p_ref)) / d_now;
                y = has_s_ ? 2.0f * (s_now - last_s_) / (2.0f * delta_) : wave_->Value(x);
                last_s_ = s_now;
                has_s_ = true;
                break;
            }
            case Method::AaIIR1:
            case Method::AaIIR2: {
                std::complex<float> acc{};
                for (int i = 0; i < iir_.GetFilter().n_pairs; ++i) {
                    const auto beta = iir_.GetFilter().pairs[i].pole;
                    const auto I = ForcingPeriodic(x_prev, x, beta, iir_.Exp(i), *wave_);
                    acc += iir_.Step(i, 2.0f * iir_.GetFilter().pairs[i].residue * I);
                }
                y = acc.real();
                break;
            }
            default:
                break;
        }
        return y;
    }

private:
    const PeriodicWave* wave_{};
    Method method_{Method::Trivial};
    AaIIRState iir_{};
    double phase_{};
    float delta_{1.0f / 48000.0f};
    float last_s_{};
    bool has_s_{};
};

// ------------------------------------------------------------
// (A) FIR 多相解析过采样（L=8）
// ------------------------------------------------------------
// 与 labs/adaa_iir/polyphase_analytic_ovs.py 的 design() 逐式一致：
//   fs_in = 48 kHz、L = 8、fs_up = 384 kHz、512 抽头、通带 [0, 0.45·fs_in]。
//   h_lp = remez(通带 [0, f_p/2]、阻带自 [fs_in − 1.5·f_p])，f_p = 0.45·fs_in；
//   h_a[n] = 2L·h_lp[n]·exp(j·2π·(f_p/2)·(n−τ)/fs_up)，τ = (512−1)/2（相位基准必须取滤波器
//           中心，否则 Re(z) 是旋转过的信号而非带限插值本身）；
//   h_d = remez(通带 [0, f_p]、阻带自 [fs_in − f_p])，通带增益 1（抽取 ÷L）。
// 系数由 numpy/scipy 生成后原样贴入（生成脚本在 temp/，不入交付）。

static constexpr int kFirOvsL = 8;                                ///< 过采样倍率 L
static constexpr int kFirOvsTaps = 512;                           ///< 三条滤波器同长
static constexpr int kFirOvsPhaseTaps = kFirOvsTaps / kFirOvsL;   ///< 每相 64 抽头

static constexpr std::array<float, kFirOvsTaps> kFirOvsDecim = {
     1.839912155e-06f, -1.917877543e-06f, -1.777212396e-06f, -1.847425387e-06f,
    -1.876513312e-06f, -1.706161430e-06f, -1.250069288e-06f, -4.964647642e-07f,
     4.964242386e-07f,  1.601373751e-06f,  2.642788060e-06f,  3.417899552e-06f,
     3.736807171e-06f,  3.455073736e-06f,  2.515116354e-06f,  9.652616578e-07f,
    -1.024061750e-06f, -3.178968976e-06f, -5.146010952e-06f, -6.553112565e-06f,
    -7.061132024e-06f, -6.444488789e-06f, -4.629272422e-06f, -1.752430546e-06f,
     1.852575254e-06f,  5.679614603e-06f,  9.091740155e-06f,  1.146107196e-05f,
     1.223588682e-05f,  1.106472085e-05f,  7.879304458e-06f,  2.952145544e-06f,
    -3.116004768e-06f, -9.450442569e-06f, -1.502282387e-05f, -1.880169372e-05f,
    -1.993278868e-05f, -1.790120436e-05f, -1.266398271e-05f, -4.709187357e-06f,
     4.960124883e-06f,  1.493926040e-05f,  2.360806737e-05f,  2.938319911e-05f,
     3.098256130e-05f,  2.768022758e-05f,  1.947947225e-05f,  7.203256095e-06f,
    -7.568824696e-06f, -2.267209316e-05f, -3.566001809e-05f, -4.418420787e-05f,
    -4.638783623e-05f, -4.126554396e-05f, -2.891777653e-05f, -1.064248959e-05f,
     1.116184760e-05f,  3.328404003e-05f,  5.214937500e-05f,  6.437372636e-05f,
     6.733756941e-05f,  5.968864468e-05f,  4.168096238e-05f,  1.528198781e-05f,
    -1.599374642e-05f, -4.751846528e-05f, -7.420881182e-05f, -9.131443202e-05f,
    -9.522341732e-05f, -8.415065607e-05f, -5.858554350e-05f, -2.141139753e-05f,
     2.236507525e-05f,  6.624230933e-05f,  1.031586198e-04f,  1.265912409e-04f,
     1.316574316e-04f,  1.160410156e-04f,  8.057687125e-05f,  2.936757515e-05f,
    -3.061902806e-05f, -9.044906499e-05f, -1.405110921e-04f, -1.720164051e-04f,
    -1.784806613e-04f, -1.569471630e-04f, -1.087316659e-04f, -3.953444289e-05f,
     4.114905230e-05f,  1.212723127e-04f,  1.879881779e-04f,  2.296536326e-04f,
     2.377904976e-04f,  2.086743675e-04f,  1.442755575e-04f,  5.234983927e-05f,
    -5.439684147e-05f, -1.599937641e-04f, -2.475382972e-04f, -3.018371827e-04f,
    -3.119560880e-04f, -2.732615686e-04f, -1.885913871e-04f, -6.830503601e-05f,
     7.086576870e-05f,  2.080642450e-04f,  3.213652581e-04f,  3.912048921e-04f,
     4.036556846e-04f,  3.530144872e-04f,  2.432430463e-04f,  8.795799023e-05f,
    -9.112225087e-05f, -2.671212613e-04f, -4.119584866e-04f, -5.007404905e-04f,
    -5.159210863e-04f, -4.505442352e-04f, -3.100039988e-04f, -1.119420225e-04f,
     1.158089972e-04f,  3.390276365e-04f,  5.221524456e-04f,  6.338448339e-04f,
     6.522112083e-04f,  5.688349710e-04f,  3.909035674e-04f,  1.409827271e-04f,
    -1.456635725e-04f, -4.259228051e-04f, -6.552093055e-04f, -7.944374293e-04f,
    -8.165208661e-04f, -7.113383790e-04f, -4.882930951e-04f, -1.759218261e-04f,
     1.815458618e-04f,  5.303051534e-04f,  8.149461305e-04f,  9.871172844e-04f,
     1.013548425e-03f,  8.821255421e-04f,  6.049541825e-04f,  2.177590804e-04f,
    -2.244778824e-04f, -6.551550410e-04f, -1.005931509e-03f, -1.217406808e-03f,
    -1.248955553e-03f, -1.086119969e-03f, -7.442632065e-04f, -2.677125416e-04f,
     2.757072706e-04f,  8.041273798e-04f,  1.233790834e-03f,  1.492132359e-03f,
     1.529769063e-03f,  1.329456309e-03f,  9.104440171e-04f,  3.273118838e-04f,
    -3.368103243e-04f, -9.818569333e-04f, -1.505687800e-03f, -1.820026065e-03f,
    -1.865017299e-03f, -1.620053403e-03f, -1.108974672e-03f, -3.985510696e-04f,
     4.098498451e-04f,  1.194439073e-03f,  1.831090392e-03f,  2.212685967e-03f,
     2.266752142e-03f,  1.968534963e-03f,  1.347241047e-03f,  4.841304979e-04f,
    -4.976381792e-04f, -1.450225342e-03f, -2.223038720e-03f, -2.686170690e-03f,
    -2.751750668e-03f, -2.389770720e-03f, -1.635638714e-03f, -5.878740142e-04f,
     6.041753334e-04f,  1.761164935e-03f,  2.700298712e-03f,  3.263724996e-03f,
     3.344440291e-03f,  2.905539776e-03f,  1.989484388e-03f,  7.154474856e-04f,
    -7.354332607e-04f, -2.145201766e-03f, -3.291229057e-03f, -3.980703398e-03f,
    -4.082219673e-03f, -3.549410646e-03f, -2.432548493e-03f, -8.757068217e-04f,
     9.008074877e-04f,  2.630779429e-03f,  4.041118104e-03f,  4.893997361e-03f,
     5.025739649e-03f,  4.376257961e-03f,  3.004006412e-03f,  1.083367180e-03f,
    -1.116057233e-03f, -3.266019463e-03f, -5.027288479e-03f, -6.101660987e-03f,
    -6.280558566e-03f, -5.482555048e-03f, -3.773454501e-03f, -1.364856474e-03f,
     1.409796553e-03f,  4.139326106e-03f,  6.393517803e-03f,  7.788326886e-03f,
     8.048055615e-03f,  7.054875749e-03f,  4.877415391e-03f,  1.772789355e-03f,
    -1.839902173e-03f, -5.432438759e-03f, -8.440482149e-03f, -1.034712211e-02f,
    -1.076522084e-02f, -9.506257058e-03f, -6.624574058e-03f, -2.428820929e-03f,
     2.543398214e-03f,  7.586608274e-03f,  1.191805086e-02f,  1.478696691e-02f,
     1.558845604e-02f,  1.396628418e-02f,  9.889588364e-03f,  3.691115856e-03f,
    -3.940843758e-03f, -1.201886718e-02f, -1.935681281e-02f, -2.470383676e-02f,
    -2.689664659e-02f, -2.501150191e-02f, -1.849664894e-02f, -7.267552255e-03f,
     8.248691278e-03f,  2.712745766e-02f,  4.801984282e-02f,  6.927832322e-02f,
     8.912142886e-02f,  1.058189027e-01f,  1.178759352e-01f,  1.241948185e-01f,
     1.241948185e-01f,  1.178759352e-01f,  1.058189027e-01f,  8.912142886e-02f,
     6.927832322e-02f,  4.801984282e-02f,  2.712745766e-02f,  8.248691278e-03f,
    -7.267552255e-03f, -1.849664894e-02f, -2.501150191e-02f, -2.689664659e-02f,
    -2.470383676e-02f, -1.935681281e-02f, -1.201886718e-02f, -3.940843758e-03f,
     3.691115856e-03f,  9.889588364e-03f,  1.396628418e-02f,  1.558845604e-02f,
     1.478696691e-02f,  1.191805086e-02f,  7.586608274e-03f,  2.543398214e-03f,
    -2.428820929e-03f, -6.624574058e-03f, -9.506257058e-03f, -1.076522084e-02f,
    -1.034712211e-02f, -8.440482149e-03f, -5.432438759e-03f, -1.839902173e-03f,
     1.772789355e-03f,  4.877415391e-03f,  7.054875749e-03f,  8.048055615e-03f,
     7.788326886e-03f,  6.393517803e-03f,  4.139326106e-03f,  1.409796553e-03f,
    -1.364856474e-03f, -3.773454501e-03f, -5.482555048e-03f, -6.280558566e-03f,
    -6.101660987e-03f, -5.027288479e-03f, -3.266019463e-03f, -1.116057233e-03f,
     1.083367180e-03f,  3.004006412e-03f,  4.376257961e-03f,  5.025739649e-03f,
     4.893997361e-03f,  4.041118104e-03f,  2.630779429e-03f,  9.008074877e-04f,
    -8.757068217e-04f, -2.432548493e-03f, -3.549410646e-03f, -4.082219673e-03f,
    -3.980703398e-03f, -3.291229057e-03f, -2.145201766e-03f, -7.354332607e-04f,
     7.154474856e-04f,  1.989484388e-03f,  2.905539776e-03f,  3.344440291e-03f,
     3.263724996e-03f,  2.700298712e-03f,  1.761164935e-03f,  6.041753334e-04f,
    -5.878740142e-04f, -1.635638714e-03f, -2.389770720e-03f, -2.751750668e-03f,
    -2.686170690e-03f, -2.223038720e-03f, -1.450225342e-03f, -4.976381792e-04f,
     4.841304979e-04f,  1.347241047e-03f,  1.968534963e-03f,  2.266752142e-03f,
     2.212685967e-03f,  1.831090392e-03f,  1.194439073e-03f,  4.098498451e-04f,
    -3.985510696e-04f, -1.108974672e-03f, -1.620053403e-03f, -1.865017299e-03f,
    -1.820026065e-03f, -1.505687800e-03f, -9.818569333e-04f, -3.368103243e-04f,
     3.273118838e-04f,  9.104440171e-04f,  1.329456309e-03f,  1.529769063e-03f,
     1.492132359e-03f,  1.233790834e-03f,  8.041273798e-04f,  2.757072706e-04f,
    -2.677125416e-04f, -7.442632065e-04f, -1.086119969e-03f, -1.248955553e-03f,
    -1.217406808e-03f, -1.005931509e-03f, -6.551550410e-04f, -2.244778824e-04f,
     2.177590804e-04f,  6.049541825e-04f,  8.821255421e-04f,  1.013548425e-03f,
     9.871172844e-04f,  8.149461305e-04f,  5.303051534e-04f,  1.815458618e-04f,
    -1.759218261e-04f, -4.882930951e-04f, -7.113383790e-04f, -8.165208661e-04f,
    -7.944374293e-04f, -6.552093055e-04f, -4.259228051e-04f, -1.456635725e-04f,
     1.409827271e-04f,  3.909035674e-04f,  5.688349710e-04f,  6.522112083e-04f,
     6.338448339e-04f,  5.221524456e-04f,  3.390276365e-04f,  1.158089972e-04f,
    -1.119420225e-04f, -3.100039988e-04f, -4.505442352e-04f, -5.159210863e-04f,
    -5.007404905e-04f, -4.119584866e-04f, -2.671212613e-04f, -9.112225087e-05f,
     8.795799023e-05f,  2.432430463e-04f,  3.530144872e-04f,  4.036556846e-04f,
     3.912048921e-04f,  3.213652581e-04f,  2.080642450e-04f,  7.086576870e-05f,
    -6.830503601e-05f, -1.885913871e-04f, -2.732615686e-04f, -3.119560880e-04f,
    -3.018371827e-04f, -2.475382972e-04f, -1.599937641e-04f, -5.439684147e-05f,
     5.234983927e-05f,  1.442755575e-04f,  2.086743675e-04f,  2.377904976e-04f,
     2.296536326e-04f,  1.879881779e-04f,  1.212723127e-04f,  4.114905230e-05f,
    -3.953444289e-05f, -1.087316659e-04f, -1.569471630e-04f, -1.784806613e-04f,
    -1.720164051e-04f, -1.405110921e-04f, -9.044906499e-05f, -3.061902806e-05f,
     2.936757515e-05f,  8.057687125e-05f,  1.160410156e-04f,  1.316574316e-04f,
     1.265912409e-04f,  1.031586198e-04f,  6.624230933e-05f,  2.236507525e-05f,
    -2.141139753e-05f, -5.858554350e-05f, -8.415065607e-05f, -9.522341732e-05f,
    -9.131443202e-05f, -7.420881182e-05f, -4.751846528e-05f, -1.599374642e-05f,
     1.528198781e-05f,  4.168096238e-05f,  5.968864468e-05f,  6.733756941e-05f,
     6.437372636e-05f,  5.214937500e-05f,  3.328404003e-05f,  1.116184760e-05f,
    -1.064248959e-05f, -2.891777653e-05f, -4.126554396e-05f, -4.638783623e-05f,
    -4.418420787e-05f, -3.566001809e-05f, -2.267209316e-05f, -7.568824696e-06f,
     7.203256095e-06f,  1.947947225e-05f,  2.768022758e-05f,  3.098256130e-05f,
     2.938319911e-05f,  2.360806737e-05f,  1.493926040e-05f,  4.960124883e-06f,
    -4.709187357e-06f, -1.266398271e-05f, -1.790120436e-05f, -1.993278868e-05f,
    -1.880169372e-05f, -1.502282387e-05f, -9.450442569e-06f, -3.116004768e-06f,
     2.952145544e-06f,  7.879304458e-06f,  1.106472085e-05f,  1.223588682e-05f,
     1.146107196e-05f,  9.091740155e-06f,  5.679614603e-06f,  1.852575254e-06f,
    -1.752430546e-06f, -4.629272422e-06f, -6.444488789e-06f, -7.061132024e-06f,
    -6.553112565e-06f, -5.146010952e-06f, -3.178968976e-06f, -1.024061750e-06f,
     9.652616578e-07f,  2.515116354e-06f,  3.455073736e-06f,  3.736807171e-06f,
     3.417899552e-06f,  2.642788060e-06f,  1.601373751e-06f,  4.964242386e-07f,
    -4.964647642e-07f, -1.250069288e-06f, -1.706161430e-06f, -1.876513312e-06f,
    -1.847425387e-06f, -1.777212396e-06f, -1.917877543e-06f,  1.839912155e-06f,
};

static constexpr std::array<float, kFirOvsTaps> kFirOvsAnalyticRe = {
    -2.063270977e-05f, -1.135590674e-05f, -1.582847067e-05f, -1.945148534e-05f,
    -2.124740217e-05f, -2.055613213e-05f, -1.689818811e-05f, -1.043429202e-05f,
    -1.645505989e-06f,  8.235974222e-06f,  1.771789340e-05f,  2.491640770e-05f,
     2.840239803e-05f,  2.679815527e-05f,  1.968753090e-05f,  7.539146479e-06f,
    -8.163585888e-06f, -2.500275192e-05f, -3.992626659e-05f, -4.968207580e-05f,
    -5.139310728e-05f, -4.304493949e-05f, -2.406586218e-05f,  4.449093218e-06f,
     3.963709209e-05f,  7.723197585e-05f,  1.118762874e-04f,  1.382570813e-04f,
     1.516621155e-04f,  1.492055458e-04f,  1.302974888e-04f,  9.726060827e-05f,
     5.514643666e-05f,  1.132360253e-05f, -2.555619108e-05f, -4.683871151e-05f,
    -4.533473866e-05f, -1.673947057e-05f,  3.935114893e-05f,  1.188350269e-04f,
     2.133993433e-04f,  3.112729327e-04f,  3.990178014e-04f,  4.633559325e-04f,
     4.936708628e-04f,  4.838149576e-04f,  4.338479581e-04f,  3.504698021e-04f,
     2.468123279e-04f,  1.407235797e-04f,  5.245102531e-05f,  1.407227108e-06f,
     2.930851643e-06f,  6.530914091e-05f,  1.876885513e-04f,  3.593572843e-04f,
     5.604545260e-04f,  7.644764042e-04f,  9.418317040e-04f,  1.064576765e-03f,
     1.110958076e-03f,  1.069743731e-03f,  9.427926229e-04f,  7.461627909e-04f,
     5.085120997e-04f,  2.677468885e-04f,  6.542904503e-05f, -5.953715426e-05f,
    -7.763008424e-05f,  2.573300660e-05f,  2.464681794e-04f,  5.608914372e-04f,
     9.279153561e-04f,  1.294184456e-03f,  1.601783410e-03f,  1.797110624e-03f,
     1.840038363e-03f,  1.711306778e-03f,  1.417484005e-03f,  9.917171168e-04f,
     4.906718776e-04f, -1.303019951e-05f, -4.411209499e-04f, -7.221002112e-04f,
    -8.031472866e-04f, -6.602349509e-04f, -3.040885184e-04f,  2.188256652e-04f,
     8.311339392e-04f,  1.434683776e-03f,  1.924745720e-03f,  2.206400383e-03f,
     2.210373582e-03f,  1.906103271e-03f,  1.309227491e-03f,  4.824991774e-04f,
    -4.711654075e-04f, -1.422196442e-03f, -2.233489465e-03f, -2.781652735e-03f,
    -2.977895882e-03f, -2.784415634e-03f, -2.224053648e-03f, -1.380589515e-03f,
    -3.895221412e-04f,  5.802046814e-04f,  1.350939934e-03f,  1.764103520e-03f,
     1.706763394e-03f,  1.132532450e-03f,  7.334438223e-05f, -1.360740823e-03f,
    -2.993948667e-03f, -4.608763567e-03f, -5.977727075e-03f, -6.899156010e-03f,
    -7.230582526e-03f, -6.915190313e-03f, -5.995864798e-03f, -4.614550219e-03f,
    -2.995508374e-03f, -1.414717041e-03f, -1.590161882e-04f,  5.187697396e-04f,
     4.420378132e-04f, -4.570234465e-04f, -2.120278218e-03f, -4.365376318e-03f,
    -6.907453225e-03f, -9.398130823e-03f, -1.147707330e-02f, -1.282808404e-02f,
    -1.323169547e-02f, -1.260506581e-02f, -1.102279976e-02f, -8.713672163e-03f,
    -6.033067867e-03f, -3.413568686e-03f, -1.300860799e-03f, -8.383792654e-05f,
    -3.003474697e-05f, -1.236625590e-03f, -3.605747939e-03f, -6.849309905e-03f,
    -1.052404299e-02f, -1.409303262e-02f, -1.700514372e-02f, -1.878108842e-02f,
    -1.909246600e-02f, -1.782140481e-02f, -1.508985461e-02f, -1.125288498e-02f,
    -6.854603251e-03f, -2.552296478e-03f,  9.814365277e-04f,  3.163743503e-03f,
     3.598038304e-03f,  2.145953719e-03f, -1.037550550e-03f, -5.513529930e-03f,
    -1.061350038e-02f, -1.553271333e-02f, -1.944960272e-02f, -2.165377176e-02f,
    -2.166305927e-02f, -1.931043255e-02f, -1.478592662e-02f, -8.624389088e-03f,
    -1.638742893e-03f,  5.193935657e-03f,  1.087746836e-02f,  1.455201092e-02f,
     1.563554376e-02f,  1.392896566e-02f,  9.666159630e-03f,  3.498808777e-03f,
    -3.585091532e-03f, -1.039933157e-02f, -1.573541945e-02f, -1.855064939e-02f,
    -1.813939680e-02f, -1.426017576e-02f, -7.196219130e-03f,  2.262851202e-03f,
     1.291957062e-02f,  2.333555466e-02f,  3.204445140e-02f,  3.777953554e-02f,
     3.968154520e-02f,  3.745298143e-02f,  3.143255857e-02f,  2.257393248e-02f,
     1.232760195e-02f,  2.439350962e-03f, -5.307346597e-03f, -9.367777582e-03f,
    -8.689189402e-03f, -2.900204885e-03f,  7.596105859e-03f,  2.164493280e-02f,
     3.748037062e-02f,  5.296839731e-02f,  6.591949571e-02f,  7.442633806e-02f,
     7.717524113e-02f,  7.368215041e-02f,  6.441240501e-02f,  5.075995845e-02f,
     3.488153795e-02f,  1.940382009e-02f,  7.041710490e-03f,  1.818495143e-04f,
     4.933857906e-04f,  8.627415141e-03f,  2.405613607e-02f,  4.508423791e-02f,
     6.904063982e-02f,  9.263154428e-02f,  1.124104210e-01f,  1.253000037e-01f,
     1.290901554e-01f,  1.228343317e-01f,  1.070785627e-01f,  8.387774610e-02f,
     5.658364018e-02f,  2.942172490e-02f,  6.906946083e-03f, -6.825021641e-03f,
    -8.676850001e-03f,  2.859368679e-03f,  2.734791722e-02f,  6.231780981e-02f,
     1.034655970e-01f,  1.451441736e-01f,  1.810876843e-01f,  2.052860857e-01f,
     2.128963083e-01f,  2.010650266e-01f,  1.695425040e-01f,  1.209886839e-01f,
     6.090866312e-02f, -2.798599294e-03f, -6.064377614e-02f, -1.026194911e-01f,
    -1.194186167e-01f, -1.036505179e-01f, -5.090076777e-02f,  3.950683543e-02f,
     1.641548600e-01f,  3.156450740e-01f,  4.831844966e-01f,  6.535908324e-01f,
     8.126100580e-01f,  9.463979244e-01f,  1.042994069e+00f,  1.093615444e+00f,
     1.093615444e+00f,  1.042994069e+00f,  9.463979244e-01f,  8.126100580e-01f,
     6.535908324e-01f,  4.831844966e-01f,  3.156450740e-01f,  1.641548600e-01f,
     3.950683543e-02f, -5.090076777e-02f, -1.036505179e-01f, -1.194186167e-01f,
    -1.026194911e-01f, -6.064377614e-02f, -2.798599294e-03f,  6.090866312e-02f,
     1.209886839e-01f,  1.695425040e-01f,  2.010650266e-01f,  2.128963083e-01f,
     2.052860857e-01f,  1.810876843e-01f,  1.451441736e-01f,  1.034655970e-01f,
     6.231780981e-02f,  2.734791722e-02f,  2.859368679e-03f, -8.676850001e-03f,
    -6.825021641e-03f,  6.906946083e-03f,  2.942172490e-02f,  5.658364018e-02f,
     8.387774610e-02f,  1.070785627e-01f,  1.228343317e-01f,  1.290901554e-01f,
     1.253000037e-01f,  1.124104210e-01f,  9.263154428e-02f,  6.904063982e-02f,
     4.508423791e-02f,  2.405613607e-02f,  8.627415141e-03f,  4.933857906e-04f,
     1.818495143e-04f,  7.041710490e-03f,  1.940382009e-02f,  3.488153795e-02f,
     5.075995845e-02f,  6.441240501e-02f,  7.368215041e-02f,  7.717524113e-02f,
     7.442633806e-02f,  6.591949571e-02f,  5.296839731e-02f,  3.748037062e-02f,
     2.164493280e-02f,  7.596105859e-03f, -2.900204885e-03f, -8.689189402e-03f,
    -9.367777582e-03f, -5.307346597e-03f,  2.439350962e-03f,  1.232760195e-02f,
     2.257393248e-02f,  3.143255857e-02f,  3.745298143e-02f,  3.968154520e-02f,
     3.777953554e-02f,  3.204445140e-02f,  2.333555466e-02f,  1.291957062e-02f,
     2.262851202e-03f, -7.196219130e-03f, -1.426017576e-02f, -1.813939680e-02f,
    -1.855064939e-02f, -1.573541945e-02f, -1.039933157e-02f, -3.585091532e-03f,
     3.498808777e-03f,  9.666159630e-03f,  1.392896566e-02f,  1.563554376e-02f,
     1.455201092e-02f,  1.087746836e-02f,  5.193935657e-03f, -1.638742893e-03f,
    -8.624389088e-03f, -1.478592662e-02f, -1.931043255e-02f, -2.166305927e-02f,
    -2.165377176e-02f, -1.944960272e-02f, -1.553271333e-02f, -1.061350038e-02f,
    -5.513529930e-03f, -1.037550550e-03f,  2.145953719e-03f,  3.598038304e-03f,
     3.163743503e-03f,  9.814365277e-04f, -2.552296478e-03f, -6.854603251e-03f,
    -1.125288498e-02f, -1.508985461e-02f, -1.782140481e-02f, -1.909246600e-02f,
    -1.878108842e-02f, -1.700514372e-02f, -1.409303262e-02f, -1.052404299e-02f,
    -6.849309905e-03f, -3.605747939e-03f, -1.236625590e-03f, -3.003474697e-05f,
    -8.383792654e-05f, -1.300860799e-03f, -3.413568686e-03f, -6.033067867e-03f,
    -8.713672163e-03f, -1.102279976e-02f, -1.260506581e-02f, -1.323169547e-02f,
    -1.282808404e-02f, -1.147707330e-02f, -9.398130823e-03f, -6.907453225e-03f,
    -4.365376318e-03f, -2.120278218e-03f, -4.570234465e-04f,  4.420378132e-04f,
     5.187697396e-04f, -1.590161882e-04f, -1.414717041e-03f, -2.995508374e-03f,
    -4.614550219e-03f, -5.995864798e-03f, -6.915190313e-03f, -7.230582526e-03f,
    -6.899156010e-03f, -5.977727075e-03f, -4.608763567e-03f, -2.993948667e-03f,
    -1.360740823e-03f,  7.334438223e-05f,  1.132532450e-03f,  1.706763394e-03f,
     1.764103520e-03f,  1.350939934e-03f,  5.802046814e-04f, -3.895221412e-04f,
    -1.380589515e-03f, -2.224053648e-03f, -2.784415634e-03f, -2.977895882e-03f,
    -2.781652735e-03f, -2.233489465e-03f, -1.422196442e-03f, -4.711654075e-04f,
     4.824991774e-04f,  1.309227491e-03f,  1.906103271e-03f,  2.210373582e-03f,
     2.206400383e-03f,  1.924745720e-03f,  1.434683776e-03f,  8.311339392e-04f,
     2.188256652e-04f, -3.040885184e-04f, -6.602349509e-04f, -8.031472866e-04f,
    -7.221002112e-04f, -4.411209499e-04f, -1.303019951e-05f,  4.906718776e-04f,
     9.917171168e-04f,  1.417484005e-03f,  1.711306778e-03f,  1.840038363e-03f,
     1.797110624e-03f,  1.601783410e-03f,  1.294184456e-03f,  9.279153561e-04f,
     5.608914372e-04f,  2.464681794e-04f,  2.573300660e-05f, -7.763008424e-05f,
    -5.953715426e-05f,  6.542904503e-05f,  2.677468885e-04f,  5.085120997e-04f,
     7.461627909e-04f,  9.427926229e-04f,  1.069743731e-03f,  1.110958076e-03f,
     1.064576765e-03f,  9.418317040e-04f,  7.644764042e-04f,  5.604545260e-04f,
     3.593572843e-04f,  1.876885513e-04f,  6.530914091e-05f,  2.930851643e-06f,
     1.407227108e-06f,  5.245102531e-05f,  1.407235797e-04f,  2.468123279e-04f,
     3.504698021e-04f,  4.338479581e-04f,  4.838149576e-04f,  4.936708628e-04f,
     4.633559325e-04f,  3.990178014e-04f,  3.112729327e-04f,  2.133993433e-04f,
     1.188350269e-04f,  3.935114893e-05f, -1.673947057e-05f, -4.533473866e-05f,
    -4.683871151e-05f, -2.555619108e-05f,  1.132360253e-05f,  5.514643666e-05f,
     9.726060827e-05f,  1.302974888e-04f,  1.492055458e-04f,  1.516621155e-04f,
     1.382570813e-04f,  1.118762874e-04f,  7.723197585e-05f,  3.963709209e-05f,
     4.449093218e-06f, -2.406586218e-05f, -4.304493949e-05f, -5.139310728e-05f,
    -4.968207580e-05f, -3.992626659e-05f, -2.500275192e-05f, -8.163585888e-06f,
     7.539146479e-06f,  1.968753090e-05f,  2.679815527e-05f,  2.840239803e-05f,
     2.491640770e-05f,  1.771789340e-05f,  8.235974222e-06f, -1.645505989e-06f,
    -1.043429202e-05f, -1.689818811e-05f, -2.055613213e-05f, -2.124740217e-05f,
    -1.945148534e-05f, -1.582847067e-05f, -1.135590674e-05f, -2.063270977e-05f,
};

static constexpr std::array<float, kFirOvsTaps> kFirOvsAnalyticIm = {
     4.846057553e-05f,  1.736191028e-05f,  1.678942486e-05f,  1.442620473e-05f,
     1.056464084e-05f,  6.015902141e-06f,  1.832000651e-06f, -7.181999416e-07f,
    -4.121777932e-07f,  3.699234065e-06f,  1.209199242e-05f,  2.443191610e-05f,
     3.991071071e-05f,  5.665993882e-05f,  7.252468313e-05f,  8.510352198e-05f,
     9.215232957e-05f,  9.210482868e-05f,  8.441699810e-05f,  6.981265992e-05f,
     5.039378470e-05f,  2.937702979e-05f,  1.080931712e-05f, -1.114439836e-06f,
    -2.728250001e-06f,  8.373029650e-06f,  3.274141228e-05f,  6.874423496e-05f,
     1.124802908e-04f,  1.582638875e-04f,  1.992102755e-04f,  2.284384894e-04f,
     2.400598765e-04f,  2.304972247e-04f,  1.991525379e-04f,  1.491204152e-04f,
     8.685576674e-05f,  2.166811241e-05f, -3.566580243e-05f, -7.443673327e-05f,
    -8.594826909e-05f, -6.509917991e-05f, -1.175544358e-05f,  6.873235982e-05f,
     1.657789788e-04f,  2.647437785e-04f,  3.489785794e-04f,  4.022822776e-04f,
     4.117814918e-04f,  3.702141132e-04f,  2.779215422e-04f,  1.433343710e-04f,
    -1.739747542e-05f, -1.825268212e-04f, -3.275569582e-04f, -4.292153136e-04f,
    -4.692363252e-04f, -4.380412787e-04f, -3.369927722e-04f, -1.793430647e-04f,
     1.090715583e-05f,  2.018884720e-04f,  3.583687062e-04f,  4.472327656e-04f,
     4.430176144e-04f,  3.328612348e-04f,  1.195705177e-04f, -1.772948447e-04f,
    -5.233395181e-04f, -8.734615283e-04f, -1.178492158e-03f, -1.392626816e-03f,
    -1.481376754e-03f, -1.427912505e-03f, -1.237440796e-03f, -9.380095709e-04f,
    -5.779559152e-04f, -2.196029408e-04f,  6.963652376e-05f,  2.278167679e-04f,
     2.089098269e-04f, -8.522664153e-06f, -4.158730912e-04f, -9.736394247e-04f,
    -1.615274354e-03f, -2.256000275e-03f, -2.804881637e-03f, -3.179185574e-03f,
    -3.318072679e-03f, -3.194179765e-03f, -2.820249823e-03f, -2.250153908e-03f,
    -1.573009078e-03f, -9.015197063e-04f, -3.554025410e-04f, -4.274361260e-05f,
    -4.173957716e-05f, -3.860690621e-04f, -1.056361844e-03f, -1.979559040e-03f,
    -3.036948375e-03f, -4.079888381e-03f, -4.951632742e-03f, -5.511622297e-03f,
    -5.659131303e-03f, -5.351749107e-03f, -4.616115610e-03f, -3.547931021e-03f,
    -2.301304034e-03f, -1.067711182e-03f, -4.797236887e-05f,  5.793528073e-04f,
     6.877684140e-04f,  2.264140354e-04f, -7.670900754e-04f, -2.167024398e-03f,
    -3.774033452e-03f, -5.342256980e-03f, -6.615417362e-03f, -7.366929944e-03f,
    -7.437491488e-03f, -6.764495696e-03f, -5.397524080e-03f, -3.497261509e-03f,
    -1.316338116e-03f,  8.352029296e-04f,  2.635916442e-03f,  3.803131928e-03f,
     4.140168133e-03f,  3.572361410e-03f,  2.166022315e-03f,  1.259434668e-04f,
    -2.229066889e-03f, -4.510164663e-03f, -6.316010900e-03f, -7.295456304e-03f,
    -7.205879021e-03f, -5.957412786e-03f, -3.635662377e-03f, -4.976602176e-04f,
     3.059215916e-03f,  6.552491380e-03f,  9.485963752e-03f,  1.142738320e-02f,
     1.207988808e-02f,  1.133615224e-02f,  9.305222857e-03f,  6.306853196e-03f,
     2.832099798e-03f, -5.250355198e-04f, -3.155870804e-03f, -4.532188204e-03f,
    -4.293635364e-03f, -2.313266687e-03f,  1.270403201e-03f,  6.061342269e-03f,
     1.145507527e-02f,  1.672284137e-02f,  2.111982664e-02f,  2.400113941e-02f,
     2.492820101e-02f,  2.374778003e-02f,  2.063039761e-02f,  1.605953166e-02f,
     1.077131401e-02f,  5.651339065e-03f,  1.603001874e-03f, -5.936230005e-04f,
    -4.104837258e-04f,  2.332885363e-03f,  7.423583718e-03f,  1.426905171e-02f,
     2.197087948e-02f,  2.945032351e-02f,  3.560807945e-02f,  3.949531297e-02f,
     4.046929143e-02f,  3.830892919e-02f,  3.326974916e-02f,  2.606715110e-02f,
     1.778668202e-02f,  9.732191822e-03f,  3.232222227e-03f, -5.668147190e-04f,
    -8.892634827e-04f,  2.529901494e-03f,  9.378042649e-03f,  1.878475405e-02f,
     2.942983967e-02f,  3.972677025e-02f,  4.805686366e-02f,  5.301997517e-02f,
     5.366371394e-02f,  4.965412953e-02f,  4.135872757e-02f,  2.982419536e-02f,
     1.664741499e-02f,  3.754119055e-03f, -6.884709041e-03f, -1.355810767e-02f,
    -1.509551496e-02f, -1.107773553e-02f, -1.942050981e-03f,  1.104010435e-02f,
     2.591611867e-02f,  4.031890829e-02f,  5.181204423e-02f,  5.826416876e-02f,
     5.819633020e-02f,  5.104736579e-02f,  3.731181665e-02f,  1.852244432e-02f,
    -2.928727964e-03f, -2.411201003e-02f, -4.198314016e-02f, -5.384848493e-02f,
    -5.780375501e-02f, -5.307742644e-02f, -4.022109176e-02f, -2.110856390e-02f,
     1.267380355e-03f,  2.318203399e-02f,  4.070206430e-02f,  5.027438626e-02f,
     4.929587576e-02f,  3.657690192e-02f,  1.262233185e-02f, -2.032413486e-02f,
    -5.849457131e-02f, -9.705622726e-02f, -1.307645720e-01f, -1.547277197e-01f,
    -1.651783530e-01f, -1.601419176e-01f, -1.398973711e-01f, -1.071499498e-01f,
    -6.687071486e-02f, -2.580161062e-02f,  8.329089124e-03f,  2.779346091e-02f,
     2.593264225e-02f, -1.830480176e-03f, -5.717278820e-02f, -1.383663663e-01f,
    -2.401724219e-01f, -3.541702791e-01f, -4.695035169e-01f, -5.739708874e-01f,
    -6.553429365e-01f, -7.027521501e-01f, -7.079901395e-01f, -6.665517180e-01f,
    -5.782927415e-01f, -4.476128822e-01f, -2.831308884e-01f, -9.688114938e-02f,
     9.688114938e-02f,  2.831308884e-01f,  4.476128822e-01f,  5.782927415e-01f,
     6.665517180e-01f,  7.079901395e-01f,  7.027521501e-01f,  6.553429365e-01f,
     5.739708874e-01f,  4.695035169e-01f,  3.541702791e-01f,  2.401724219e-01f,
     1.383663663e-01f,  5.717278820e-02f,  1.830480176e-03f, -2.593264225e-02f,
    -2.779346091e-02f, -8.329089124e-03f,  2.580161062e-02f,  6.687071486e-02f,
     1.071499498e-01f,  1.398973711e-01f,  1.601419176e-01f,  1.651783530e-01f,
     1.547277197e-01f,  1.307645720e-01f,  9.705622726e-02f,  5.849457131e-02f,
     2.032413486e-02f, -1.262233185e-02f, -3.657690192e-02f, -4.929587576e-02f,
    -5.027438626e-02f, -4.070206430e-02f, -2.318203399e-02f, -1.267380355e-03f,
     2.110856390e-02f,  4.022109176e-02f,  5.307742644e-02f,  5.780375501e-02f,
     5.384848493e-02f,  4.198314016e-02f,  2.411201003e-02f,  2.928727964e-03f,
    -1.852244432e-02f, -3.731181665e-02f, -5.104736579e-02f, -5.819633020e-02f,
    -5.826416876e-02f, -5.181204423e-02f, -4.031890829e-02f, -2.591611867e-02f,
    -1.104010435e-02f,  1.942050981e-03f,  1.107773553e-02f,  1.509551496e-02f,
     1.355810767e-02f,  6.884709041e-03f, -3.754119055e-03f, -1.664741499e-02f,
    -2.982419536e-02f, -4.135872757e-02f, -4.965412953e-02f, -5.366371394e-02f,
    -5.301997517e-02f, -4.805686366e-02f, -3.972677025e-02f, -2.942983967e-02f,
    -1.878475405e-02f, -9.378042649e-03f, -2.529901494e-03f,  8.892634827e-04f,
     5.668147190e-04f, -3.232222227e-03f, -9.732191822e-03f, -1.778668202e-02f,
    -2.606715110e-02f, -3.326974916e-02f, -3.830892919e-02f, -4.046929143e-02f,
    -3.949531297e-02f, -3.560807945e-02f, -2.945032351e-02f, -2.197087948e-02f,
    -1.426905171e-02f, -7.423583718e-03f, -2.332885363e-03f,  4.104837258e-04f,
     5.936230005e-04f, -1.603001874e-03f, -5.651339065e-03f, -1.077131401e-02f,
    -1.605953166e-02f, -2.063039761e-02f, -2.374778003e-02f, -2.492820101e-02f,
    -2.400113941e-02f, -2.111982664e-02f, -1.672284137e-02f, -1.145507527e-02f,
    -6.061342269e-03f, -1.270403201e-03f,  2.313266687e-03f,  4.293635364e-03f,
     4.532188204e-03f,  3.155870804e-03f,  5.250355198e-04f, -2.832099798e-03f,
    -6.306853196e-03f, -9.305222857e-03f, -1.133615224e-02f, -1.207988808e-02f,
    -1.142738320e-02f, -9.485963752e-03f, -6.552491380e-03f, -3.059215916e-03f,
     4.976602176e-04f,  3.635662377e-03f,  5.957412786e-03f,  7.205879021e-03f,
     7.295456304e-03f,  6.316010900e-03f,  4.510164663e-03f,  2.229066889e-03f,
    -1.259434668e-04f, -2.166022315e-03f, -3.572361410e-03f, -4.140168133e-03f,
    -3.803131928e-03f, -2.635916442e-03f, -8.352029296e-04f,  1.316338116e-03f,
     3.497261509e-03f,  5.397524080e-03f,  6.764495696e-03f,  7.437491488e-03f,
     7.366929944e-03f,  6.615417362e-03f,  5.342256980e-03f,  3.774033452e-03f,
     2.167024398e-03f,  7.670900754e-04f, -2.264140354e-04f, -6.877684140e-04f,
    -5.793528073e-04f,  4.797236887e-05f,  1.067711182e-03f,  2.301304034e-03f,
     3.547931021e-03f,  4.616115610e-03f,  5.351749107e-03f,  5.659131303e-03f,
     5.511622297e-03f,  4.951632742e-03f,  4.079888381e-03f,  3.036948375e-03f,
     1.979559040e-03f,  1.056361844e-03f,  3.860690621e-04f,  4.173957716e-05f,
     4.274361260e-05f,  3.554025410e-04f,  9.015197063e-04f,  1.573009078e-03f,
     2.250153908e-03f,  2.820249823e-03f,  3.194179765e-03f,  3.318072679e-03f,
     3.179185574e-03f,  2.804881637e-03f,  2.256000275e-03f,  1.615274354e-03f,
     9.736394247e-04f,  4.158730912e-04f,  8.522664153e-06f, -2.089098269e-04f,
    -2.278167679e-04f, -6.963652376e-05f,  2.196029408e-04f,  5.779559152e-04f,
     9.380095709e-04f,  1.237440796e-03f,  1.427912505e-03f,  1.481376754e-03f,
     1.392626816e-03f,  1.178492158e-03f,  8.734615283e-04f,  5.233395181e-04f,
     1.772948447e-04f, -1.195705177e-04f, -3.328612348e-04f, -4.430176144e-04f,
    -4.472327656e-04f, -3.583687062e-04f, -2.018884720e-04f, -1.090715583e-05f,
     1.793430647e-04f,  3.369927722e-04f,  4.380412787e-04f,  4.692363252e-04f,
     4.292153136e-04f,  3.275569582e-04f,  1.825268212e-04f,  1.739747542e-05f,
    -1.433343710e-04f, -2.779215422e-04f, -3.702141132e-04f, -4.117814918e-04f,
    -4.022822776e-04f, -3.489785794e-04f, -2.647437785e-04f, -1.657789788e-04f,
    -6.873235982e-05f,  1.175544358e-05f,  6.509917991e-05f,  8.594826909e-05f,
     7.443673327e-05f,  3.566580243e-05f, -2.166811241e-05f, -8.685576674e-05f,
    -1.491204152e-04f, -1.991525379e-04f, -2.304972247e-04f, -2.400598765e-04f,
    -2.284384894e-04f, -1.992102755e-04f, -1.582638875e-04f, -1.124802908e-04f,
    -6.874423496e-05f, -3.274141228e-05f, -8.373029650e-06f,  2.728250001e-06f,
     1.114439836e-06f, -1.080931712e-05f, -2.937702979e-05f, -5.039378470e-05f,
    -6.981265992e-05f, -8.441699810e-05f, -9.210482868e-05f, -9.215232957e-05f,
    -8.510352198e-05f, -7.252468313e-05f, -5.665993882e-05f, -3.991071071e-05f,
    -2.443191610e-05f, -1.209199242e-05f, -3.699234065e-06f,  4.121777932e-07f,
     7.181999416e-07f, -1.832000651e-06f, -6.015902141e-06f, -1.056464084e-05f,
    -1.442620473e-05f, -1.678942486e-05f, -1.736191028e-05f, -4.846057553e-05f,
};

/**
 * FIR 多相解析过采样（L=8）。
 *
 * 链路：x → 多相 FIR 解析上采样（复）→ 复整形 → Re → 实数 remez 抗混叠 → 抽取 ÷8。
 * 一次实低通 remez + 复频移给出单边（解析）滤波器：正频 [0, f_p] 通、镜像与负频落阻带
 * （h_a 的频移量 f_p/2，见 polyphase_analytic_ovs.py 与笔记 §10）。
 *
 * 整形两条路径（见 Cubic）：
 *   * 「多项式」shape：复整形 H(z)=z+z²+z³，输出 Re(z + d·z² + d²·z³)（= Re(H(d·z))/d）；
 *   * 其它 shape：实数 Pwl 作用在 Re(z) 上（解析上采样仍保证 Re(z) 是带限插值）。
 */
class AnalyticOvsFir {
public:
    void Reset() noexcept {
        x_ring_.fill(0.0f);
        w_ring_.fill(0.0f);
        x_head_ = 0;
        w_head_ = 0;
    }

    /// 选用「多项式」shape（复整形），d 为驱动量。
    void SetPolyDrive(float d) noexcept {
        d_ = d;
        d2_ = d * d;
        poly_ = true;
        shape_ = nullptr;
    }

    /// 选用分段线性 shape（作用在 Re(z) 上）。
    void SetPwlShape(const Pwl* f) noexcept {
        shape_ = f;
        poly_ = false;
    }

    /// 推入一个输入样本，返回一个输出样本。
    [[nodiscard]] float Process(float x) noexcept {
        x_head_ = (x_head_ + 1) & (kFirOvsPhaseTaps - 1);
        x_ring_[static_cast<size_t>(x_head_)] = x;
        for (int p = 0; p < kFirOvsL; ++p) {
            float re = 0.0f;
            float im = 0.0f;
            for (int k = 0; k < kFirOvsPhaseTaps; ++k) {
                const int idx = p + kFirOvsL * k;
                const float xv = x_ring_[static_cast<size_t>(
                    (x_head_ - k) & (kFirOvsPhaseTaps - 1))];
                re += kFirOvsAnalyticRe[idx] * xv;
                im += kFirOvsAnalyticIm[idx] * xv;
            }
            float w;
            if (poly_) {
                const std::complex<float> z{re, im};
                const std::complex<float> z2 = z * z;
                const std::complex<float> z3 = z2 * z;
                // Re((d·z + (d·z)² + (d·z)³))/d = Re(z + d·z² + d²·z³)
                w = re + d_ * z2.real() + d2_ * z3.real();
            }
            else {
                w = shape_->Value(re);
            }
            w_head_ = (w_head_ + 1) & (kFirOvsTaps - 1);
            w_ring_[static_cast<size_t>(w_head_)] = w;
        }
        // 抽取：y[m] = Σ_j h_d[j]·w[mL−j]；刚写入的最新 w 在位置 mL+L−1。
        float acc = 0.0f;
        for (int j = 0; j < kFirOvsTaps; ++j) {
            acc += kFirOvsDecim[j] * w_ring_[static_cast<size_t>(
                (w_head_ - (kFirOvsL - 1) - j) & (kFirOvsTaps - 1))];
        }
        return acc;
    }

private:
    std::array<float, kFirOvsPhaseTaps> x_ring_{};
    std::array<float, kFirOvsTaps> w_ring_{};
    int x_head_{};
    int w_head_{};
    float d_{1.0f};
    float d2_{1.0f};
    bool poly_{true};
    const Pwl* shape_{};
};

// ------------------------------------------------------------
// (B) IIR 全通和多相解析过采样（L=8，并行全通和）
// ------------------------------------------------------------
// 与 labs/adaa_iir/polyphase_allpass_iir.py 逐式一致：
//   半带：椭圆半带 N=19、rs=100 dB、Wp=0.4726712（通带边沿 = 该级 fs/4）；
//   节参数 a = ρ²（按 |p| 升序逐节交替分到两条链：链0 5 节、链1 4 节）；
//   低速率全通链 A(ζ) = Π (ζ⁻¹ + a_k)/(1 + a_k ζ⁻¹)，ζ⁻¹ = z⁻²（分支跑该级半速率）；
//   半带插值 H_hb = ½[A0(ζ) + z⁻¹·A1(ζ)]；
//   8× = **三级级联**半带插值（48k→96k→192k→384k，三级同一组系数）：
//     限制通带的只有第一级，后两级交叉点 48/96 kHz 远高于信号带，只压各自镜像；
//   解析（单边）链 = 半带插值 ∘ 解析全通 H_an = ½[A0(−ζ) + j·z⁻¹·A1(−ζ)]，
//     位置与 2× 版一致（最后一级之后、抽取之前），抽取端同样三级级联。

static constexpr std::array<float, 5> kOvsIirChain0 = {
    3.993399361e-02f, 2.953322197e-01f, 5.952248425e-01f, 8.146834120e-01f, 9.656342782e-01f,
};
static constexpr std::array<float, 4> kOvsIirChain1 = {
    1.478918562e-01f, 4.516440986e-01f, 7.163745215e-01f, 8.951804120e-01f,
};

/**
 * 低速率一阶全通链（ζ⁻¹ = z⁻²）。
 *
 * `rotated` 表示 A(−ζ)（解析/单边路径）：节 H(ζ⁻¹) = (ζ⁻¹ + a)/(1 + a·ζ⁻¹)，
 * ζ → −ζ 后奇次项取负，等价递推 out = a·in − s（未旋转为 out = a·in + s）。
 */
template <int N>
struct AllpassChain {
    void Init(const std::array<float, N>& a, bool rotated) noexcept {
        a_ = a;
        rotated_ = rotated;
        s_.fill(0.0f);
    }

    void Reset() noexcept { s_.fill(0.0f); }

    [[nodiscard]] float Tick(float x) noexcept {
        float y = x;
        for (int i = 0; i < N; ++i) {
            const float a = a_[static_cast<size_t>(i)];
            const float s = s_[static_cast<size_t>(i)];
            float in;
            if (rotated_) {
                in = y + a * s;
                y = a * in - s;
            }
            else {
                in = y - a * s;
                y = a * in + s;
            }
            s_[static_cast<size_t>(i)] = in;
        }
        return y;
    }

private:
    std::array<float, N> a_{};
    std::array<float, N> s_{};
    bool rotated_{};
};

/**
 * IIR 全通和多相解析过采样（L=8）。
 *
 * 链路：x →（半带插值三级级联 48k→96k→192k→384k ∘ 解析全通）解析上采样 → 复整形 → Re
 *       → 半带抗混叠三级级联 ÷8 → y @48k。
 * 各级的两条全通链都跑在该级输入速率（级1 48k、级2 96k、级3 192k、解析级 192k）；
 * 整形路径同 AnalyticOvsFir。每输入采样滤波器 MAC：198（实链路 126）。
 */
class AnalyticOvsIir {
public:
    AnalyticOvsIir() noexcept {
        up0_a0_.Init(kOvsIirChain0, false);
        up0_a1_.Init(kOvsIirChain1, false);
        up1_a0_.Init(kOvsIirChain0, false);
        up1_a1_.Init(kOvsIirChain1, false);
        up2_a0_.Init(kOvsIirChain0, false);
        up2_a1_.Init(kOvsIirChain1, false);
        an_ve_re_.Init(kOvsIirChain0, true);
        an_vo_re_.Init(kOvsIirChain0, true);
        an_ve_im_.Init(kOvsIirChain1, true);
        an_vo_im_.Init(kOvsIirChain1, true);
        dec0_a0_.Init(kOvsIirChain0, false);
        dec0_a1_.Init(kOvsIirChain1, false);
        dec1_a0_.Init(kOvsIirChain0, false);
        dec1_a1_.Init(kOvsIirChain1, false);
        dec2_a0_.Init(kOvsIirChain0, false);
        dec2_a1_.Init(kOvsIirChain1, false);
    }

    void Reset() noexcept {
        up0_a0_.Reset();
        up0_a1_.Reset();
        up1_a0_.Reset();
        up1_a1_.Reset();
        up2_a0_.Reset();
        up2_a1_.Reset();
        an_ve_re_.Reset();
        an_vo_re_.Reset();
        an_ve_im_.Reset();
        an_vo_im_.Reset();
        dec0_a0_.Reset();
        dec0_a1_.Reset();
        dec1_a0_.Reset();
        dec1_a1_.Reset();
        dec2_a0_.Reset();
        dec2_a1_.Reset();
        a1vo_prev_ = 0.0f;
        co0_prev_ = 0.0f;
        co1_prev_ = 0.0f;
        co2_prev_ = 0.0f;
    }

    /// 选用「多项式」shape（复整形），d 为驱动量。
    void SetPolyDrive(float d) noexcept {
        d_ = d;
        d2_ = d * d;
        poly_ = true;
        shape_ = nullptr;
    }

    /// 选用分段线性 shape（作用在 Re(z) 上）。
    void SetPwlShape(const Pwl* f) noexcept {
        shape_ = f;
        poly_ = false;
    }

    [[nodiscard]] float Process(float x) noexcept {
        // 半带插值三级级联：每级输入一个样本、输出两个（偶 = A0、奇 = A1）。
        const float u1[2] = {up0_a0_.Tick(x), up0_a1_.Tick(x)};
        float u2[4];
        for (int i = 0; i < 2; ++i) {
            u2[2 * i] = up1_a0_.Tick(u1[i]);
            u2[2 * i + 1] = up1_a1_.Tick(u1[i]);
        }
        float u3[8];
        for (int i = 0; i < 4; ++i) {
            u3[2 * i] = up2_a0_.Tick(u2[i]);
            u3[2 * i + 1] = up2_a1_.Tick(u2[i]);
        }
        // 解析全通：u3 的偶/奇各半；y[2m]=A0r(ve)+j·A1r(vo)[m−1]、y[2m+1]=A0r(vo)+j·A1r(ve)。
        float w[8];
        for (int m = 0; m < 4; ++m) {
            const float ve = u3[2 * m];
            const float vo = u3[2 * m + 1];
            const float re_e = an_ve_re_.Tick(ve);
            const float re_o = an_vo_re_.Tick(vo);
            const float im_e = an_ve_im_.Tick(ve);
            const float im_o = an_vo_im_.Tick(vo);
            w[2 * m] = Shaped(re_e, a1vo_prev_);
            w[2 * m + 1] = Shaped(re_o, im_e);
            a1vo_prev_ = im_o;
        }
        // 半带抽取三级级联：z[m] = ½[A0(w 偶)[m] + A1(w 奇)[m−1]]。
        float d1[4];
        for (int m = 0; m < 4; ++m) {
            const float ce = dec0_a0_.Tick(w[2 * m]);
            const float co = dec0_a1_.Tick(w[2 * m + 1]);
            d1[m] = 0.5f * (ce + co0_prev_);
            co0_prev_ = co;
        }
        float d2[2];
        for (int m = 0; m < 2; ++m) {
            const float ce = dec1_a0_.Tick(d1[2 * m]);
            const float co = dec1_a1_.Tick(d1[2 * m + 1]);
            d2[m] = 0.5f * (ce + co1_prev_);
            co1_prev_ = co;
        }
        const float ce = dec2_a0_.Tick(d2[0]);
        const float co = dec2_a1_.Tick(d2[1]);
        const float out = 0.5f * (ce + co2_prev_);
        co2_prev_ = co;
        return out;
    }

private:
    [[nodiscard]] float Shaped(float re, float im) const noexcept {
        if (!poly_) { return shape_->Value(re); }
        const std::complex<float> z{re, im};
        const std::complex<float> z2 = z * z;
        const std::complex<float> z3 = z2 * z;
        // Re((d·z + (d·z)² + (d·z)³))/d = Re(z + d·z² + d²·z³)
        return re + d_ * z2.real() + d2_ * z3.real();
    }

    AllpassChain<5> up0_a0_;
    AllpassChain<4> up0_a1_;
    AllpassChain<5> up1_a0_;
    AllpassChain<4> up1_a1_;
    AllpassChain<5> up2_a0_;
    AllpassChain<4> up2_a1_;
    AllpassChain<5> an_ve_re_;
    AllpassChain<5> an_vo_re_;
    AllpassChain<4> an_ve_im_;
    AllpassChain<4> an_vo_im_;
    AllpassChain<5> dec0_a0_;
    AllpassChain<4> dec0_a1_;
    AllpassChain<5> dec1_a0_;
    AllpassChain<4> dec1_a1_;
    AllpassChain<5> dec2_a0_;
    AllpassChain<4> dec2_a1_;
    float d_{1.0f};
    float d2_{1.0f};
    bool poly_{true};
    const Pwl* shape_{};
    float a1vo_prev_{};
    float co0_prev_{};
    float co1_prev_{};
    float co2_prev_{};
};

} // namespace adaa
