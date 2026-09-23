#pragma once
#include "iir_design.hpp"
#include <algorithm>
#include <array>

namespace qwqdsp_filter {
/**
 * @brief shelving 滤波器原型设计(highshelf)
 *
 * 每个函数输出"中点在 (1)rad/sec 的 highshelf 原型", 之后照常经过
 * IIRDesign::ProtyleToLowpass 做频率缩放、IIRDesign::Bilinear 离散化。
 *
 * 幅度语义(与常见 EQ 的 highshelf 一致):
 * - 低频端 |H| = 0dB
 * - 中频(1)rad/sec 处 |H| = boost/2 dB
 * - 高频端 |H| = boost dB
 *
 * 构造用书里 (10.12)/(10.13a) 的混合式 |H|^2 = (1 + b^2*F^2)/(1 + F^2/b^2),
 * 其中 F 是原型的归一化特征函数(|G|^2 = 1/(1+F^2), F(1) = 1):
 * H 的极点取"特征 F/b 的原型"的极点, 零点取"特征 b*F 的原型"的极点,
 * 于是零点与极点互为倒易、H 最小相位, 低频/高频增益由特征函数的首末值决定。
 * Butterworth 下这个式子与低通比值法 (10.6) 等价(书 10.2 有说明); 椭圆则直接是
 * (10.20)/(10.21): 分子分母两个椭圆低通的传输零点完全相同, 对消后不出现。
 *
 * @note 输出节数仍是 num_filter, 每节都带有限零点(不再有"零点在无穷远"的节)
 * @note 与 IIRDesign 的原型一样, 这里所有频率都以 (1)rad/sec 为参考
 * @ref Vadim Zavalishin, "The Art of VA Filter Design", 2.1.2, Chapter 10
 */
struct IIRDesignShelf {
    using ZPK = IIRDesign::ZPK;
    static constexpr auto pi = IIRDesign::pi;

    /// 单次设计支持的极点对数上限(临时缓冲大小, 只影响断言)
    static constexpr size_t kMaxPrototypePairs = 64;

    /// |boost| 小于这个值时按直通处理(boost = 0 时下面几个公式是 0/0)
    static constexpr double kMinBoostDb = 1.0e-9;
private:
    /**
     * @brief 巴特沃斯单位圆极点(偶数阶)
     * @param ret 输出, 至少 num_filter 项, 只写 p
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     */
    static void UnitButterworthPoles(std::span<ZPK> ret, size_t num_filter) noexcept {
        size_t const n = 2 * num_filter;
        for (size_t k = 1; k <= num_filter; ++k) {
            double const phi = (2.0 * static_cast<double>(k) - 1.0) * pi / (2.0 * static_cast<double>(n));
            ret[k - 1].p = std::complex{-std::sin(phi), std::cos(phi)};
        }
    }

    /**
     * @brief 按 |H(0)| = 1 定整体增益(写进 section 0), 并返回它
     *
     * 逐节 k = 1 时 H(0) = prod|z|^2/prod|p|^2, s->inf 时 H(inf) = prod k,
     * 所以把 gain = prod|p|^2/prod|z|^2 放在 section 0 上, 既标定了直流增益,
     * 又让 H(inf) 正好等于同一个 gain。
     *
     * @param ret 设计好的节(每节都必须带有限零点)
     * @param num_filter 节数
     * @return 整体增益, 也就是 |H(inf)|
     */
    static double NormalizeDcGain(std::span<ZPK> ret, size_t num_filter) noexcept {
        double num = 1.0;
        double den = 1.0;
        for (size_t i = 0; i < num_filter; ++i) {
            assert(ret[i].z);
            num *= std::norm(ret[i].p);
            den *= std::norm(*ret[i].z);
        }
        double const gain = num / den;
        for (size_t i = 0; i < num_filter; ++i) {
            ret[i].k = 1.0;
        }
        ret[0].k = gain;
        return gain;
    }

    /**
     * @brief 把一整串节都设成"直通"(0dB)的退化设计
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 节数
     * @note 零极点重合在 s = -1, 每节 H(z) 恒等于 1; 用于 boost = 0 的退化输入
     */
    static void MakeBypass(std::span<ZPK> ret, size_t num_filter) noexcept {
        for (size_t i = 0; i < num_filter; ++i) {
            ret[i].z = std::complex{-1.0, 0.0};
            ret[i].p = std::complex{-1.0, 0.0};
            ret[i].k = 1.0;
        }
    }

    /**
     * @brief 把 eps 换算成原型函数要的"边沿处衰减(dB)"
     * @param eps 特征函数的幅度参数
     * @return 10*log10(1 + eps^2)
     */
    static double DbOfEps(double eps) noexcept {
        return 10.0 * std::log10(1.0 + eps * eps);
    }

    /**
     * @brief DbOfEps 的逆
     * @param db 边沿处衰减(dB)
     * @return sqrt(10^(db/10) - 1)
     */
    static double EpsOfDb(double db) noexcept {
        return std::sqrt(std::pow(10.0, db / 10.0) - 1.0);
    }

    /**
     * @brief 复算椭圆原型的设计模数 k(阻带边沿在 1/k)
     *
     * IIRDesign::Elliptic 内部用同一套 cephes 工具由 nome 反解 k 但不对外暴露,
     * 这里复算一遍, 好让 H 的零点组与极点组共享同一个频率尺度。
     *
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param ripple_db 通带纹波(dB)
     * @param stopband_db 阻带衰减(dB)
     * @return 设计模数 k; 规格过陡时可能返回 >= 1 的无效值, 由调用者判
     */
    static double DesignModulus(size_t num_filter, double ripple_db, double stopband_db) noexcept {
        double const eps_p = EpsOfDb(ripple_db);
        double const eps_s = EpsOfDb(stopband_db);
        double const m1 = eps_p / eps_s;
        double const m1_sq = m1 * m1;
        size_t const n = 2 * num_filter;
        double const kk1 = IIRDesign::EllipticHelperCephes::Ellpk(1.0 - m1_sq);
        double const kpk1 = IIRDesign::EllipticHelperCephes::Ellpk(m1_sq);
        double const q = std::exp(-pi * kpk1 / (static_cast<double>(n) * kk1));
        return IIRDesign::EllipticHelperCephes::Cay(q);
    }

public:
    /**
     * @brief 巴特沃斯 highshelf 原型
     *
     * 特征函数 F = w^N(N = 2*num_filter); "特征缩放 b 倍"等价于极点半径乘
     * b^(-1/N)。巴特沃斯有 F(0) = 0, 直流增益本来就落在 1, 因此整体增益直接取
     * b^2, 正好等于幅度棚高。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param boost_db 高频端相对低频端的提升(dB), 可正可负(负值即反向 shelf),
     *        为 0 时输出直通
     * @note 每节的零点与极点互为倒易(比值 b^2), 频响在对数频率轴上严格对称
     */
    static void Butterworth(std::span<ZPK> ret, size_t num_filter, double boost_db) {
        assert(ret.size() >= num_filter);
        assert(num_filter <= kMaxPrototypePairs);

        if (std::abs(boost_db) < kMinBoostDb) {
            MakeBypass(ret, num_filter);
            return;
        }
        size_t const n = 2 * num_filter;
        double const amp = std::pow(10.0, boost_db / 20.0);  // 幅度棚高 |H(inf)|
        double const b = std::sqrt(amp);
        double const s = std::pow(b, 1.0 / static_cast<double>(n));

        UnitButterworthPoles(ret, num_filter);
        for (size_t i = 0; i < num_filter; ++i) {
            // 零点取"特征 b*F"的原型极点(半径 b^(-1/N)), 极点取"特征 F/b"的(半径 b^(1/N))
            ret[i].z = ret[i].p / s;
            ret[i].p *= s;
            ret[i].k = 1.0;
        }
        ret[0].k = amp;
    }

    /**
     * @brief 切比雪夫 I 型 highshelf 原型
     *
     * 特征函数 F = eps*T_N(w/lam)。偶数阶时 F(0) = eps != 0, 低频端会跟着抬起来,
     * 所以这里把 ripple_db 定义成"成品 shelf 通带的实际纹波", 由它反解 eps 与 lam,
     * 使 0dB / boost/2 dB / boost dB 三个电平同时精确成立。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param boost_db 高频端相对低频端的提升(dB), 可正可负(负值即反向 shelf),
     *        为 0 时输出直通
     * @param ripple_db 通带纹波(dB, >0); 提升时纹波落在 0dB 参考线下方, 衰减时落上方
     * @note ripple_db 趋近 0 时退化成巴特沃斯形状
     */
    static void Chebyshev1(std::span<ZPK> ret, size_t num_filter, double boost_db, double ripple_db) {
        assert(ret.size() >= num_filter);
        assert(num_filter <= kMaxPrototypePairs);

        if (std::abs(boost_db) < kMinBoostDb) {
            MakeBypass(ret, num_filter);
            return;
        }
        size_t const n = 2 * num_filter;
        double const g = std::pow(10.0, boost_db / 10.0);   // 棚高(功率比)
        double const d = std::pow(10.0, ripple_db / 10.0);  // 纹波深度(功率比)
        // 通带参考点(直流, F = eps 处)的原始电平: 提升时直流取纹波上沿, 衰减时取下沿
        double const a0 = (boost_db > 0.0) ? d : 1.0 / d;
        double const b2 = std::sqrt(g * a0);                // b^2, 且 b^4 = g*a0 是高频端原始电平
        double const b = std::sqrt(b2);
        double const eps2 = (1.0 - a0) * b2 / (a0 - b2 * b2);
        double const eps = std::sqrt(eps2 > 0.0 ? eps2 : 0.0);

        // 中点的原始电平: |H(1)|^2/a0 应为 sqrt(g)(即 boost/2 dB), 由此解特征值 f
        double const a = std::sqrt(g) * a0;
        double const f2 = (a - 1.0) / (b2 - a / b2);
        double const f = (f2 > 0.0) ? std::sqrt(f2) : 0.0;
        // 频率标定: eps*T_N(1/lam) = f
        double lam = 1.0;
        if (eps > 0.0 && f > eps) {
            lam = 1.0 / std::cosh(std::acosh(f / eps) / static_cast<double>(n));
        }

        std::array<ZPK, kMaxPrototypePairs> scratch{};
        // 零点: "特征 b*F"的原型极点
        IIRDesign::Chebyshev1(std::span{scratch}.first(num_filter), num_filter, DbOfEps(eps * b), false);
        for (size_t i = 0; i < num_filter; ++i) {
            ret[i].z = scratch[i].p * lam;
        }
        // 极点: "特征 F/b"的原型极点
        IIRDesign::Chebyshev1(std::span{scratch}.first(num_filter), num_filter, DbOfEps(eps / b), false);
        for (size_t i = 0; i < num_filter; ++i) {
            ret[i].p = scratch[i].p * lam;
            ret[i].k = 1.0;
        }
        NormalizeDcGain(ret, num_filter);
    }

    /**
     * @brief 切比雪夫 II 型 highshelf 原型
     *
     * 特征函数 F(w) = 1/(eps*T~_N(m/w)), 于是 |H|^2 = 1/(1+F^2) 与
     * IIRDesign::Chebyshev2Natural 一致。要拼成 shelf, H1/H2 只差一个 eps:
     * 两个原型的传输零点必须逐点相同才能相消, 这正是自然归一化提供的东西
     * (IIRDesign::Chebyshev2 的零点会随 ripple 移动, 拿来做 shelf 会在零点处炸出尖峰)。
     *
     * 记 q = T~(m/w)^2, 则 |H|^2 = (eps1^2/eps2^2)*(1+eps2^2 q)/(1+eps1^2 q), 于是
     * - 棚高 (q -> 0):  (eps1/eps2)^2 = beta^4
     * - 棚带起伏 (q 在 0..1 之间摆动): (1+eps2^2)/(1+beta^4 eps2^2)
     * - 中点 (w = 1):  T~(m) = 1/(eps2 beta)
     * (eps1, eps2, m) 三个未知量正好由上面三个条件定出, 其中 ripple_db 是输入旋钮。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param boost_db 高频端相对低频端的提升(dB), 可正可负(负值即反向 shelf),
     *        为 0 时输出直通
     * @param ripple_db 棚带等波纹深度(dB), 会被夹到 (0, |boost_db|) 内
     * @return 成功返回 true
     * @retval false 参数退化(boost 为 0 不算, 那种情况输出直通)
     * @note 低频端"最平坦", 等波纹落在**棚带**里(与 Chebyshev1 shelf 相反, 那个的
     *       纹波在低频通带)。ripple_db 越小, 棚带越平, 但过渡带越宽(m 越大)。
     * @note 用偶数阶修正原型, 否则高频端会饱和在半高(|H(inf)| = beta 而不是 beta^2)
     * @note 频率尺度 m 的实现方式是把原型的极点与零点整体乘 m(传输零点一起搬, 所以
     *       两组零点仍然逐点相同; 无穷远处的那个零点保持不动)
     */
    [[nodiscard]] static bool Chebyshev2(std::span<ZPK> ret, size_t num_filter, double boost_db,
                                         double ripple_db) {
        assert(ret.size() >= num_filter);
        assert(num_filter <= kMaxPrototypePairs);

        if (std::abs(boost_db) < kMinBoostDb) {
            MakeBypass(ret, num_filter);
            return true;
        }
        // 衰减方向 = 提升方向的倒数: 先按 |boost| 设计, 最后交换零点与极点
        bool const invert = boost_db < 0.0;
        double const beta2 = std::pow(10.0, std::abs(boost_db) / 20.0);   // 幅度棚高
        double const beta4 = beta2 * beta2;
        // ripple -> 0 时过渡带无限宽, ripple -> |boost| 时退化成半高, 两端都夹一下
        double const ripple = std::clamp(ripple_db, 1.0e-3, std::abs(boost_db) * 0.999);
        double const rho = std::pow(10.0, -ripple / 10.0);
        // (1+eps2^2)/(1+beta^4 eps2^2) = rho
        double const eps2_sq = (1.0 - rho) / (rho * beta4 - 1.0);
        if (!(eps2_sq > 0.0)) {
            return false;
        }
        double const eps2 = std::sqrt(eps2_sq);
        double const eps1 = beta2 * eps2;                                  // (eps1/eps2)^2 = beta^4

        // 频率尺度 m: 解 T~(m) = 1/(eps2*beta), 使中点落在 (1)rad/sec
        size_t const n = 2 * num_filter;
        double const c = std::cos(pi * (static_cast<double>(n) - 1.0) / (2.0 * static_cast<double>(n)));
        double const r = 1.0 / (eps2 * std::sqrt(beta2));
        double u = (r >= 1.0) ? std::cosh(std::acosh(r) / static_cast<double>(n))
                              : std::cos(std::acos(std::min(r, 1.0)) / static_cast<double>(n));
        u = std::max(u, c);
        double const m = std::sqrt((u * u - c * c) / (1.0 - c * c));
        if (!(m > 0.0) || !std::isfinite(m)) {
            return false;
        }
        // eps -> Chebyshev2Natural 要的阻带衰减 dB
        auto const rs_of_eps = [](double eps) {
            return 10.0 * std::log10(eps * eps / (1.0 + eps * eps));
        };

        std::array<ZPK, kMaxPrototypePairs> scratch{};
        // 第一批: eps1 的原型, 极点乘 m 后当 H 的极点(衰减方向则当零点)
        IIRDesign::Chebyshev2Natural(std::span{scratch}.first(num_filter), num_filter, rs_of_eps(eps1), true);
        for (size_t i = 0; i < num_filter; ++i) {
            if (invert) {
                ret[i].z = scratch[i].p * m;
            }
            else {
                ret[i].p = scratch[i].p * m;
            }
        }
        // 第二批: eps2 的原型, 极点乘 m 后当 H 的零点(衰减方向则当极点)
        IIRDesign::Chebyshev2Natural(std::span{scratch}.first(num_filter), num_filter, rs_of_eps(eps2), true);
        for (size_t i = 0; i < num_filter; ++i) {
            if (invert) {
                ret[i].p = scratch[i].p * m;
            }
            else {
                ret[i].z = scratch[i].p * m;
            }
            ret[i].k = 1.0;
        }
        NormalizeDcGain(ret, num_filter);
        return true;
    }

    /**
     * @brief 椭圆(考尔) highshelf 原型
     *
     * 直接对应书里的 (10.20)/(10.21): 分子分母各是一个椭圆低通, 两者的传输零点
     * (虚轴上那对零点)完全相同, 对消后消失, 只剩两组左半平面极点分别当 H 的零点
     * 与极点, 因此不会出现"分子零点处塌成 0、分母零点处冲成无穷"的毛病。
     *
     * 椭圆的设计模数 k(即过渡带宽度)只由纹波与阶数决定, 没有切比雪夫那样的额外
     * 频率标定自由度: 书里 (10.14) 注明偶数阶 EMQF 的参考点不是 w = 0 与 w = inf,
     * 而是特征函数取 0 与 inf 的那两个点, 因此这里只有直流端精确, 高频端的幅度会比
     * |boost| 少一点(纹波"向内", 即书 Fig. 10.13 的形状)。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param boost_db 高频端相对低频端的提升(dB), 可正可负(负值即反向 shelf),
     *        为 0 时输出直通
     * @param ripple_db 通带纹波(dB, >0), 即椭圆原型的通带边沿电平
     * @param stopband_db 阻带衰减(dB, >ripple_db), 决定过渡带陡度
     * @return 成功返回 true
     * @retval false 规格过陡(阻带边沿与通带边沿重合, 模数在 double 下就是 1),
     *         与 IIRDesign::Elliptic 的失败条件一致
     * @note 高频端的缺口(boost 与高频端幅度之差)随 ripple_db 增大、stopband_db 减小
     *       而变大。实测(阶数 2~16, |boost| <= 24dB):
     *       - stopband >= 40dB 时 缺口 <= 0.4dB, 中点偏差 <= 0.25dB;
     *       - stopband = 30dB 时 缺口 <= 1.2dB, 中点偏差 <= 0.6dB;
     *       - 阻带很浅且 boost 很大时(例如 ripple=6dB / stopband=10dB / boost=40dB)
     *         缺口可超过 30dB, 这种参数下应当改用巴特沃斯或切比雪夫
     */
    [[nodiscard]] static bool Elliptic(std::span<ZPK> ret, size_t num_filter, double boost_db,
                                       double ripple_db, double stopband_db) {
        assert(ret.size() >= num_filter);
        assert(num_filter <= kMaxPrototypePairs);

        if (std::abs(boost_db) < kMinBoostDb) {
            MakeBypass(ret, num_filter);
            return true;
        }
        double const eps_p = EpsOfDb(ripple_db);
        double const eps_s = EpsOfDb(stopband_db);
        if (!(eps_s > eps_p)) {
            return false;
        }
        double const m1 = eps_p / eps_s;
        double const root_m1 = std::sqrt(m1);
        // 提升与衰减是同一个方程的两个方向: 用 |boost| 定 b, 再按需要交换零点与极点
        double const g = std::pow(10.0, std::abs(boost_db) / 10.0);
        double const k_ref = DesignModulus(num_filter, ripple_db, stopband_db);
        if (!(k_ref > 0.0 && k_ref < 1.0)) {
            return false;
        }
        // 解 b: b^4 = g * (1 + b^2*m1)/(1 + m1/b^2), 即 b^4 + m1*b^2*(1-g) - g = 0
        double const c = m1 * (1.0 - g);
        double const y = (std::sqrt(c * c + 4.0 * g) - c) / 2.0;  // y = b^2
        if (!(y > 0.0)) {
            return false;
        }
        double const b = std::sqrt(y);
        // 衰减方向: 把"特征 b*F"与"特征 F/b"的角色互换
        double const e_num = (boost_db > 0.0) ? root_m1 * b : root_m1 / b;
        double const e_den = (boost_db > 0.0) ? root_m1 / b : root_m1 * b;

        std::array<ZPK, kMaxPrototypePairs> scratch{};
        double const scale = std::sqrt(k_ref);
        // 零点: "特征 e_num"的椭圆原型极点
        if (!IIRDesign::Elliptic(
                std::span{scratch}.first(num_filter), num_filter, DbOfEps(e_num), DbOfEps(e_num / m1))) {
            return false;
        }
        for (size_t i = 0; i < num_filter; ++i) {
            ret[i].z = scratch[i].p * scale;
        }
        // 极点: "特征 e_den"的椭圆原型极点
        if (!IIRDesign::Elliptic(
                std::span{scratch}.first(num_filter), num_filter, DbOfEps(e_den), DbOfEps(e_den / m1))) {
            return false;
        }
        for (size_t i = 0; i < num_filter; ++i) {
            ret[i].p = scratch[i].p * scale;
            ret[i].k = 1.0;
        }
        NormalizeDcGain(ret, num_filter);
        return true;
    }
};
} // namespace qwqdsp_filter
