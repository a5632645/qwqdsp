#pragma once
#include "iir_design.hpp"
#include <algorithm>
#include <array>
#include <limits>

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
 * - 幅度极值 = boost dB; 未修正的 Chebyshev2 / 椭圆在 w -> inf 处到不了它(极值落在棚带
 *   内的传输零点上), 其余原型在 w -> inf 处就是极值
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

    /**
     * @brief 由原型节复算特征函数 F(w) = sqrt(1/|G(jw)|^2 - 1)
     * @param sections 原型节(库的原型函数输出, 零点为空表示无穷远)
     * @param num_filter 节数
     * @param w 频率(rad/sec), 原型以 (1)rad/sec 为通带边沿
     * @return F(w); |G| = 1 处为 0
     * @note 就是原型定义式 |G|^2 = 1/(1+F^2) 的反解, 用来查形状函数
     */
    static double PrototypeChar(std::span<ZPK const> sections, size_t num_filter, double w) noexcept {
        std::complex<double> const jw{0.0, w};
        std::complex<double> h{1.0, 0.0};
        for (size_t i = 0; i < num_filter; ++i) {
            ZPK const& s = sections[i];
            std::complex<double> num{1.0, 0.0};
            if (s.z) {
                num = (jw - *s.z) * (jw - std::conj(*s.z));
            }
            h *= s.k * num / ((jw - s.p) * (jw - std::conj(s.p)));
        }
        double const v = 1.0 / std::norm(h) - 1.0;
        return v > 0.0 ? std::sqrt(v) : 0.0;
    }

    /**
     * @brief 中点条件在"低频通带电平偏离 0dB 不超过 3dB"时需要的 R(1/lambda)
     *
     * 未修正的椭圆 shelf 里 R(1/lambda) 越大, 低频通带的直流电平 C 越靠近 1; 提升时 R 小到
     * 1 + amp 的平方根、衰减时 R 大到 1 + 1/amp 的平方根, C 就会塌到 0 或冲成无穷。这里
     * 取使 C = 1/2(提升) / 2(衰减) 的值, 也就是低频通带偏离 0dB 不超过 3dB。
     *
     * @param amp 中点电平(幅度), = 10^(boost/40)
     * @return 目标形状函数值
     */
    static double RhoForBandDeviation(double amp) noexcept {
        return std::sqrt((amp + 1.0) * std::max(2.0 - 1.0 / amp, 2.0 / amp - 1.0));
    }

    /**
     * @brief 把频率尺度推进椭圆原型阻带的单调段, 求形状函数达到目标值的频率
     *
     * R(1/lambda) 在 [1/sqrt(k), 1/k] 上从 R(1/sqrt(k)) 单调升到 1/m1, 在 [1/k, 第一个
     * 传输零点] 上再从 1/m1 单调升到无穷, 目标值落在哪一段就二分哪一段。
     *
     * @param ref 参考原型(eps = 1)的节
     * @param num_filter 节数
     * @param m1 eps_p/eps_s, 1/m1 是 R 在阻带边沿的值
     * @param w_start 搜索起点(1/sqrt(k)), 该处形状函数值必须小于目标值
     * @param rho_target 目标形状函数值
     * @return 刚达到目标值的频率; 区间退化时返回起点
     */
    static double ShapeInverse(std::span<ZPK const> ref, size_t num_filter, double m1, double w_start,
                               double rho_target) noexcept {
        double w_lo = w_start;
        double w_hi;
        if (rho_target <= 1.0 / m1) {
            w_hi = w_start * w_start; // R(1/k) = 1/m1
        }
        else {
            // 比阻带边沿还大的值只能落在第一个传输零点之前的那一段
            w_hi = std::numeric_limits<double>::max();
            for (size_t i = 0; i < num_filter; ++i) {
                if (ref[i].z) {
                    w_hi = std::min(w_hi, std::abs(*ref[i].z));
                }
            }
        }
        if (!(w_hi > w_lo)) {
            return w_lo;
        }
        for (int i = 0; i < 80; ++i) {
            double const mid = std::sqrt(w_lo * w_hi);
            if (PrototypeChar(ref, num_filter, mid) < rho_target) {
                w_lo = mid;
            }
            else {
                w_hi = mid;
            }
        }
        return w_hi;
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
        double const amp = std::pow(10.0, boost_db / 20.0); // 幅度棚高 |H(inf)|
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
     * 特征函数 F = eps*Phi(w/lam), 未修正时 Phi = T_N, 修正时 Phi = T_N(sqrt((1-c^2)w^2+c^2))
     * (c 是最小的那个切比雪夫节点, 修正把直流钉到 0dB)。三次电平(直流 / 中点 / 棚顶)
     * 由 (eps, lam) 精确满足, ripple_db 决定 eps 也就是成品 shelf 低频通带的实际纹波:
     * - 未修正: F(0) = eps != 0, 直流是纹波的一侧(提升时通带落在 0dB 下方, 衰减时上方);
     * - 修正后: F(0) = 0(直流定死 0dB), 通带只能向另一侧鼓出(提升时向 0dB 上方,
     *   衰减时下方), 因此要它落在 0dB 下方就用不了修正版。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param boost_db 高频端相对低频端的提升(dB), 可正可负(负值即反向 shelf),
     *        为 0 时输出直通
     * @param ripple_db 通带纹波(dB, >0), 含义见上
     * @param even_modify 用偶数阶修正原型; 修正后 ripple_db 会被夹到 (0, |boost_db|)
     * @note ripple_db 趋近 0 时退化成巴特沃斯形状
     */
    static void Chebyshev1(std::span<ZPK> ret, size_t num_filter, double boost_db, double ripple_db, bool even_modify) {
        assert(ret.size() >= num_filter);
        assert(num_filter <= kMaxPrototypePairs);

        if (std::abs(boost_db) < kMinBoostDb) {
            MakeBypass(ret, num_filter);
            return;
        }
        size_t const n = 2 * num_filter;
        double const g = std::pow(10.0, boost_db / 10.0); // 棚高(功率比)
        // 修正后 eps 的分母是 a0 - b^4 = a0 - g, 要它为正就得 ripple < |boost|
        double const ripple = even_modify ? std::clamp(ripple_db, 1.0e-3, std::abs(boost_db) * 0.999) : ripple_db;
        double const d = std::pow(10.0, ripple / 10.0); // 纹波深度(功率比)
        // 通带参考点(直流, F = eps 处)的原始电平: 提升时直流取纹波上沿, 衰减时取下沿
        double const a0 = (boost_db > 0.0) ? d : 1.0 / d;
        // b^2: 棚带极值(功率) = C*b^4 = g, C 是直流电平。未修正时 C = 1/a0;
        // 修正后直流处的 F = 0 使 C 只能是 1, 于是 b^4 = g 而不是 g*a0
        double const b2 = even_modify ? std::sqrt(g) : std::sqrt(g * a0);
        double const b = std::sqrt(b2);
        double const eps2 = (1.0 - a0) * b2 / (a0 - b2 * b2);
        double const eps = std::sqrt(eps2 > 0.0 ? eps2 : 0.0);

        // 频率标定 lam: 让中点落在 (1)rad/sec
        double lam = 1.0;
        if (even_modify) {
            // 中点要求 F(1) = 1, 即 eps*Phi(1/lam) = 1; Phi(w) = T_N(sqrt((1-c^2)w^2+c^2)),
            // 反解出 Phi 的宗量 u, 再由 sqrt((1-c^2)/lam^2+c^2) = u 得 lam
            if (eps > 0.0) {
                double const c = std::cos(pi * (static_cast<double>(n) - 1.0) / (2.0 * static_cast<double>(n)));
                // T_N(u) = 1/eps: 1/eps < 1 时 u 落在 [c, 1], 否则落在 [1, inf)
                double const r = 1.0 / eps;
                double const u = (r >= 1.0) ? std::cosh(std::acosh(r) / static_cast<double>(n))
                                            : std::cos(std::acos(r) / static_cast<double>(n));
                lam = std::sqrt((1.0 - c * c) / (u * u - c * c));
            }
        }
        else {
            // 中点的原始电平: |H(1)|^2/a0 应为 sqrt(g)(即 boost/2 dB), 由此解特征值 f
            double const a = std::sqrt(g) * a0;
            double const f2 = (a - 1.0) / (b2 - a / b2);
            double const f = (f2 > 0.0) ? std::sqrt(f2) : 0.0;
            // eps*T_N(1/lam) = f
            if (eps > 0.0 && f > eps) {
                lam = 1.0 / std::cosh(std::acosh(f / eps) / static_cast<double>(n));
            }
        }

        std::array<ZPK, kMaxPrototypePairs> scratch{};
        // 零点: "特征 b*F"的原型极点
        IIRDesign::Chebyshev1(std::span{scratch}.first(num_filter), num_filter, DbOfEps(eps * b), even_modify);
        for (size_t i = 0; i < num_filter; ++i) {
            ret[i].z = scratch[i].p * lam;
        }
        // 极点: "特征 F/b"的原型极点
        IIRDesign::Chebyshev1(std::span{scratch}.first(num_filter), num_filter, DbOfEps(eps / b), even_modify);
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
     * @note 不用偶数阶修正时高频端会落在棚带的低谷处(|H(inf)| = boost - ripple_db dB, 而不是
     *       棚高), 但棚带的极值仍精确等于 beta^2(在有限频率处达到)
     * @note 频率尺度 m 的实现方式是把原型的极点与零点整体乘 m(传输零点一起搬, 所以
     *       两组零点仍然逐点相同; 无穷远处的那个零点保持不动)。未修正时特征函数是
     *       1/(eps*T~_N(m/w)), 直接解 T~_N(m) = 1/(eps2*beta); 修正后特征函数的宗量
     *       多一层 sqrt((1-c^2)(m/w)^2+c^2) 的映射, 要把它反解回 m
     */
    [[nodiscard]] static bool Chebyshev2(std::span<ZPK> ret, size_t num_filter, double boost_db, double ripple_db,
                                         bool even_modify) {
        assert(ret.size() >= num_filter);
        assert(num_filter <= kMaxPrototypePairs);

        if (std::abs(boost_db) < kMinBoostDb) {
            MakeBypass(ret, num_filter);
            return true;
        }
        // 衰减方向 = 提升方向的倒数: 先按 |boost| 设计, 最后交换零点与极点
        bool const invert = boost_db < 0.0;
        double const beta2 = std::pow(10.0, std::abs(boost_db) / 20.0); // 幅度棚高
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
        double const eps1 = beta2 * eps2; // (eps1/eps2)^2 = beta^4

        // 频率尺度 m: 让中点落在 (1)rad/sec, 即在 w = 1 处 F_p = 1/b(见类注释的三电平条件)
        size_t const n = 2 * num_filter;
        double const r = 1.0 / (eps2 * std::sqrt(beta2));
        double u = (r >= 1.0) ? std::cosh(std::acosh(r) / static_cast<double>(n))
                              : std::cos(std::acos(std::min(r, 1.0)) / static_cast<double>(n));
        // 未修正的特征函数是 1/(eps*T~_N(m/w)): 直接 T~_N(m) = r, 所以 m 就是 u。
        // 修正后是 1/(eps*T~_N(sqrt((1-c^2)(m/w)^2+c^2))): u 是那个 sqrt 里的值, 要反解
        double m = u;
        if (even_modify) {
            double const c = std::cos(pi * (static_cast<double>(n) - 1.0) / (2.0 * static_cast<double>(n)));
            u = std::max(u, c);
            m = std::sqrt((u * u - c * c) / (1.0 - c * c));
        }
        if (!(m > 0.0) || !std::isfinite(m)) {
            return false;
        }
        // eps -> Chebyshev2Natural 要的阻带衰减 dB
        auto const rs_of_eps = [](double eps) { return 10.0 * std::log10(eps * eps / (1.0 + eps * eps)); };

        std::array<ZPK, kMaxPrototypePairs> scratch{};
        // 第一批: eps1 的原型, 极点乘 m 后当 H 的极点(衰减方向则当零点)
        IIRDesign::Chebyshev2Natural(std::span{scratch}.first(num_filter), num_filter, rs_of_eps(eps1), even_modify);
        for (size_t i = 0; i < num_filter; ++i) {
            if (invert) {
                ret[i].z = scratch[i].p * m;
            }
            else {
                ret[i].p = scratch[i].p * m;
            }
        }
        // 第二批: eps2 的原型, 极点乘 m 后当 H 的零点(衰减方向则当极点)
        IIRDesign::Chebyshev2Natural(std::span{scratch}.first(num_filter), num_filter, rs_of_eps(eps2), even_modify);
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
     * 形状函数 R(只由阶数与 m1 = eps_p/eps_s 决定)配两个 eps 一起定出三个电平:
     * 直流 0dB、中点 boost/2 dB、棚带极值 boost dB。未修正时高频端本身到不了那个极值
     * (R(inf) 是有限的), 极值出现在棚带的传输零点上; 修正后 R(0) = 0 且 R(inf) = inf,
     * 直流电平就是 1、w -> inf 就是极值。频率尺度取 sqrt(k)(通带边沿 1 与阻带边沿 1/k
     * 的几何中点); 未修正时若该尺度下三个电平不可兼得(低频通带电平会塌成深陷波),
     * 就把尺度推进原型的阻带, 直到低频通带电平偏离 0dB 不超过 3dB。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param boost_db 高频端相对低频端的提升(dB), 可正可负(负值即反向 shelf),
     *        为 0 时输出直通
     * @param ripple_db 通带纹波(dB, >0), 与 stopband_db 一起定出形状函数
     * @param stopband_db 阻带衰减(dB, >ripple_db), 决定过渡带陡度
     * @param even_modify 用偶数阶修正原型: 直流精确 0dB, 且 w -> inf 精确到棚高
     * @return 成功返回 true
     * @retval false 参数退化: 阻带 <= 通带纹波, 或规格过陡(模数在 double 下就是 1),
     *         与 IIRDesign::Elliptic 的失败条件一致
     * @note 三个电平是按参考原型(eps = 1)的形状函数解出来的, 因此会继承 IIRDesign::Elliptic
     *       形状函数的数值误差(阶数越高、阻带越浅越大): 常规规格下偏差在 1e-9 dB 内, 16 阶
     *       配 20dB 阻带约 1e-4 dB, 16 阶配 10dB 阻带可到 0.3dB(此时原型自身已接近退化)
     */
    [[nodiscard]] static bool Elliptic(std::span<ZPK> ret, size_t num_filter, double boost_db, double ripple_db,
                                       double stopband_db, bool even_modify) {
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
        double const k_ref = DesignModulus(num_filter, ripple_db, stopband_db);
        if (!(k_ref > 0.0 && k_ref < 1.0)) {
            return false;
        }
        // 形状函数 R 只由阶数与 m1 决定, 与 eps 的绝对大小无关: 用 eps = 1 的参考原型查值
        std::array<ZPK, kMaxPrototypePairs> ref{};
        if (!IIRDesign::Elliptic(std::span{ref}.first(num_filter), num_filter, DbOfEps(1.0), DbOfEps(1.0 / m1),
                                 even_modify)) {
            return false;
        }
        std::span<ZPK const> const ref_view{ref.data(), num_filter};

        double const amp = std::pow(10.0, boost_db / 20.0);            // 中点目标电平(幅度), 也记作 b
        double const peak = amp * amp;                                 // 棚带极值(功率比), 可正可负
        double scale = std::sqrt(k_ref);                               // 通带边沿(1)与阻带边沿(1/k)的几何中点
        double rho = PrototypeChar(ref_view, num_filter, 1.0 / scale); // R(1/scale)
        if (!(rho > 0.0) || !std::isfinite(rho)) {
            return false;
        }
        double b2;    // b^2: 零点特征函数与极点特征函数的尺度比
        double e_den; // 极点原型的 eps
        if (even_modify) {
            // 修正原型: R(0) = 0 => 直流电平 C = 1, 棚带极值就是 b^4; R(inf) = inf 让它
            // 恰好落在 w -> inf。中点要求 F_p(1) = e_den*R(1/scale) = 1/b
            b2 = amp;
            e_den = 1.0 / (std::sqrt(b2) * rho);
        }
        else {
            // 未修正原型: R(0) = 1, 直流电平 C = (1+e_den^2)/(1+b^4*e_den^2) 由 b 与 e_den 定,
            // 于是"棚带极值 = C*b^4 = peak"与"中点电平 = amp"正好解出这两个量 —— 前提是
            // 尺度给的 R(1/scale) 足够大; 不足时 C 会塌到 0(低频通带变成深陷波), 把 1/scale
            // 推进原型阻带的单调段再解一次
            double const rho_min = RhoForBandDeviation(amp);
            if (rho < rho_min) {
                scale = 1.0 / ShapeInverse(ref_view, num_filter, m1, 1.0 / scale, rho_min);
                rho = PrototypeChar(ref_view, num_filter, 1.0 / scale);
            }
            double const rho2 = rho * rho;
            double const C = (rho2 * (amp - 1.0) - (peak - 1.0)) / (rho2 * (amp - 1.0) - amp * (peak - 1.0) / peak);
            double const e_den2 = (1.0 - C) / (peak - 1.0);
            if (!(C > 0.0) || !(e_den2 > 0.0) || !std::isfinite(rho)) {
                return false;
            }
            b2 = std::sqrt(peak / C);
            e_den = std::sqrt(e_den2);
        }
        double const e_num = b2 * e_den;

        std::array<ZPK, kMaxPrototypePairs> scratch{};
        // 只用原型的极点(当 shelf 的零点/极点用), 传输零点对消后不出现
        // 零点: "特征 e_num"的椭圆原型极点
        if (!IIRDesign::Elliptic(std::span{scratch}.first(num_filter), num_filter, DbOfEps(e_num), DbOfEps(e_num / m1),
                                 even_modify)) {
            return false;
        }
        for (size_t i = 0; i < num_filter; ++i) {
            ret[i].z = scratch[i].p * scale;
        }
        // 极点: "特征 e_den"的椭圆原型极点
        if (!IIRDesign::Elliptic(std::span{scratch}.first(num_filter), num_filter, DbOfEps(e_den), DbOfEps(e_den / m1),
                                 even_modify)) {
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
