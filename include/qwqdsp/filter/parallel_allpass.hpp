#pragma once
#include "allpass.hpp"
#include "iir_design.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <optional>
#include <vector>

namespace qwqdsp_filter {
namespace detail {
// ------------------------------------------------------------
// 奇数阶模拟低通原型（边沿归一化到 1 rad/sec）
//
// 每个原型给出 1 个实极点 + (order-1)/2 个上半平面的共轭对代表。
// 落到数字域只需先乘 tan(w/2)(双线性预畸变, T=2), 再做 z = (1+s)/(1-s)。
// ------------------------------------------------------------

/**
 * @brief 奇数阶原型的极点
 */
struct OddOrderPrototype {
    double real_pole{};
    std::vector<std::complex<double>> pairs; // 每项代表一对共轭极点(取上半个)
};

/**
 * @brief 巴特沃斯原型，极点在单位圆上，实极点 s = -1
 * @param order 奇数阶
 */
[[nodiscard]] inline OddOrderPrototype ButterworthPrototype(size_t order) {
    auto const n = static_cast<double>(order);
    size_t const half = (order - 1) / 2;
    OddOrderPrototype proto;
    for (size_t k = 0; k <= half; ++k) {
        double const phi = (2.0 * static_cast<double>(k) + 1.0) * std::numbers::pi_v<double> / (2.0 * n);
        if (k == half) {
            proto.real_pole = -1.0;
        }
        else {
            proto.pairs.emplace_back(-std::sin(phi), std::cos(phi));
        }
    }
    return proto;
}

/**
 * @brief 切比雪夫 I 型原型，通带等波纹
 * @param order 奇数阶
 * @param ripple_db 通带纹波(dB, >0)
 */
[[nodiscard]] inline OddOrderPrototype Chebyshev1Prototype(size_t order, double ripple_db) {
    auto const n = static_cast<double>(order);
    double const eps = std::sqrt(std::pow(10.0, ripple_db / 10.0) - 1.0);
    double const gain = std::asinh(1.0 / eps) / n;
    double const sh = std::sinh(gain);
    double const ch = std::cosh(gain);
    size_t const half = (order - 1) / 2;
    OddOrderPrototype proto;
    for (size_t k = 0; k <= half; ++k) {
        double const phi = (2.0 * static_cast<double>(k) + 1.0) * std::numbers::pi_v<double> / (2.0 * n);
        if (k == half) {
            proto.real_pole = -sh;
        }
        else {
            proto.pairs.emplace_back(-std::sin(phi) * sh, std::cos(phi) * ch);
        }
    }
    return proto;
}

/**
 * @brief 切比雪夫 II 型(逆切比雪夫)原型，阻带等波纹，阻带边沿在 (1)rad/sec
 * @param order 奇数阶
 * @param stop_db 阻带衰减(dB, >0)
 * @note 极点 = 切比雪夫 I 型极点的倒数，所以实极点 = -1/sinh(A)
 */
[[nodiscard]] inline OddOrderPrototype Chebyshev2Prototype(size_t order, double stop_db) {
    auto const n = static_cast<double>(order);
    double const eps = 1.0 / std::sqrt(std::pow(10.0, stop_db / 10.0) - 1.0);
    double const gain = std::asinh(1.0 / eps) / n;
    double const sh = std::sinh(gain);
    double const ch = std::cosh(gain);
    size_t const half = (order - 1) / 2;
    OddOrderPrototype proto;
    for (size_t k = 0; k <= half; ++k) {
        double const phi = (2.0 * static_cast<double>(k) + 1.0) * std::numbers::pi_v<double> / (2.0 * n);
        if (k == half) {
            proto.real_pole = -1.0 / sh;
        }
        else {
            auto const q = std::complex{-std::sin(phi) * sh, std::cos(phi) * ch};
            proto.pairs.push_back(std::conj(std::complex{1.0, 0.0} / q));
        }
    }
    return proto;
}

/**
 * @brief 椭圆(考尔)原型，通带阻带都等波纹，通带边沿在 (1)rad/sec
 *
 * 与 IIRDesign::Elliptic 同一套公式(Gray & Markel 的实值形式)，区别只是这里按
 * **奇数阶**取极点：循环到 i = (order-1)/2，此时 arg = 0，sn = 0，公式自然退化成
 * 实极点 -sn1/cn1，无需特判。模数 k 与 u 仍用 cephes 的 Ellpk/Cay/Ellik/Ellpj。
 *
 * @param order 奇数阶
 * @param pass_db 通带纹波(dB, >0)
 * @param stop_db 阻带衰减(dB, >0)，需要大于 pass_db
 * @return 规格非法时返回 nullopt(与 IIRDesign::Elliptic 的失败条件一致)
 */
[[nodiscard]] inline std::optional<OddOrderPrototype> EllipticPrototype(size_t order, double pass_db,
                                                                        double stop_db) {
    double const eps_passband = std::sqrt(std::pow(10.0, pass_db / 10.0) - 1.0);
    double const eps_stopband = std::sqrt(std::pow(10.0, stop_db / 10.0) - 1.0);
    if (!(eps_stopband > eps_passband)) {
        return std::nullopt;
    }
    auto const n = static_cast<double>(order);
    double const m1 = eps_passband / eps_stopband;
    double const m1_sq = m1 * m1;
    double const Kk1 = IIRDesign::EllipticHelperCephes::Ellpk(1.0 - m1_sq);
    double const Kpk1 = IIRDesign::EllipticHelperCephes::Ellpk(m1_sq);
    double const q = std::exp(-std::numbers::pi_v<double> * Kpk1 / (n * Kk1));
    double const k = IIRDesign::EllipticHelperCephes::Cay(q);
    if (!(k > 0.0 && k < 1.0)) {
        return std::nullopt;
    }
    double const m = k * k;
    double const Kk = IIRDesign::EllipticHelperCephes::Ellpk(1.0 - m);
    double const u = IIRDesign::EllipticHelperCephes::Ellik(std::atan(1.0 / eps_passband), 1.0 - m1_sq) * Kk
                   / (n * Kk1);
    double sn1 = 0.0;
    double cn1 = 0.0;
    double dn1 = 0.0;
    IIRDesign::EllipticHelperCephes::Ellpj(u, 1.0 - m, sn1, cn1, dn1);

    OddOrderPrototype proto;
    size_t const half = (order - 1) / 2;
    for (size_t i = 0; i <= half; ++i) {
        double const arg = static_cast<double>(order - 1 - 2 * i) * Kk / n;
        double sn = 0.0;
        double cn = 0.0;
        double dn = 0.0;
        IIRDesign::EllipticHelperCephes::Ellpj(arg, m, sn, cn, dn);
        double const r = k * sn * sn1;
        double const den = cn1 * cn1 + r * r;
        auto const pole = std::complex{-cn * dn * sn1 * cn1 / den, sn * dn1 / den};
        if (i == half) {
            // arg = 0 -> sn = 0，公式退化为实极点
            proto.real_pole = pole.real();
        }
        else {
            proto.pairs.push_back(pole);
        }
    }
    return proto;
}
} // namespace detail

/**
 * @brief 双路全通滤波器组成的低通/高通滤波器，只能为奇数阶
 *
 * H(z) = 1/2 (A0(z) + A1(z)) 为低通，1/2 (A0(z) - A1(z)) 为功率互补高通。
 * 极点按模长 |p| 升序交替分配到两条链（实测对 butter/cheby1/cheby2/ellip 成立，
 * 注意**不是**按幅角交替）。系数全部由模拟原型公式 + 双线性变换直接生成，不经
 * 求根，窄带高阶时比先求 A 的根再回代更稳。
 *
 * @note up + down = LP, up - down = HP
 * @ref 纯极点 https://radiosystemdesign.com/assets/pdf/downloads/Reducing_IIR_Comp_Workload_Lyons.pdf
 * @ref 可零点
 * https://www.researchgate.net/publication/278320928_A_Most_Efficient_Digital_Filter_The_Two-Path_Recursive_All-Pass_Filter
 * @ref 分解方法与数值验证见 labs/tf2ca
 */
class ParallelAllpass {
public:
    void Reset() noexcept {}

    /**
     * @brief 巴特沃斯低通，-3.01dB 点在 w
     * @param order 奇数阶
     * @param w 数字截止频率(rad/sample, 0 < w < pi)
     */
    void BuildButterworth(size_t order, float w) {
        Build(order, w, detail::ButterworthPrototype(order));
    }

    /**
     * @brief 切比雪夫 I 型低通，通带边沿在 w
     * @param order 奇数阶
     * @param w 数字通带边沿(rad/sample)
     * @param ripple_db 通带纹波(dB, >0)
     */
    void BuildChebyshev1(size_t order, float w, float ripple_db) {
        Build(order, w, detail::Chebyshev1Prototype(order, static_cast<double>(ripple_db)));
    }

    /**
     * @brief 切比雪夫 II 型(逆切比雪夫)低通，阻带边沿在 w
     * @param order 奇数阶
     * @param w 数字阻带边沿(rad/sample)
     * @param stop_db 阻带衰减(dB, >0)
     */
    void BuildChebyshev2(size_t order, float w, float stop_db) {
        Build(order, w, detail::Chebyshev2Prototype(order, static_cast<double>(stop_db)));
    }

    /**
     * @brief 椭圆低通，通带边沿在 w（阻带边沿由阶数与两个纹波决定）
     * @param order 奇数阶
     * @param w 数字通带边沿(rad/sample)
     * @param pass_db 通带纹波(dB, >0)
     * @param stop_db 阻带衰减(dB, >0)
     * @return 规格过陡导致模数无效时返回 false，此时不修改滤波器
     */
    [[nodiscard]] bool BuildElliptic(size_t order, float w, float pass_db, float stop_db) {
        auto proto = detail::EllipticPrototype(order, static_cast<double>(pass_db), static_cast<double>(stop_db));
        if (!proto) {
            return false;
        }
        Build(order, w, *proto);
        return true;
    }

    /**
     * @note up + down = LP, up - down = HP
     * @return {up, down}
     */
    std::pair<float, float> Tick(float x) noexcept {
        // up
        float up = x;
        up = allpass1_.Tick(up);
        for (size_t i = 0; i < num_up2_; ++i) {
            up = allpass_[i].Tick(up);
        }
        // down
        float down = x;
        for (size_t i = 0; i < num_down2_; ++i) {
            down = allpass_[i + num_up2_].Tick(down);
        }
        return {up, down};
    }

    /**
     * @note up + down = LP, up - down = HP
     * @return {up, down}
     */
    std::pair<float, float> Tick(float up, float down) noexcept {
        // up
        up = allpass1_.Tick(up);
        for (size_t i = 0; i < num_up2_; ++i) {
            up = allpass_[i].Tick(up);
        }
        // down
        for (size_t i = 0; i < num_down2_; ++i) {
            down = allpass_[i + num_up2_].Tick(down);
        }
        return {up, down};
    }

    /**
     * @return {up, down}
     */
    std::pair<std::complex<float>, std::complex<float>> GetResponce(std::complex<float> z) noexcept {
        // up
        std::complex<float> up = allpass1_.GetResponce(z);
        for (size_t i = 0; i < num_up2_; ++i) {
            up *= allpass_[i].GetResponce(z);
        }
        // down
        std::complex<float> down = 1.0f;
        for (size_t i = 0; i < num_down2_; ++i) {
            down *= allpass_[i + num_up2_].GetResponce(z);
        }
        return {up, down};
    }
private:
    /**
     * @brief 由奇数阶原型装配两条全通链
     *
     * 数字极点 z = (1 + s*tan(w/2)) / (1 - s*tan(w/2))，按 |z| 升序交替分配；
     * 交换两条链不影响和，所以固定让实极点所在链为 up。
     */
    void Build(size_t order, float w, detail::OddOrderPrototype const& proto) {
        assert(order % 2 == 1);
        assert(proto.pairs.size() + 1 == (order + 1) / 2);

        double const t = std::tan(static_cast<double>(w) / 2.0);
        auto to_digital = [t](std::complex<double> s) {
            s *= t;
            return (1.0 + s) / (1.0 - s);
        };

        struct Node {
            double radius;
            std::complex<double> pole;
            bool first_order;
        };
        std::vector<Node> nodes;
        nodes.reserve(proto.pairs.size() + 1);
        {
            auto const z = to_digital(proto.real_pole);
            nodes.push_back({std::abs(z), z, true});
        }
        for (auto const& p : proto.pairs) {
            auto const z = to_digital(p);
            nodes.push_back({std::abs(z), z, false});
        }
        std::sort(nodes.begin(), nodes.end(), [](Node const& a, Node const& b) {
            if (a.radius != b.radius) {
                return a.radius < b.radius;
            }
            return a.first_order && !b.first_order;
        });

        size_t real_rank = 0;
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].first_order) {
                real_rank = i;
                break;
            }
        }

        std::vector<std::array<float, 2>> up;
        std::vector<std::array<float, 2>> down;
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].first_order) {
                continue;
            }
            auto& target = ((i % 2) == (real_rank % 2)) ? up : down;
            target.push_back({
                static_cast<float>(-2.0 * nodes[i].pole.real()),
                static_cast<float>(std::norm(nodes[i].pole)),
            });
        }

        allpass1_.SetA(static_cast<float>(-nodes[real_rank].pole.real()));
        num_up2_ = up.size();
        num_down2_ = down.size();
        allpass_.resize(num_up2_ + num_down2_);
        for (size_t i = 0; i < num_up2_; ++i) {
            allpass_[i].SetA1(up[i][0]);
            allpass_[i].SetA2(up[i][1]);
        }
        for (size_t i = 0; i < num_down2_; ++i) {
            allpass_[num_up2_ + i].SetA1(down[i][0]);
            allpass_[num_up2_ + i].SetA2(down[i][1]);
        }
    }

    size_t num_up2_{};
    size_t num_down2_{};
    std::vector<qwqdsp_filter::AllpassOrder2> allpass_;
    qwqdsp_filter::AllpassOrder1 allpass1_;
};
} // namespace qwqdsp_filter
