#pragma once
#include "iir_design.hpp"

namespace qwqdsp_filter {
class IIRDesignExtra {
public:
    using ZPK = IIRDesign::ZPK;
    static constexpr auto pi = IIRDesign::pi;

    /**
     * @brief 把"截止频率处的线性幅度"换算成 eps^2
     * @param atten 截止频率(1rad/sec)处的线性幅度(0, 1)
     * @return eps^2 = 1 / atten^2 - 1, 使 |H(j*1)| = atten
     */
    static constexpr double AttenGain2SquareEpsi(double atten) noexcept {
        return (1.0f - atten * atten) / (atten * atten);
    }

    /**
     * @brief 把"截止频率处的衰减(dB)"换算成 eps^2
     * @param atten 截止频率(1rad/sec)处的衰减(dB, >0)
     * @return eps^2 = 10^(atten/10) - 1, 使 |H(j*1)| = -atten dB
     */
    static double AttenDb2SquareEpsi(double atten) noexcept {
        return std::pow(10.0, atten / 10.0) - 1.0;
    }

    /**
     * @brief 巴特沃斯原型, 指定 (1)rad/sec 处的线性幅度
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param atten 截止频率处的线性幅度(0, 1), 例如 0.5 就是 -6.02dB
     * @note 零点全在无穷远处, 直流增益 0dB
     * @note -3.01dB 点位于 (1/atten^2 - 1)^(-1/(4*num_filter)) rad/sec
     */
    static void ButterworthAttenGain(std::span<ZPK> ret, size_t num_filter, double atten) {
        return ButterworthAtten(ret, num_filter, AttenGain2SquareEpsi(atten));
    }

    /**
     * @brief 巴特沃斯原型, 指定 (1)rad/sec 处的衰减
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param atten 截止频率处的衰减(dB, >0), 例如 40 就是 -40dB
     * @note 零点全在无穷远处, 直流增益 0dB
     * @note -3.01dB 点位于 (10^(atten/10) - 1)^(-1/(4*num_filter)) rad/sec
     */
    static void ButterworthAttenDb(std::span<ZPK> ret, size_t num_filter, double atten) {
        return ButterworthAtten(ret, num_filter, AttenDb2SquareEpsi(atten));
    }

    /**
     * @brief 切比雪夫 I 型原型, 通带等波纹
     *
     * 通带在 0dB 与 -ripple dB 之间等波纹, (1)rad/sec 处正好是 -atten dB,
     * 也就是说 atten 决定"截止频率落在过渡带的哪一点"。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param ripple 通带纹波(dB, >0), 需要 <= atten
     * @param atten 截止频率(1rad/sec)处的衰减(dB, >0), 需要 >= ripple
     * @param even_pole_modify 偶数阶修正: 把通带参考电平从 -ripple dB 抬到 0dB
     *        (直流增益 1); (1)rad/sec 处仍是 -atten dB
     * @return 所有节增益之积(逐节 k 也已写入 ret, 与 IIRDesign 的约定一致)
     * @note 零点全在无穷远处; 不修正时 section 0 会再除以 sqrt(1 + eps^2)
     * @see IIRDesign::Chebyshev1
     */
    static double Chebyshev1(std::span<ZPK> ret, size_t num_filter, double ripple, double atten,
                             bool even_pole_modify) {
        assert(ret.size() >= num_filter);
        assert(atten >= ripple);

        size_t const n = 2 * num_filter;
        double first_pole = std::cos(pi * (static_cast<double>(n) - 1.0) / (2.0 * static_cast<double>(n)));
        first_pole = first_pole * first_pole;

        double scale = 0.0;
        if (!even_pole_modify) {
            if (atten >= ripple) {
                scale = 1.0
                      / std::cosh(std::acosh(std::sqrt((std::pow(10.0, atten / 10.0) - 1.0)
                                                       / (std::pow(10.0, ripple / 10.0) - 1.0)))
                                  / static_cast<double>(n));
            }
            else {
                scale = 1.0
                      / std::cos(std::acos(std::sqrt((std::pow(10.0, atten / 10.0) - 1.0)
                                                     / (std::pow(10.0, ripple / 10.0) - 1.0)))
                                 / static_cast<double>(n));
            }
        }
        else {
            if (atten >= ripple) {
                scale = std::cosh(
                    std::acosh(std::sqrt((std::pow(10.0, atten / 10.0) - 1.0) / (std::pow(10.0, ripple / 10.0) - 1.0)))
                    / static_cast<double>(n));
            }
            else {
                scale = std::cos(
                    std::acos(std::sqrt((std::pow(10.0, atten / 10.0) - 1.0) / (std::pow(10.0, ripple / 10.0) - 1.0)))
                    / static_cast<double>(n));
            }
            scale *= scale;
            scale = std::sqrt((1.0 - first_pole) / (scale - first_pole));
        }

        double const eps = std::sqrt(std::pow(10.0, ripple / 10.0) - 1.0);
        double const A = 1.0 / static_cast<double>(n) * std::asinh(1.0 / eps);
        double const k_re = std::sinh(A);
        double const k_im = std::cosh(A);

        double gain = 1.0;
        size_t i = 0;
        for (size_t k = 1; k <= num_filter; ++k) {
            double phi = (2.0 * static_cast<double>(k) - 1.0) * pi / (2.0 * static_cast<double>(n));
            if (even_pole_modify) {
                auto pole = std::complex{-std::sin(phi) * k_re, std::cos(phi) * k_im};
                ret[i].p = scale * std::sqrt((pole * pole + first_pole) / (1.0 - first_pole));
            }
            else {
                ret[i].p = scale * std::complex{-std::sin(phi) * k_re, std::cos(phi) * k_im};
            }
            // 与 IIRDesign::Chebyshev1 一致: 逐节 k = |p|^2, 通带跌落因子放到 section 0
            ret[i].k = std::norm(ret[i].p);
            gain *= std::norm(ret[i].p);
            ++i;
        }

        if (!even_pole_modify) {
            ret[0].k /= std::sqrt(1.0f + eps * eps);
            gain /= std::sqrt(1.0f + eps * eps);
        }
        return gain;
    }

    /**
     * @brief 切比雪夫 II 型(逆切比雪夫)原型, 阻带等波纹
     *
     * 通带最平坦(直流增益 0dB), 阻带等波纹深度为 -ripple dB,
     * (1)rad/sec 处正好是 -atten dB, 也就是说 atten 决定"截止频率落在过渡带的哪一点"。
     *
     * @warning 注意两个参数的含义与 Chebyshev1 正好相反:
     *          ripple 是**阻带**涟漪(dB, >0), atten 是截止频率处的幅度(dB, >0)。
     *          习惯写法是 atten < ripple, 例如 atten=3, ripple=40
     *          (atten=3.01, ripple=40 就退化成 IIRDesign::Chebyshev2 的形状)。
     *
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param ripple 阻带等波纹深度(dB, >0)
     * @param atten 截止频率(1rad/sec)处的衰减(dB, >0)
     * @param even_order_modify 偶数阶修正: 让最后一对极点不再有有限零点
     * @return 所有节增益之积(逐节 k 也已写入 ret, 与 IIRDesign 的约定一致)
     * @note 零点在虚轴上(有限频率)
     * @see IIRDesign::Chebyshev2
     */
    static double Chebyshev2(std::span<ZPK> ret, size_t num_filter, double ripple, double atten,
                             bool even_order_modify) {
        assert(ret.size() >= num_filter);

        size_t n = 2 * num_filter;
        double first_pole = std::cos(pi * (static_cast<double>(n) - 1.0) / (2.0 * static_cast<double>(n)));
        first_pole = first_pole * first_pole;

        double scale = 0.0;
        if (!even_order_modify) {
            if (atten < ripple) {
                scale = std::cosh(
                    std::acosh(std::sqrt((std::pow(10.0, ripple / 10.0) - 1.0) / (std::pow(10.0, atten / 10.0) - 1.0)))
                    / static_cast<double>(n));
            }
            else {
                scale = std::cos(
                    std::acos(std::sqrt((std::pow(10.0, ripple / 10.0) - 1.0) / (std::pow(10.0, atten / 10.0) - 1.0)))
                    / static_cast<double>(n));
            }
        }
        else {
            if (atten < ripple) {
                scale = std::cosh(
                    std::acosh(std::sqrt((std::pow(10.0, ripple / 10.0) - 1.0) / (std::pow(10.0, atten / 10.0) - 1.0)))
                    / static_cast<double>(n));
            }
            else {
                scale = std::cos(
                    std::acos(std::sqrt((std::pow(10.0, ripple / 10.0) - 1.0) / (std::pow(10.0, atten / 10.0) - 1.0)))
                    / static_cast<double>(n));
            }
            scale *= scale;
            scale = std::sqrt((scale - first_pole) / (1.0 - first_pole));
        }

        size_t i = 0;
        double eps = 1.0 / std::sqrt(std::pow(10.0, ripple / 10.0) - 1.0);
        double A = 1.0 / static_cast<double>(n) * std::asinh(1.0 / eps);
        double k_re = std::sinh(A);
        double k_im = std::cosh(A);
        double k = 1.0;

        // 最接近0的零点
        double const first_zero = std::cos(
            (static_cast<double>(n) / 2.0 - 1.0 + 0.5) * std::numbers::pi_v<double> / static_cast<double>(n));
        for (size_t kk = 1; kk <= num_filter; ++kk) {
            double phi = (2.0 * static_cast<double>(kk) - 1.0) * pi / (2.0 * static_cast<double>(n));
            if (!even_order_modify) {
                ret[i].z = scale / std::complex{0.0, std::cos(phi)};
                ret[i].p = scale / std::complex{-std::sin(phi) * k_re, std::cos(phi) * k_im};
            }
            else {
                auto pole = std::complex{-std::sin(phi) * k_re, std::cos(phi) * k_im};
                ret[i].p = scale / std::sqrt((pole * pole + first_pole) / (1.0 - first_pole));
                if (kk != num_filter) {
                    // 最靠近0的切比雪夫多项式的零点被映射到0，所以零点在无穷远处不赋值
                    double const zero = std::cos(phi);
                    double const tt = std::sqrt(
                        std::max(0.0, (zero * zero - first_zero * first_zero) / (1.0 - first_zero * first_zero)));
                    ret[i].z = scale / std::complex{0.0, tt};
                }
            }
            if (ret[i].z) {
                ret[i].k = std::norm(ret[i].p) / std::norm(*ret[i].z);
            }
            else {
                ret[i].k = std::norm(ret[i].p);
            }
            k *= ret[i].k;
            ++i;
        }
        return k;
    }
private:
    /**
     * @brief 巴特沃斯原型的共同实现, 直接给 eps^2
     * @param ret 输出的零极点节, 至少 num_filter 个
     * @param num_filter 极点对数, 阶数 = 2 * num_filter
     * @param square_epsi eps^2, 决定 (1)rad/sec 处的幅度: |H(j*1)| = 1/sqrt(1 + eps^2)
     * @note 极点按 g = eps^(-1/(2*num_filter*2)) 缩放(即 -3.01dB 点落在 g 处),
     *       section 0 的 k 取 1/sqrt(eps^2), 使 (1)rad/sec 处正好达到目标幅度
     */
    static void ButterworthAtten(std::span<ZPK> ret, size_t num_filter, double square_epsi) {
        assert(ret.size() >= num_filter);

        double const g = 1.0 / std::pow(square_epsi, 0.25 / static_cast<double>(num_filter));
        size_t const n = 2 * num_filter;
        size_t i = 0;
        for (size_t k = 1; k <= num_filter; ++k) {
            double phi = (2.0 * static_cast<double>(k) - 1.0) * pi / (2.0 * static_cast<double>(n));
            ret[i].p = g * std::complex{-std::sin(phi), std::cos(phi)};
            ret[i].k = 1.0;
            ++i;
        }
        ret[0].k = 1.0 / std::sqrt(square_epsi);
    }
};
} // namespace qwqdsp_filter
