#pragma once

#include <cmath>

namespace qwqdsp_cephes {

namespace detail {

/**
 * @brief 切比雪夫级数求值（单精度）
 */
static inline float chbevlf(float x, const float* array, int n) noexcept {
    float b0, b1, b2;
    const float* p = array;
    b0 = *p++;
    b1 = 0.0f;
    int i = n - 1;
    do {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + *p++;
    } while (--i);
    return 0.5f * (b0 - b2);
}

} // namespace detail

/**
 * @brief 整数阶修正贝塞尔函数（单精度）
 */
struct Besself {
    /**
     * @brief I₀(x) — 零阶修正贝塞尔函数（单精度）
     * @details 相对误差: IEEE 域 [0,30], 峰值 4.0e-7, RMS 7.9e-8
     */
    static float i0(float x) noexcept {
        // Chebyshev coefficients for exp(-x) I0(x) in [0, 8]
        static constexpr float kA[] = {
            -1.30002500998624804212E-8f, 6.04699502254191894932E-8f,  -2.67079385394061173391E-7f,
            1.11738753912010371815E-6f,  -4.41673835845875056359E-6f, 1.64484480707288970893E-5f,
            -5.75419501008210370398E-5f, 1.88502885095841655729E-4f,  -5.76375574538582365885E-4f,
            1.63947561694133579842E-3f,  -4.32430999505057594430E-3f, 1.05464603945949983183E-2f,
            -2.37374148058994688156E-2f, 4.93052842396707084878E-2f,  -9.49010970480476444210E-2f,
            1.71620901522208775349E-1f,  -3.04682672343198398683E-1f, 6.76795274409476084995E-1f,
        };
        // Chebyshev coefficients for exp(-x) sqrt(x) I0(x) in [8, ∞)
        static constexpr float kB[] = {
            3.39623202570838634515E-9f, 2.26666899049817806459E-8f, 2.04891858946906374183E-7f,
            2.89137052083475648297E-6f, 6.88975834691682398426E-5f, 3.36911647825569408990E-3f,
            8.04490411014108831608E-1f,
        };

        float y;
        if (x < 0.0f)
            x = -x;
        if (x <= 8.0f) {
            y = 0.5f * x - 2.0f;
            return std::exp(x) * detail::chbevlf(y, kA, 18);
        }
        return std::exp(x) * detail::chbevlf(32.0f / x - 2.0f, kB, 7) / std::sqrt(x);
    }

    /**
     * @brief I₁(x) — 一阶修正贝塞尔函数（单精度）
     * @details 相对误差: IEEE 域 [0,30], 峰值 1.5e-6, RMS 1.6e-7
     */
    static float i1(float x) noexcept {
        // Chebyshev coefficients for exp(-x) I1(x) / x in [0, 8]
        static constexpr float kA[] = {
            9.38153738649577178388E-9f,  -4.44505912879632808065E-8f, 2.00329475355213526229E-7f,
            -8.56872026469545474066E-7f, 3.47025130813767847674E-6f,  -1.32731636560394358279E-5f,
            4.78156510755005422638E-5f,  -1.61760815825896745588E-4f, 5.12285956168575772895E-4f,
            -1.51357245063125314899E-3f, 4.15642294431288815669E-3f,  -1.05640848946261981558E-2f,
            2.47264490306265168283E-2f,  -5.29459812080949914269E-2f, 1.02643658689847095384E-1f,
            -1.76416518357834055153E-1f, 2.52587186443633654823E-1f,
        };
        // Chebyshev coefficients for exp(-x) sqrt(x) I1(x) in [8, ∞)
        static constexpr float kB[] = {
            -3.83538038596423702205E-9f, -2.63146884688951950684E-8f, -2.51223623787020892529E-7f,
            -3.88256480887769039346E-6f, -1.10588938762623716291E-4f, -9.76109749136146840777E-3f,
            7.78576235018280120474E-1f,
        };

        float y, z;
        z = std::abs(x);
        if (z <= 8.0f) {
            y = 0.5f * z - 2.0f;
            z = detail::chbevlf(y, kA, 17) * z * std::exp(z);
        }
        else {
            z = std::exp(z) * detail::chbevlf(32.0f / z - 2.0f, kB, 7) / std::sqrt(z);
        }
        if (x < 0.0f)
            z = -z;
        return z;
    }
};

} // namespace qwqdsp_cephes
