#pragma once
#include <cassert>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>

//    Some software in this archive may be from the book _Methods and
// Programs for Mathematical Functions_ (Prentice-Hall or Simon & Schuster
// International, 1989) or from the Cephes Mathematical Library, a
// commercial product. In either event, it is copyrighted by the author.
// What you see here may be used freely but it comes with no support or
// guarantee.

//    The two known misprints in the book are repaired here in the
// source listings for the gamma function and the incomplete beta
// integral.

//    Stephen L. Moshier
//    moshier@na-net.ornl.gov
namespace cephes {
static constexpr double A[] = {
    -4.41534164647933937950E-18, 3.33079451882223809783E-17, -2.43127984654795469359E-16, 1.71539128555513303061E-15,
    -1.16853328779934516808E-14, 7.67618549860493561688E-14, -4.85644678311192946090E-13, 2.95505266312963983461E-12,
    -1.72682629144155570723E-11, 9.67580903537323691224E-11, -5.18979560163526290666E-10, 2.65982372468238665035E-9,
    -1.30002500998624804212E-8,  6.04699502254191894932E-8,  -2.67079385394061173391E-7,  1.11738753912010371815E-6,
    -4.41673835845875056359E-6,  1.64484480707288970893E-5,  -5.75419501008210370398E-5,  1.88502885095841655729E-4,
    -5.76375574538582365885E-4,  1.63947561694133579842E-3,  -4.32430999505057594430E-3,  1.05464603945949983183E-2,
    -2.37374148058994688156E-2,  4.93052842396707084878E-2,  -9.49010970480476444210E-2,  1.71620901522208775349E-1,
    -3.04682672343198398683E-1,  6.76795274409476084995E-1,
};
static constexpr double B[] = {
    -7.23318048787475395456E-18, -4.83050448594418207126E-18, 4.46562142029675999901E-17,  3.46122286769746109310E-17,
    -2.82762398051658348494E-16, -3.42548561967721913462E-16, 1.77256013305652638360E-15,  3.81168066935262242075E-15,
    -9.55484669882830764870E-15, -4.15056934728722208663E-14, 1.54008621752140982691E-14,  3.85277838274214270114E-13,
    7.18012445138366623367E-13,  -1.79417853150680611778E-12, -1.32158118404477131188E-11, -3.14991652796324136454E-11,
    1.18891471078464383424E-11,  4.94060238822496958910E-10,  3.39623202570838634515E-9,   2.26666899049817806459E-8,
    2.04891858946906374183E-7,   2.89137052083475648297E-6,   6.88975834691682398426E-5,   3.36911647825569408990E-3,
    8.04490411014108831608E-1,
};

static inline double chbevl(double x, const double* array, int N) {
    double b0, b1, b2;
    int i;

    const double* p = array;
    b0 = *p++;
    b1 = 0.0;
    i = N - 1;

    do {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + *p++;
    } while (--i);

    return (0.5 * (b0 - b2));
}

static inline double i0(double x) {
    double y;

    if (x < 0) x = -x;
    if (x <= 8.0) {
        y = (x / 2.0) - 2.0;
        return (std::exp(x) * chbevl(y, A, 30));
    }

    return (std::exp(x) * chbevl(32.0 / x - 2.0, B, 25) / std::sqrt(x));
}

static constexpr double A2[] = {
    2.77791411276104639959E-18, -2.11142121435816608115E-17, 1.55363195773620046921E-16, -1.10559694773538630805E-15,
    7.60068429473540693410E-15, -5.04218550472791168711E-14, 3.22379336594557470981E-13, -1.98397439776494371520E-12,
    1.17361862988909016308E-11, -6.66348972350202774223E-11, 3.62559028155211703701E-10, -1.88724975172282928790E-9,
    9.38153738649577178388E-9,  -4.44505912879632808065E-8,  2.00329475355213526229E-7,  -8.56872026469545474066E-7,
    3.47025130813767847674E-6,  -1.32731636560394358279E-5,  4.78156510755005422638E-5,  -1.61760815825896745588E-4,
    5.12285956168575772895E-4,  -1.51357245063125314899E-3,  4.15642294431288815669E-3,  -1.05640848946261981558E-2,
    2.47264490306265168283E-2,  -5.29459812080949914269E-2,  1.02643658689847095384E-1,  -1.76416518357834055153E-1,
    2.52587186443633654823E-1,
};
static constexpr double B2[] = {
    7.51729631084210481353E-18,  4.41434832307170791151E-18,  -4.65030536848935832153E-17, -3.20952592199342395980E-17,
    2.96262899764595013876E-16,  3.30820231092092828324E-16,  -1.88035477551078244854E-15, -3.81440307243700780478E-15,
    1.04202769841288027642E-14,  4.27244001671195135429E-14,  -2.10154184277266431302E-14, -4.08355111109219731823E-13,
    -7.19855177624590851209E-13, 2.03562854414708950722E-12,  1.41258074366137813316E-11,  3.25260358301548823856E-11,
    -1.89749581235054123450E-11, -5.58974346219658380687E-10, -3.83538038596423702205E-9,  -2.63146884688951950684E-8,
    -2.51223623787020892529E-7,  -3.88256480887769039346E-6,  -1.10588938762623716291E-4,  -9.76109749136146840777E-3,
    7.78576235018280120474E-1,
};

static inline double i1(double x) {
    double y, z;

    z = fabs(x);
    if (z <= 8.0) {
        y = (z / 2.0) - 2.0;
        z = chbevl(y, A2, 29) * z * exp(z);
    }
    else {
        z = exp(z) * chbevl(32.0 / z - 2.0, B2, 25) / sqrt(z);
    }
    if (x < 0.0) {
        z = -z;
    }
    return z;
}
} // namespace cephes

namespace qwqdsp_window {
struct Kaiser {
    // 和分析有关的
    // f = width / N
    static float MainLobeWidth(float beta) noexcept {
        float a = beta / std::numbers::pi_v<float>;
        return std::sqrt(1.0f + a * a);
    }
    static constexpr float kSidelobeRolloff = -6.0f;
    // 和滤波器设计有关的
    // 卷积之后第一个旁瓣的大小
    // qwqfixme: 可能需要补充信息
    // static constexpr float kStopband = -53.0f;
    // static constexpr float kTransmit = 3.3f;

    static void Window(std::span<float> window, float beta, bool for_analyze_not_fir) noexcept {
        const size_t N = window.size();
        if (for_analyze_not_fir) {
            auto down = 1.0f / cephes::i0(beta);
            for (size_t i = 0; i < N; ++i) {
                auto t = static_cast<float>(i) / static_cast<float>(N);
                t = 2 * t - 1;
                auto arg = std::sqrt(1.0f - t * t);
                window[i] = cephes::i0(beta * arg) * down;
            }
        }
        else {
            auto down = 1.0f / cephes::i0(beta);
            for (size_t i = 0; i < N; ++i) {
                auto t = static_cast<float>(i) / (static_cast<float>(N) - 1.0f);
                t = 2 * t - 1;
                auto arg = std::sqrt(1.0f - t * t);
                window[i] = cephes::i0(beta * arg) * down;
            }
        }
    }

    static void ApplyWindow(std::span<float> x, float beta, bool for_analyze_not_fir) noexcept {
        const size_t N = x.size();
        if (for_analyze_not_fir) {
            auto down = 1.0f / cephes::i0(beta);
            for (size_t i = 0; i < N; ++i) {
                auto t = static_cast<float>(i) / static_cast<float>(N);
                t = 2 * t - 1;
                auto arg = std::sqrt(1.0f - t * t);
                x[i] *= cephes::i0(beta * arg) * down;
            }
        }
        else {
            auto down = 1.0f / cephes::i0(beta);
            for (size_t i = 0; i < N; ++i) {
                auto t = static_cast<float>(i) / (static_cast<float>(N) - 1.0f);
                t = 2 * t - 1;
                auto arg = std::sqrt(1.0f - t * t);
                x[i] *= cephes::i0(beta * arg) * down;
            }
        }
    }

    /**
     * @note 此方法只会生成周期性，也就是for analyze
     */
    static void Window(std::span<float> window, std::span<float> dwindow, float beta) noexcept {
        constexpr auto kTimeDelta = 0.001f;
        const size_t N = window.size();

        auto down = 1.0f / cephes::i0(beta);
        for (size_t i = 0; i < N; ++i) {
            auto t = static_cast<float>(i) / static_cast<float>(N);
            t = 2 * t - 1;

            auto arg = std::sqrt(1.0f - t * t);
            window[i] = cephes::i0(beta * arg) * down;
            if (i == 0) {
                dwindow.front() =
                    (cephes::i0(beta * std::sqrt(1.0f - (t + kTimeDelta) * (t + kTimeDelta))) * down - window.front())
                    / kTimeDelta;
            }
            else if (i == N - 1) {
                dwindow.back() =
                    (cephes::i0(beta * std::sqrt(1.0f - (t - kTimeDelta) * (t - kTimeDelta))) * down - window.back())
                    / -kTimeDelta;
            }
            else {
                dwindow[i] = cephes::i1(beta * arg) * beta * (-t / arg) * down;
            }
        }
    }

    /**
     * @param side_lobe >0
     * @ref https://ww2.mathworks.cn/help/signal/ref/kaiser.html
     */
    static float Beta(float side_lobe) noexcept {
        assert(side_lobe > 0);
        if (side_lobe < 21.0f) {
            return 0.0f;
        }
        else if (side_lobe <= 50.0f) {
            return 0.5842f * std::pow(side_lobe - 21.0f, 0.4f) + 0.07886f * (side_lobe - 21.0f);
        }
        else {
            return 0.1102f * (side_lobe - 8.7f);
        }
    }
};
} // namespace qwqdsp_window
