#pragma once
#include <complex>

namespace qwqdsp_filter {
class ChebyHilbert {
public:
    static constexpr int kState = 8;

    ChebyHilbert() noexcept {
        for (int i = 0; i < kState; ++i) {
            auto& f = filters_[i];
            f.a_ = a[i];
            f.b_ = b[i];
            f.c_ = c[i];
            f.d_ = d[i];
        }
    }

    void Reset() noexcept {
        for (auto& f : filters_) {
            f.Reset();
        }
    }

    std::complex<float> Tick(float x) noexcept {
        std::complex<float> y = x * direct_;
        for (auto& f : filters_) {
            y += f.Tick(x);
        }
        return y;
    }
private:
    // class Filter {
    // public:
    //     void Reset() noexcept {
    //         x1_ = 0;
    //         yre1_ = 0;
    //         yre2_ = 0;
    //         yim1_ = 0;
    //         yim2_ = 0;
    //     }

    //     std::complex<float> Tick(float x) noexcept {
    //         float yre = a_ * x - c_ * yim1_ + d_ * yre2_;
    //         float yim = -b_ * x1_ + c_ * yre1_ + d_ * yim2_;

    //         x1_ = x;
    //         yre2_ = yre1_;
    //         yre1_ = yre;
    //         yim2_ = yim1_;
    //         yim1_ = yim;

    //         return {yre, yim};
    //     }

    //     float a_{};
    //     float b_{};
    //     float c_{};
    //     float d_{};
    // private:
    //     float x1_{};
    //     float yre1_{};
    //     float yre2_{};
    //     float yim1_{};
    //     float yim2_{};
    // };

    class Filter {
    public:
        void Reset() noexcept {
            sre1_ = 0;
            sre2_ = 0;
            sim1_ = 0;
            sim2_ = 0;
        }

        std::complex<float> Tick(float x) noexcept {
            float yre = x * a_ + sre1_;
            float yim = sim1_;
            sre1_ = -c_ * yim + sre2_;
            sim1_ = -b_ * x + yre * c_ + sim2_;
            sre2_ = d_ * yre;
            sim2_ = d_ * yim;
            return {yre, yim};
        }

        float a_{};
        float b_{};
        float c_{};
        float d_{};
    private:
        float sre1_{};
        float sre2_{};
        float sim1_{};
        float sim2_{};
    };

    static constexpr float a[8]{-2.3907175666e+02f, -8.2565184991e+01f,
                                4.2140279038e+01f,  4.8262957055e+01f,
                                8.7986278635e+00f,  -5.2244294420e+00f,
                                -1.6655014413e+00f, 5.8373371598e-02f};

    static constexpr float b[8]{5.6924986659e+01f,  4.8269780980e+01f,
                                2.3090074198e+01f,  -5.8144834213e+00f,
                                -1.2236600752e+01f, -3.4469509136e+00f,
                                7.4799516590e-01f,  2.3831392174e-01f};

    static constexpr float c[8]{-4.3091312308e-01f, -3.9923317942e-01f,
                                -3.4107482346e-01f, -2.6546447571e-01f,
                                -1.8325915687e-01f, -1.0549696976e-01f,
                                -4.2721753637e-02f, -5.4974465931e-03f};

    static constexpr float d[8]{5.1700107953e-02f, 8.6329023551e-02f,
                                1.5308354230e-01f, 2.4805334332e-01f,
                                3.6764840562e-01f, 5.1042339746e-01f,
                                6.7881544631e-01f, 8.8114801219e-01f};
    static constexpr float direct_ = 229.2709967;
    Filter filters_[kState];
};
}
