#pragma once
#include <complex>

namespace qwqdsp_filter {
class ChebyHilbert {
public:
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

    std::complex<float> Tick(std::complex<float> x) noexcept {
        auto y = x * direct_;
        for (auto& f : filters_) {
            y += f.Tick(x);
        }
        return y;
    }
private:
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

        std::complex<float> Tick(std::complex<float> x) noexcept {
            float xre = x.real();
            float xim = x.imag();
            float yre = xre * a_ + sre1_;
            float yim = xim * a_ + sim1_;
            sre1_ = b_ * xim - c_ * yim + sre2_;
            sim1_ = -b_ * xre + yre * c_ + sim2_;
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

    static constexpr float a[10]{1.6574956863e-01f,  -4.4946462011e-01f,
                                 5.1626530259e-01f,  -3.9732789152e-01f,
                                 2.3194735135e-01f,  -1.0092960601e-01f,
                                 2.3796249319e-02f,  8.1745353265e-03f,
                                 -1.1271962217e-02f, 3.8697627207e-03f};

    static constexpr float b[10]{-2.4695571071e-01f, 3.9926539486e-01f,
                                 -2.3138154268e-01f, 4.3887683681e-02f,
                                 5.8384540151e-02f,  -8.2074422876e-02f,
                                 6.3154090954e-02f,  -3.2417582617e-02f,
                                 8.6884116301e-03f,  -1.1365023712e-04f};

    static constexpr float c[10]{1.1786887458e+00f, 9.4547600008e-01f,
                                 6.4253963534e-01f, 3.9284635087e-01f,
                                 2.2382264913e-01f, 1.1968608563e-01f,
                                 5.8667826590e-02f, 2.4256723221e-02f,
                                 6.0386857175e-03f, -1.8047429878e-03f};

    static constexpr float d[10]{3.8609229860e-01f, 5.0509252022e-01f,
                                 6.5969622012e-01f, 7.8719096750e-01f,
                                 8.7362328210e-01f, 9.2711452523e-01f,
                                 9.5889624968e-01f, 9.7761400992e-01f,
                                 9.8898685400e-01f, 9.9669232049e-01f};

    static constexpr float direct_ = 0.00994509f;
    static constexpr int kState = 10;
    Filter filters_[kState];
};
}
