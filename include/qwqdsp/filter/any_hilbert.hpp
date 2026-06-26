#pragma once
#include <complex>

namespace qwqdsp_filter {
template <typename T>
class AnyHilbert {
public:
    AnyHilbert() noexcept {
        for (int i = 0; i < kState; ++i) {
            auto& f = filters_[i];
            f.a_ = static_cast<T>(a[i]);
            f.b_ = static_cast<T>(b[i]);
            f.c_ = static_cast<T>(c[i]);
            f.d_ = static_cast<T>(d[i]);
        }
    }

    void Reset() noexcept {
        for (auto& f : filters_) {
            f.Reset();
        }
    }

    std::complex<T> Tick(T x) noexcept {
        std::complex<T> y = x * static_cast<T>(direct_);
        for (auto& f : filters_) {
            y += f.Tick(x);
        }
        return y;
    }

    std::complex<T> Tick(std::complex<T> x) noexcept {
        auto y = x * static_cast<T>(direct_);
        for (auto& f : filters_) {
            y += f.Tick(x);
        }
        return y;
    }
private:
    class Filter {
    public:
        void Reset() noexcept {
            sre1_ = T(0);
            sre2_ = T(0);
            sim1_ = T(0);
            sim2_ = T(0);
        }

        std::complex<T> Tick(T x) noexcept {
            T yre = x * a_ + sre1_;
            T yim = sim1_;
            sre1_ = -c_ * yim + sre2_;
            sim1_ = -b_ * x + yre * c_ + sim2_;
            sre2_ = d_ * yre;
            sim2_ = d_ * yim;
            return {yre, yim};
        }

        std::complex<T> Tick(std::complex<T> x) noexcept {
            T xre = x.real();
            T xim = x.imag();
            T yre = xre * a_ + sre1_;
            T yim = xim * a_ + sim1_;
            sre1_ = b_ * xim - c_ * yim + sre2_;
            sim1_ = -b_ * xre + yre * c_ + sim2_;
            sre2_ = d_ * yre;
            sim2_ = d_ * yim;
            return {yre, yim};
        }

        T a_{};
        T b_{};
        T c_{};
        T d_{};
    private:
        T sre1_{};
        T sre2_{};
        T sim1_{};
        T sim2_{};
    };

    // 使用 notebook/any_hilbert.ipynb 导出
    static constexpr double a[16]{
        1.5117691950e-01,  -4.0029448985e-01, 4.9774027954e-01,  -4.4633141833e-01,
        3.2473980264e-01,  -1.9996267924e-01, 1.0336160540e-01,  -4.0935386859e-02,
        6.8327132448e-03,  7.9749493092e-03,  -1.1561746436e-02, 9.7259201265e-03,
        -6.0695498296e-03, 2.6117239910e-03,  -4.4278853955e-04, -1.1149996062e-04,
    };

    static constexpr double b[16]{
        -1.6150960039e-01, 3.1285144583e-01,  -2.3897509472e-01, 9.8686295334e-02,  1.0576186212e-02, -6.5895779752e-02,
        7.8355485382e-02,  -6.6911037420e-02, 4.7011218501e-02,  -2.7827970196e-02, 1.3290298175e-02, -4.0972773302e-03,
        -5.6086795657e-04, 1.9432727706e-03,  -1.3746116533e-03, 3.8518024932e-04,
    };

    static constexpr double c[16]{
        1.3052689240e+00, 1.1135166236e+00, 8.3791368666e-01, 5.8031118445e-01,  3.8228975485e-01, 2.4442612941e-01,
        1.5330819938e-01, 9.4761364871e-02, 5.7737148988e-02, 3.4545671410e-02,  2.0113644317e-02, 1.1189489279e-02,
        5.7275727638e-03, 2.4630086794e-03, 6.3520964134e-04, -1.8305928139e-04,
    };

    static constexpr double d[16]{
        4.5547132049e-01, 5.3524281193e-01, 6.4989752517e-01, 7.5706439789e-01, 8.3944569600e-01, 8.9680198411e-01,
        9.3471366069e-01, 9.5907883314e-01, 9.7449591774e-01, 9.8416742751e-01, 9.9020946329e-01, 9.9398370257e-01,
        9.9635568114e-01, 9.9787510426e-01, 9.9889639359e-01, 9.9965877371e-01,
    };

    static constexpr double direct_ = 0.0016725282372669745;
    static constexpr int kState = 16;
    Filter filters_[kState];
};
} // namespace qwqdsp_filter
