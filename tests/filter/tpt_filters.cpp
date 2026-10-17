// TPT 滤波器族的频响测试。
//
// 每个滤波器输入单位脉冲，做 FFT 得到数字域 |H(e^jw)|；
// 理想响应取模拟原型在双线性变换后的频率处求值：
//   s = j * tan(w/2) / g
// 其中 g 是该滤波器积分器使用的预扭曲截止增益。失败返回非零退出码。
#include <qwqdsp/filter/iir_design.hpp>
#include <qwqdsp/filter/tpt/ladder.hpp>
#include <qwqdsp/filter/tpt/ladder_8pole.hpp>
#include <qwqdsp/filter/tpt/one_pole_tpt.hpp>
#include <qwqdsp/filter/tpt/onepole_tpt_shelf.hpp>
#include <qwqdsp/filter/tpt/ota_one_pole.hpp>
#include <qwqdsp/filter/tpt/svf_tpt.hpp>
#include <qwqdsp/filter/tpt/svf_tpt_bandpass.hpp>
#include <qwqdsp/filter/tpt/svf_tpt_shelf.hpp>
#include <qwqdsp/filter/tpt/transpose_sallen_key.hpp>
#include <qwqdsp/filter/tpt/transpose_svf_tpt.hpp>
#include <qwqdsp/filter/tpt/transpose_svf_tpt_bandpass.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <iostream>
#include <numbers>
#include <string_view>
#include <vector>

namespace {
using Complex = std::complex<double>;

constexpr double kPi = std::numbers::pi_v<double>;
constexpr double kTwoPi = 2.0 * kPi;
constexpr size_t kFftSize = 8192;

int g_failures = 0;

void Check(bool ok, std::string_view what, double got, double limit) {
    if (!ok) {
        ++g_failures;
        std::cout << std::format("  FAIL {}: got {:.3e}, limit {:.3e}\n", what, got, limit);
    }
}

void Fft(std::vector<Complex>& a) {
    size_t const n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    for (size_t len = 2; len <= n; len <<= 1) {
        double const angle = -kTwoPi / static_cast<double>(len);
        Complex const wlen{std::cos(angle), std::sin(angle)};
        for (size_t i = 0; i < n; i += len) {
            Complex w{1.0, 0.0};
            for (size_t j = 0; j < len / 2; ++j) {
                Complex const u = a[i + j];
                Complex const v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

Complex S(double w, double g) {
    return {0.0, std::tan(w / 2.0) / g};
}

template <typename Tick, typename Ideal>
void CheckResponse(std::string_view tag, Tick&& tick, Ideal&& ideal, double amplitude = 1.0,
                   double rel_abs_tol = 1e-4, double rel_tol = 2e-3, bool check_complex = true) {
    std::vector<Complex> spectrum(kFftSize);
    for (size_t n = 0; n < kFftSize; ++n) {
        double const x = n == 0 ? amplitude : 0.0;
        spectrum[n] = static_cast<double>(tick(x));
    }
    Fft(spectrum);

    double max_ideal = 0.0;
    for (size_t k = 1; k < kFftSize / 2; ++k) {
        double const w = kTwoPi * static_cast<double>(k) / static_cast<double>(kFftSize);
        max_ideal = std::max(max_ideal, std::abs(ideal(w) * amplitude));
    }

    double max_abs = 0.0;
    double max_rel = 0.0;
    double max_complex = 0.0;
    for (size_t k = 1; k < kFftSize / 2; ++k) {
        double const w = kTwoPi * static_cast<double>(k) / static_cast<double>(kFftSize);
        Complex const expected = ideal(w) * amplitude;
        double const abs_err = std::abs(std::abs(spectrum[k]) - std::abs(expected));
        max_abs = std::max(max_abs, abs_err);
        if (std::abs(expected) > 1e-3 * max_ideal) {
            max_rel = std::max(max_rel, abs_err / std::abs(expected));
        }
        if (check_complex) {
            max_complex = std::max(max_complex, std::abs(spectrum[k] - expected));
        }
    }

    double const scale = std::max(1.0, max_ideal);
    Check(max_abs <= rel_abs_tol * scale, std::format("{}: |H| abs error", tag), max_abs / scale, rel_abs_tol);
    Check(max_rel <= rel_tol, std::format("{}: |H| rel error", tag), max_rel, rel_tol);
    if (check_complex) {
        Check(max_complex <= rel_tol * scale, std::format("{}: complex error", tag), max_complex / scale, rel_tol);
    }
}

void CheckOnePole() {
    constexpr double w0 = 0.2 * kPi;
    double const g = std::tan(w0 / 2.0);

    {
        qwqdsp_filter::OnePoleTPT f;
        f.MakeLowpass(static_cast<float>(w0));
        CheckResponse("OnePoleTPT LP", [&](double x) { return f.TickLowpass(static_cast<float>(x)); },
                      [&](double w) { return 1.0 / (1.0 + S(w, g)); });
    }
    {
        qwqdsp_filter::OnePoleTPT f;
        f.MakeLowpass(static_cast<float>(w0));
        CheckResponse("OnePoleTPT HP", [&](double x) { return f.TickHighpass(static_cast<float>(x)); },
                      [&](double w) { Complex const s = S(w, g); return s / (1.0 + s); });
    }
    {
        qwqdsp_filter::OnePoleTPT f;
        f.MakeLowpass(static_cast<float>(w0));
        CheckResponse("OnePoleTPT AP", [&](double x) { return f.TickAllpass(static_cast<float>(x)); },
                      [&](double w) { Complex const s = S(w, g); return (1.0 - s) / (1.0 + s); });
    }
}

void CheckOnePoleShelf() {
    constexpr double w0 = 0.2 * kPi;
    constexpr double db = 12.0;
    double const m = std::pow(10.0, db / 40.0);
    double const g = std::tan(w0 / 2.0) * m;
    double const m2 = m * m;
    double const invm2 = 1.0 / m2;

    qwqdsp_filter::OnepoleTPTShelf f;
    f.Set(static_cast<float>(w0), static_cast<float>(m));

    f.Reset();
    CheckResponse("OnepoleTPTShelf low", [&](double x) { return f.TickLowshelf(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const l = 1.0 / (1.0 + S(w, g));
                      return 1.0 + (invm2 - 1.0) * l;
                  });
    f.Reset();
    CheckResponse("OnepoleTPTShelf high", [&](double x) { return f.TickHighshelf(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const l = 1.0 / (1.0 + S(w, g));
                      return m2 + (1.0 - m2) * l;
                  });
    f.Reset();
    CheckResponse("OnepoleTPTShelf tilt", [&](double x) { return f.TickTiltshelf(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const l = 1.0 / (1.0 + S(w, g));
                      return m + (1.0 / m - m) * l;
                  });
}

void CheckOtaOnePole() {
    constexpr double w0 = 0.2 * kPi;
    double const g = std::tan(w0 / 2.0);

    qwqdsp_filter::OTAOnePole f;
    f.Set(static_cast<float>(w0));
    CheckResponse("OTAOnePole LP (small signal)", [&](double x) { return f.Tick(static_cast<float>(x)); },
                  [&](double w) { return 1.0 / (1.0 + S(w, g)); }, 1e-3, 1e-4, 5e-3, false);
}

void CheckSvf() {
    constexpr double w0 = 0.2 * kPi;
    constexpr double q = 2.0;
    double const g = std::tan(w0 / 2.0);
    double const r2 = 1.0 / q;

    qwqdsp_filter::SvfTPT f;
    f.SetCoeffQ(static_cast<float>(w0), static_cast<float>(q));

    f.Reset();
    CheckResponse("SvfTPT LP", [&](double x) { return f.TickLowpass(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const s = S(w, g);
                      return 1.0 / (s * s + r2 * s + 1.0);
                  });
    f.Reset();
    CheckResponse("SvfTPT HP", [&](double x) { return f.TickHighpass(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const s = S(w, g);
                      return s * s / (s * s + r2 * s + 1.0);
                  });
    f.Reset();
    CheckResponse("SvfTPT BP", [&](double x) { return f.TickBandpass(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const s = S(w, g);
                      return s / (s * s + r2 * s + 1.0);
                  });
    f.Reset();
    CheckResponse("SvfTPT norm BP", [&](double x) { return f.TickNormBandpass(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const s = S(w, g);
                      return r2 * s / (s * s + r2 * s + 1.0);
                  });
    f.Reset();
    CheckResponse("SvfTPT notch", [&](double x) { return f.TickNotch(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const s = S(w, g);
                      return (s * s + 1.0) / (s * s + r2 * s + 1.0);
                  });
    f.Reset();
    CheckResponse("SvfTPT allpass", [&](double x) { return f.TickAllpass(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const s = S(w, g);
                      return (s * s - r2 * s + 1.0) / (s * s + r2 * s + 1.0);
                  });
}

void CheckSvfShelf() {
    constexpr double w0 = 0.2 * kPi;
    constexpr double q = 0.70710678;
    constexpr double db = 12.0;
    double const m = std::pow(10.0, db / 80.0);
    double const m2 = m * m;
    double const invm2 = 1.0 / m2;
    double const g = std::tan(w0 / 2.0) * m;
    double const r2 = 1.0 / q;

    qwqdsp_filter::SvfTPTShelf f;
    f.SetCoeffQ(static_cast<float>(w0), static_cast<float>(q), static_cast<float>(db));

    f.Reset();
    CheckResponse("SvfTPTShelf low", [&](double x) { return f.TickLowshelf(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const s = S(w, g);
                      Complex const d = s * s + r2 * s + 1.0;
                      Complex const hp = s * s / d;
                      Complex const bp = s / d;
                      Complex const lp = 1.0 / d;
                      return lp * invm2 * invm2 + r2 * bp * invm2 + hp;
                  });
    f.Reset();
    CheckResponse("SvfTPTShelf high", [&](double x) { return f.TickHighshelf(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const s = S(w, g);
                      Complex const d = s * s + r2 * s + 1.0;
                      Complex const hp = s * s / d;
                      Complex const bp = s / d;
                      Complex const lp = 1.0 / d;
                      return lp + r2 * bp * m2 + hp * m2 * m2;
                  });
    f.Reset();
    CheckResponse("SvfTPTShelf tilt", [&](double x) { return f.TickTiltshelf(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const s = S(w, g);
                      Complex const d = s * s + r2 * s + 1.0;
                      Complex const hp = s * s / d;
                      Complex const bp = s / d;
                      Complex const lp = 1.0 / d;
                      return lp * invm2 + r2 * bp + hp * m2;
                  });
}

void CheckTransposeSallenKey() {
    constexpr double w0 = 0.2 * kPi;
    constexpr double k = 1.0;
    double const g = std::tan(w0 / 2.0);

    qwqdsp_filter::TransposeSallenKey f;
    f.Set(static_cast<float>(w0), static_cast<float>(k));

    f.Reset();
    CheckResponse("TransposeSallenKey LP", [&](double x) { return f.Tick(static_cast<float>(x)).lp2; },
                  [&](double w) {
                      Complex const s = S(w, g);
                      return 1.0 / (s * s + (2.0 - k) * s + 1.0);
                  });
    f.Reset();
    CheckResponse("TransposeSallenKey BP", [&](double x) { return f.Tick(static_cast<float>(x)).bp2; },
                  [&](double w) {
                      Complex const s = S(w, g);
                      return s / (s * s + (2.0 - k) * s + 1.0);
                  });
    f.Reset();
    CheckResponse("TransposeSallenKey HP", [&](double x) { return f.Tick(static_cast<float>(x)).hp2; },
                  [&](double w) {
                      Complex const s = S(w, g);
                      return s * s / (s * s + (2.0 - k) * s + 1.0);
                  });
}

void CheckTransposeSvf() {
    constexpr double w0 = 0.2 * kPi;
    constexpr double q = 2.0;
    double const g = std::tan(w0 / 2.0);
    double const r2 = 1.0 / q;

    auto ideal = [&](double w, double hp_mix, double bp_mix, double lp_mix) {
        Complex const s = S(w, g);
        return (hp_mix * s * s + bp_mix * s + lp_mix) / (s * s + r2 * s + 1.0);
    };

    qwqdsp_filter::TransposeSvfTPT f;
    f.SetCoeffQ(static_cast<float>(w0), static_cast<float>(q));

    f.SetInputMix(0.0f, 0.0f, 1.0f);
    f.Reset();
    CheckResponse("TransposeSvfTPT LP", [&](double x) { return f.Tick(static_cast<float>(x)); },
                  [&](double w) { return ideal(w, 0.0, 0.0, 1.0); });

    f.SetInputMix(1.0f, 0.0f, 0.0f);
    f.Reset();
    CheckResponse("TransposeSvfTPT HP", [&](double x) { return f.Tick(static_cast<float>(x)); },
                  [&](double w) { return ideal(w, 1.0, 0.0, 0.0); });

    f.SetInputMix(0.0f, 1.0f, 0.0f);
    f.Reset();
    CheckResponse("TransposeSvfTPT BP", [&](double x) { return f.Tick(static_cast<float>(x)); },
                  [&](double w) { return ideal(w, 0.0, 1.0, 0.0); });

    f.SetInputMix(1.0f, 0.0f, 1.0f);
    f.Reset();
    CheckResponse("TransposeSvfTPT notch", [&](double x) { return f.Tick(static_cast<float>(x)); },
                  [&](double w) { return ideal(w, 1.0, 0.0, 1.0); });

    f.SetInputMix(1.0f, static_cast<float>(-r2), 1.0f);
    f.Reset();
    CheckResponse("TransposeSvfTPT allpass", [&](double x) { return f.Tick(static_cast<float>(x)); },
                  [&](double w) { return ideal(w, 1.0, -r2, 1.0); });

    f.SetInputMix(0.5f, 0.3f, 1.2f);
    f.Reset();
    CheckResponse("TransposeSvfTPT mix", [&](double x) { return f.Tick(static_cast<float>(x)); },
                  [&](double w) { return ideal(w, 0.5, 0.3, 1.2); });

    // 文章里的极端混合系数：LP 系数约 2.17e4，验证转置结构的系数合并。
    constexpr double w_low = 0.01 * kPi;
    double const g_low = std::tan(w_low / 2.0);
    f.SetCoeffQ(static_cast<float>(w_low), static_cast<float>(0.70710678));
    f.SetInputMix(0.0275746f, 0.0f, 21690.9f);
    f.Reset();
    CheckResponse("TransposeSvfTPT extreme mix",
                  [&](double x) { return f.Tick(static_cast<float>(x)); },
                  [&](double w) {
                      Complex const s = S(w, g_low);
                      double const r2_low = 1.0 / 0.70710678;
                      return (0.0275746 * s * s + 21690.9) / (s * s + r2_low * s + 1.0);
                  }, 1.0, 2e-4, 5e-3);
}

void CheckTransposeSvfMultiInput() {
    constexpr float w0 = 0.25f * static_cast<float>(kPi);
    constexpr float q = 1.3f;

    qwqdsp_filter::TransposeSvfTPT multi;
    multi.SetCoeffQ(w0, q);

    qwqdsp_filter::SvfTPT ref_hp;
    qwqdsp_filter::SvfTPT ref_bp;
    qwqdsp_filter::SvfTPT ref_lp;
    ref_hp.SetCoeffQ(w0, q);
    ref_bp.SetCoeffQ(w0, q);
    ref_lp.SetCoeffQ(w0, q);

    float max_diff = 0.0f;
    for (size_t n = 0; n < kFftSize; ++n) {
        float const fn = static_cast<float>(n);
        float const x_hp = std::sin(0.010f * fn);
        float const x_bp = std::cos(0.023f * fn);
        float const x_lp = std::sin(0.031f * fn);

        float const got = multi.Tick(x_hp, x_bp, x_lp);
        float const expected = ref_hp.TickHighpass(x_hp)
                             + ref_bp.TickBandpass(x_bp)
                             + ref_lp.TickLowpass(x_lp);
        max_diff = std::max(max_diff, std::abs(got - expected));
    }
    Check(max_diff <= 1e-5f, "TransposeSvfTPT three-input", static_cast<double>(max_diff), 1e-5);
}

void CheckSvfBandpass() {
    using ZPK = qwqdsp_filter::IIRDesign::ZPK;

    std::array<ZPK, 1> prototype{};
    bool const ok = qwqdsp_filter::IIRDesign::Elliptic(prototype, 1, 0.5, 40.0, false);
    Check(ok, "SvfTPTBandpass prototype design", ok ? 1.0 : 0.0, 1.0);
    Check(prototype[0].z.has_value(), "SvfTPTBandpass prototype imaginary zero",
          prototype[0].z.has_value() ? 1.0 : 0.0, 1.0);
    if (!ok || !prototype[0].z) {
        return;
    }

    ZPK const& zpk = prototype[0];

    auto run = [&](std::string_view tag, float w1, float w2) {
        double const aw1 = std::tan(static_cast<double>(w1) / 2.0);
        double const aw2 = std::tan(static_cast<double>(w2) / 2.0);
        double const bw = aw2 - aw1;
        double const w0sq = aw1 * aw2;

        auto ideal = [&](double w) {
            double const omega = std::tan(w / 2.0);
            double const u = (omega * omega - w0sq) / (bw * omega);
            Complex const s{0.0, u};
            return zpk.k * (s * s + std::norm(*zpk.z))
                 / (s * s - 2.0 * zpk.p.real() * s + std::norm(zpk.p));
        };

        qwqdsp_filter::SvfTPTBandpass direct(zpk.p, *zpk.z, zpk.k);
        direct.Set(w1, w2);
        CheckResponse(std::format("{} direct", tag),
                      [&](double x) { return direct.Tick(static_cast<float>(x)); }, ideal);

        qwqdsp_filter::TransposeSvfTPTBandpass transposed(zpk.p, *zpk.z, zpk.k);
        transposed.Set(w1, w2);
        CheckResponse(std::format("{} transposed", tag),
                      [&](double x) { return transposed.Tick(static_cast<float>(x)); }, ideal);
    };

    run("SvfTPTBandpass 0.20-0.35pi", 0.20f * static_cast<float>(kPi), 0.35f * static_cast<float>(kPi));
    run("SvfTPTBandpass 0.10-0.13pi", 0.10f * static_cast<float>(kPi), 0.13f * static_cast<float>(kPi));
    run("SvfTPTBandpass 0.200-0.205pi", 0.200f * static_cast<float>(kPi), 0.205f * static_cast<float>(kPi));

    {
        constexpr float w1 = 0.20f * static_cast<float>(kPi);
        constexpr float w2 = 0.35f * static_cast<float>(kPi);
        qwqdsp_filter::SvfTPTBandpass direct(zpk.p, *zpk.z, zpk.k);
        qwqdsp_filter::TransposeSvfTPTBandpass transposed(zpk.p, *zpk.z, zpk.k);
        direct.Set(w1, w2);
        transposed.Set(w1, w2);

        float max_diff = 0.0f;
        for (size_t n = 0; n < kFftSize; ++n) {
            float const x = n == 0 ? 1.0f : 0.0f;
            max_diff = std::max(max_diff, std::abs(direct.Tick(x) - transposed.Tick(x)));
        }
        Check(max_diff <= 1e-5f, "SvfTPTBandpass direct vs transposed",
              static_cast<double>(max_diff), 1e-5);
    }
}

void CheckLadder() {
    constexpr double w0 = 0.2 * kPi;
    constexpr double k = 1.0;
    constexpr double q = 0.70710678;
    double const g = std::tan(w0 / 2.0);
    double const r2 = 1.0 / q;

    auto ideal = [&](double w, int mode) {
        Complex const s = S(w, g);
        Complex const one_plus_s = 1.0 + s;
        Complex const d4 = one_plus_s * one_plus_s * one_plus_s * one_plus_s + k;
        if (mode == 0) {
            return 1.0 / d4;
        }
        if (mode == 1) {
            return s * s * s * s / d4;
        }
        return s * s / d4;
    };

    qwqdsp_filter::Ladder f;
    f.Set(static_cast<float>(w0), static_cast<float>(k), static_cast<float>(q));

    f.Reset();
    CheckResponse("Ladder LP", [&](double x) { return f.Tick(static_cast<float>(x)); },
                  [&](double w) { return ideal(w, 0); });
    f.Reset();
    CheckResponse("Ladder HP", [&](double x) { return f.TickHighpass(static_cast<float>(x)); },
                  [&](double w) { return ideal(w, 1); });
    f.Reset();
    CheckResponse("Ladder BP", [&](double x) { return f.TickBandpass(static_cast<float>(x)); },
                  [&](double w) { return ideal(w, 2); });
    f.Reset();
    CheckResponse("Ladder true HP", [&](double x) { return f.TickTrueHighpass(static_cast<float>(x)).hp4; },
                  [&](double w) {
                      Complex const s = S(w, g);
                      Complex const one_plus_s = 1.0 + s;
                      Complex const s2 = s * s;
                      Complex const s4 = s2 * s2;
                      Complex const d4 = one_plus_s * one_plus_s * one_plus_s * one_plus_s;
                      return s4 / (d4 + k * s4);
                  });

    auto true_bp = [&](double w, int output) {
        Complex const s = S(w, g);
        Complex const d = s * s + r2 * s + 1.0;
        Complex const denom = d * d - k * 0.25 * (1.0 / (q * q)) * s * s;
        if (output == 0) {
            return 0.5 / q * s * d / denom;
        }
        return 0.25 / (q * q) * s * s / denom;
    };
    f.Reset();
    CheckResponse("Ladder true BPn1", [&](double x) { return f.TickTrueBandpass(static_cast<float>(x)).BPn1; },
                  [&](double w) { return true_bp(w, 0); });
    f.Reset();
    CheckResponse("Ladder true BPn2", [&](double x) { return f.TickTrueBandpass(static_cast<float>(x)).BPn2; },
                  [&](double w) { return true_bp(w, 1); });
}

void CheckLadder8Pole() {
    constexpr double w0 = 0.2 * kPi;
    constexpr double k = 0.5;
    double const g = std::tan(w0 / 2.0) * (1.0 + std::numbers::sqrt2_v<double>);

    qwqdsp_filter::Ladder8Pole f;
    f.Set(static_cast<float>(w0), static_cast<float>(k));

    auto ideal = [&](double w, int order) {
        Complex const s = S(w, g);
        Complex const one_plus_s = 1.0 + s;
        Complex denominator = one_plus_s;
        for (int i = 1; i < 8; ++i) {
            denominator *= one_plus_s;
        }
        denominator += k;
        Complex numerator{1.0, 0.0};
        for (int i = 0; i < 8 - order; ++i) {
            numerator *= one_plus_s;
        }
        return numerator / denominator;
    };

    f.Reset();
    CheckResponse("Ladder8Pole lp8", [&](double x) { return f.Tick(static_cast<float>(x)).lp8; },
                  [&](double w) { return ideal(w, 8); });
    f.Reset();
    CheckResponse("Ladder8Pole lp4", [&](double x) { return f.Tick(static_cast<float>(x)).lp4; },
                  [&](double w) { return ideal(w, 4); });
}
} // namespace

int main() {
    CheckOnePole();
    CheckOnePoleShelf();
    CheckOtaOnePole();
    CheckSvf();
    CheckSvfShelf();
    CheckTransposeSallenKey();
    CheckTransposeSvf();
    CheckTransposeSvfMultiInput();
    CheckSvfBandpass();
    CheckLadder();
    CheckLadder8Pole();

    if (g_failures != 0) {
        std::cout << std::format("{} failure(s)\n", g_failures);
        return 1;
    }
    std::cout << "PASS\n";
    return 0;
}
