// Throughput comparison of three LP2BP realizations:
//   1. direct structural LP2BP
//   2. transposed structural LP2BP
//   3. cascade of two transposed TPT SVFs from LP->BP root mapping
#include <benchmark/benchmark.h>

#include <qwqdsp/filter/iir_design.hpp>
#include <qwqdsp/filter/tpt/svf_tpt_bandpass.hpp>
#include <qwqdsp/filter/tpt/transpose_svf_tpt.hpp>
#include <qwqdsp/filter/tpt/transpose_svf_tpt_bandpass.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <utility>

namespace {
using ZPK = qwqdsp_filter::IIRDesign::ZPK;

struct Config {
    ZPK zpk;
    float w1;
    float w2;

    Config() {
        std::array<ZPK, 1> p{};
        if (!qwqdsp_filter::IIRDesign::Elliptic(p, 1, 0.5, 40.0, false) || !p[0].z) {
            std::abort();
        }
        zpk = p[0];
        constexpr float pi = 3.14159265358979323846f;
        w1 = 0.2f * pi;
        w2 = 0.2001f * pi;
    }
};

Config const& Cfg() {
    static Config const c;
    return c;
}

struct Direct {
    qwqdsp_filter::SvfTPTBandpass f;

    Direct() : f(Cfg().zpk.p, *Cfg().zpk.z, Cfg().zpk.k) { f.Set(Cfg().w1, Cfg().w2); }

    float Tick(float x) noexcept { return f.Tick(x); }
};

struct Transposed {
    qwqdsp_filter::TransposeSvfTPTBandpass f;

    Transposed() : f(Cfg().zpk.p, *Cfg().zpk.z, Cfg().zpk.k) { f.Set(Cfg().w1, Cfg().w2); }

    float Tick(float x) noexcept { return f.Tick(x); }
};

struct Cascade {
    qwqdsp_filter::TransposeSvfTPT f1;
    qwqdsp_filter::TransposeSvfTPT f2;

    Cascade() {
        auto const& c = Cfg();
        double const aw1 = std::tan(c.w1 / 2.0);
        double const aw2 = std::tan(c.w2 / 2.0);
        double const bw = aw2 - aw1;
        double const w1w2 = aw1 * aw2;

        auto split = [&](std::complex<double> root) {
            std::complex<double> const delta = std::sqrt(root * root * bw * bw - 4.0 * w1w2);
            return std::pair{(root * bw + delta) * 0.5, (root * bw - delta) * 0.5};
        };
        auto const [p1, p2] = split(c.zpk.p);
        auto const [z1, z2] = split(*c.zpk.z);
        double const ks = std::sqrt(c.zpk.k);

        auto set = [&](qwqdsp_filter::TransposeSvfTPT& f, std::complex<double> p, std::complex<double> z) {
            double const g = std::abs(p);
            double const r2 = -2.0 * p.real() / g;
            f.SetAnalogCoeff(static_cast<float>(g), static_cast<float>(r2));
            f.SetInputMix(static_cast<float>(ks), 0.0f,
                          static_cast<float>(ks * std::norm(z) / (g * g)));
        };
        set(f1, p1, z1);
        set(f2, p2, z2);
    }

    float Tick(float x) noexcept { return f2.Tick(f1.Tick(x)); }
};

struct Tdf2Cascade {
    std::array<qwqdsp_filter::BiquadCoeff, 2> bq{};
    std::array<float, 4> st{};

    Tdf2Cascade() {
        auto const& c = Cfg();
        double const aw1 = std::tan(c.w1 / 2.0);
        double const aw2 = std::tan(c.w2 / 2.0);
        std::array<ZPK, 2> bp{};
        bp[0] = c.zpk;
        qwqdsp_filter::IIRDesign::ProtyleToBandpass2(bp, 1, aw1, aw2);
        qwqdsp_filter::IIRDesign::Bilinear(bp, 0.5);
        qwqdsp_filter::IIRDesign::NormalizeSectionGains(bp);
        qwqdsp_filter::IIRDesign::TfToBiquad(bp, bq);
    }

    float Tick(float x) noexcept {
        float const y1 = bq[0].b0 * x + st[0];
        st[0] = bq[0].b1 * x - bq[0].a1 * y1 + st[1];
        st[1] = bq[0].b2 * x - bq[0].a2 * y1;

        float const y2 = bq[1].b0 * y1 + st[2];
        st[2] = bq[1].b1 * y1 - bq[1].a1 * y2 + st[3];
        st[3] = bq[1].b2 * y1 - bq[1].a2 * y2;
        return y2;
    }
};

constexpr int kBlockSize = 512;

template <typename Filter>
void BM_Filter(benchmark::State& state) {
    Filter f;
    std::array<float, kBlockSize> input{};
    for (int i = 0; i < kBlockSize; ++i) {
        input[i] = std::sin(0.013f * static_cast<float>(i));
    }

    float acc = 0.0f;
    for (auto _ : state) {
        for (float x : input) {
            acc += f.Tick(x);
        }
        benchmark::DoNotOptimize(acc);
    }
    state.SetItemsProcessed(state.iterations() * kBlockSize);
}
} // namespace

BENCHMARK_TEMPLATE(BM_Filter, Direct)->Name("direct_LP2BP");
BENCHMARK_TEMPLATE(BM_Filter, Transposed)->Name("transposed_LP2BP");
BENCHMARK_TEMPLATE(BM_Filter, Cascade)->Name("cascade_transposed_SVF");
BENCHMARK_TEMPLATE(BM_Filter, Tdf2Cascade)->Name("tdf2_cascade");

BENCHMARK_MAIN();
