// Compare fixed-point state roundoff noise and coefficient sensitivity of:
//   1. direct structural LP2BP (two bandpass integrators, nested ZDF)
//   2. transposed structural LP2BP
//   3. cascade of two transposed TPT SVFs obtained by LP->BP root mapping
#include <qwqdsp/filter/iir_design.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

namespace {
struct Section {
    double g;
    double g2;
    double r2;
    double gr2g2;
    double d;
    double ff;
    double g2_lp;
};

struct BiquadD {
    double b0;
    double b1;
    double b2;
    double a1;
    double a2;
};

struct Coeffs {
    double g;
    double s;
    double r2;
    double r2s;
    double a;
    double b;
    double d;
    double inv_d_direct;
    double at;
    double bt;
    double dt;
    double inv_d_trans;
    double hp;
    double bp;
    double lp;
    Section c1;
    Section c2;
    BiquadD bq_raw1;
    BiquadD bq_raw2;
    BiquadD bq_norm1;
    BiquadD bq_norm2;
};

void QuantizeStates(std::array<double, 4>& st, std::array<double, 4> const& step) {
    for (size_t i = 0; i < 4; ++i) {
        st[i] = std::round(st[i] / step[i]) * step[i];
    }
}

struct Direct {
    Coeffs c;
    std::array<double, 4> st{};

    void Reset() noexcept { st = {}; }

    double Tick(double x, bool quant, std::array<double, 4> const& step) noexcept {
        double const u = c.g * x - c.g * st[1] + st[0];
        double const v = st[2] - c.g * st[3];
        double const b1 = (c.d * u - c.b * v) * c.inv_d_direct;
        double const b2 = (c.a * v + c.b * u) * c.inv_d_direct;
        double const l1 = c.g * b1 + st[1];
        double const l2 = c.g * b2 + st[3];

        st[0] = 2 * b1 - st[0];
        st[1] = 2 * l1 - st[1];
        st[2] = 2 * b2 - st[2];
        st[3] = 2 * l2 - st[3];
        if (quant) {
            QuantizeStates(st, step);
        }

        double const hp = x - c.r2s * b1 - c.s * b2;
        double const bp = c.s * b1;
        double const lp = c.s * b2;
        return c.hp * hp + c.bp * bp + c.lp * lp;
    }
};

struct Transposed {
    Coeffs c;
    std::array<double, 4> st{};

    void Reset() noexcept { st = {}; }

    double Tick(double x, bool quant, std::array<double, 4> const& step) noexcept {
        double const x_hp = c.hp * x;
        double const x_bp = c.bp * x;
        double const x_lp = c.lp * x;

        double const u = c.g * (x_lp - x_hp) - c.g * st[1] + st[0];
        double const v = c.g * (x_bp - c.r2 * x_hp) - c.g * st[3] + st[2];
        double const b1 = (c.dt * u - c.bt * v) * c.inv_d_trans;
        double const b2 = (c.at * v + c.bt * u) * c.inv_d_trans;
        double const l1 = c.g * b1 + st[1];
        double const l2 = c.g * b2 + st[3];

        st[0] = 2 * b1 - st[0];
        st[1] = 2 * l1 - st[1];
        st[2] = 2 * b2 - st[2];
        st[3] = 2 * l2 - st[3];
        if (quant) {
            QuantizeStates(st, step);
        }

        return x_hp + c.s * b2;
    }
};

struct Cascade {
    Coeffs c;
    std::array<double, 4> st{}; // s1_1, s2_1, s1_2, s2_2

    void Reset() noexcept { st = {}; }

    double Tick(double x, bool quant, std::array<double, 4> const& step) noexcept {
        double const y1 = (c.c1.ff * x + st[0] + st[1]) * c.c1.d;
        double const t1 = st[0];
        st[0] = t1 + 2 * (c.c1.g2_lp * x - c.c1.g2 * y1);
        st[1] += 2 * (c.c1.g2_lp * x - c.c1.gr2g2 * y1 + t1);
        if (quant) {
            st[0] = std::round(st[0] / step[0]) * step[0];
            st[1] = std::round(st[1] / step[1]) * step[1];
        }

        double const y2 = (c.c2.ff * y1 + st[2] + st[3]) * c.c2.d;
        double const t2 = st[2];
        st[2] = t2 + 2 * (c.c2.g2_lp * y1 - c.c2.g2 * y2);
        st[3] += 2 * (c.c2.g2_lp * y1 - c.c2.gr2g2 * y2 + t2);
        if (quant) {
            st[2] = std::round(st[2] / step[2]) * step[2];
            st[3] = std::round(st[3] / step[3]) * step[3];
        }

        return y2;
    }
};

template <bool kNormalized>
struct Tdf2 {
    Coeffs c;
    std::array<double, 4> st{}; // s1_1, s2_1, s1_2, s2_2

    void Reset() noexcept { st = {}; }

    double Tick(double x, bool quant, std::array<double, 4> const& step) noexcept {
        BiquadD const& b1 = kNormalized ? c.bq_norm1 : c.bq_raw1;
        BiquadD const& b2 = kNormalized ? c.bq_norm2 : c.bq_raw2;

        double const y1 = b1.b0 * x + st[0];
        st[0] = b1.b1 * x - b1.a1 * y1 + st[1];
        st[1] = b1.b2 * x - b1.a2 * y1;
        if (quant) {
            st[0] = std::round(st[0] / step[0]) * step[0];
            st[1] = std::round(st[1] / step[1]) * step[1];
        }

        double const y2 = b2.b0 * y1 + st[2];
        st[2] = b2.b1 * y1 - b2.a1 * y2 + st[3];
        st[3] = b2.b2 * y1 - b2.a2 * y2;
        if (quant) {
            st[2] = std::round(st[2] / step[2]) * step[2];
            st[3] = std::round(st[3] / step[3]) * step[3];
        }
        return y2;
    }
};

using Tdf2Raw = Tdf2<false>;
using Tdf2Norm = Tdf2<true>;

Coeffs MakeCoeffs(qwqdsp_filter::IIRDesign::ZPK const& zpk, double w1, double w2) {
    double const wc = std::abs(zpk.p);
    double const r2 = -2.0 * zpk.p.real() / wc;
    double const wz = std::abs(*zpk.z);

    double const aw1 = std::tan(w1 / 2.0);
    double const aw2 = std::tan(w2 / 2.0);
    double const bw = aw2 - aw1;
    double const w0 = std::sqrt(aw1 * aw2);

    Coeffs c{};
    c.g = w0;
    c.s = wc * bw / w0;
    c.r2 = r2;
    c.r2s = r2 * c.s;
    c.a = 1 + c.r2s * c.g + c.g * c.g;
    c.b = c.s * c.g;
    c.d = 1 + c.g * c.g;
    c.inv_d_direct = 1 / (c.a * c.d + c.b * c.b);
    c.at = 1 + c.g * c.g;
    c.bt = c.s * c.g;
    c.dt = 1 + c.g * c.g + c.r2s * c.g;
    c.inv_d_trans = 1 / (c.at * c.dt + c.bt * c.bt);
    c.hp = zpk.k;
    c.bp = 0;
    c.lp = zpk.k * wz * wz / (wc * wc);

    auto split = [&](std::complex<double> root) {
        std::complex<double> const delta = std::sqrt(root * root * bw * bw - 4.0 * aw1 * aw2);
        return std::pair{(root * bw + delta) * 0.5, (root * bw - delta) * 0.5};
    };
    auto const [p1, p2] = split(zpk.p);
    auto const [z1, z2] = split(*zpk.z);

    double const ks = std::sqrt(zpk.k);
    auto make_section = [&](std::complex<double> p, std::complex<double> z) {
        double const g = std::abs(p);
        double const r = -2.0 * p.real() / g;
        double const g2 = g * g;
        double const lp_mix = ks * std::norm(z) / g2;
        Section s{};
        s.g = g;
        s.g2 = g2;
        s.r2 = r;
        s.gr2g2 = g * r + g2;
        s.d = 1 / (1 + r * g + g2);
        s.g2_lp = g2 * lp_mix;
        s.ff = ks + s.g2_lp;
        return s;
    };
    c.c1 = make_section(p1, z1);
    c.c2 = make_section(p2, z2);

    std::array<qwqdsp_filter::IIRDesign::ZPK, 2> bp{};
    bp[0] = zpk;
    qwqdsp_filter::IIRDesign::ProtyleToBandpass2(bp, 1, aw1, aw2);
    qwqdsp_filter::IIRDesign::Bilinear(bp, 0.5);
    auto bp_raw = bp;
    qwqdsp_filter::IIRDesign::NormalizeSectionGains(bp);

    std::array<qwqdsp_filter::BiquadCoeff, 2> bq_raw{};
    std::array<qwqdsp_filter::BiquadCoeff, 2> bq_norm{};
    qwqdsp_filter::IIRDesign::TfToBiquad(bp_raw, bq_raw);
    qwqdsp_filter::IIRDesign::TfToBiquad(bp, bq_norm);
    c.bq_raw1 = {bq_raw[0].b0, bq_raw[0].b1, bq_raw[0].b2, bq_raw[0].a1, bq_raw[0].a2};
    c.bq_raw2 = {bq_raw[1].b0, bq_raw[1].b1, bq_raw[1].b2, bq_raw[1].a1, bq_raw[1].a2};
    c.bq_norm1 = {bq_norm[0].b0, bq_norm[0].b1, bq_norm[0].b2, bq_norm[0].a1, bq_norm[0].a2};
    c.bq_norm2 = {bq_norm[1].b0, bq_norm[1].b1, bq_norm[1].b2, bq_norm[1].a1, bq_norm[1].a2};
    return c;
}

struct Stats {
    double err_rms = 0;
    double err_max = 0;
    double snr_db = 0;
};

template <typename F>
std::vector<double> Reference(F& f, Coeffs const& c, std::vector<double> const& in, size_t n,
                              std::array<double, 4>& peak) {
    f.c = c;
    f.Reset();
    peak = {0, 0, 0, 0};
    std::vector<double> y(n);
    for (size_t i = 0; i < n; ++i) {
        y[i] = f.Tick(in[i], false, {});
        for (size_t j = 0; j < 4; ++j) {
            peak[j] = std::max(peak[j], std::abs(f.st[j]));
        }
    }
    return y;
}

template <typename F>
Stats RunQuantized(F& f, Coeffs const& c, std::vector<double> const& in, size_t n,
                   std::vector<double> const& ref, std::array<double, 4> const& peak, int bits) {
    f.c = c;
    f.Reset();
    std::array<double, 4> step{};
    for (size_t i = 0; i < 4; ++i) {
        step[i] = std::max(peak[i] * 1.1, 1e-12) / static_cast<double>(1ull << bits);
    }

    size_t const warmup = n / 5;
    double sum_ref2 = 0;
    double sum_err2 = 0;
    double max_err = 0;
    for (size_t i = 0; i < n; ++i) {
        double const y = f.Tick(in[i], true, step);
        if (i >= warmup) {
            double const e = y - ref[i];
            sum_ref2 += ref[i] * ref[i];
            sum_err2 += e * e;
            max_err = std::max(max_err, std::abs(e));
        }
    }

    Stats st;
    size_t const count = n - warmup;
    double const ref_rms = std::sqrt(sum_ref2 / static_cast<double>(count));
    st.err_rms = std::sqrt(sum_err2 / static_cast<double>(count));
    st.err_max = max_err;
    st.snr_db = 20.0 * std::log10(ref_rms / std::max(st.err_rms, 1e-30));
    return st;
}

template <typename F>
Stats RunCoeffQuantized(F& f, Coeffs const& c, std::vector<double> const& in, size_t n,
                        std::vector<double> const& ref) {
    f.c = c;
    f.Reset();
    size_t const warmup = n / 5;
    double sum_ref2 = 0;
    double sum_err2 = 0;
    double max_err = 0;
    for (size_t i = 0; i < n; ++i) {
        double const y = f.Tick(in[i], false, {});
        if (i >= warmup) {
            double const e = y - ref[i];
            sum_ref2 += ref[i] * ref[i];
            sum_err2 += e * e;
            max_err = std::max(max_err, std::abs(e));
        }
    }
    Stats st;
    size_t const count = n - warmup;
    double const ref_rms = std::sqrt(sum_ref2 / static_cast<double>(count));
    st.err_rms = std::sqrt(sum_err2 / static_cast<double>(count));
    st.err_max = max_err;
    st.snr_db = 20.0 * std::log10(ref_rms / std::max(st.err_rms, 1e-30));
    return st;
}

Coeffs QuantizeCoeffs(Coeffs c, int bits) {
    double const step = std::ldexp(1.0, -(bits - 1));
    auto q = [step](double v) { return std::round(v / step) * step; };
    c.g = q(c.g);
    c.s = q(c.s);
    c.r2 = q(c.r2);
    c.r2s = q(c.r2s);
    c.a = q(c.a);
    c.b = q(c.b);
    c.d = q(c.d);
    c.inv_d_direct = q(c.inv_d_direct);
    c.at = q(c.at);
    c.bt = q(c.bt);
    c.dt = q(c.dt);
    c.inv_d_trans = q(c.inv_d_trans);
    c.hp = q(c.hp);
    c.bp = q(c.bp);
    c.lp = q(c.lp);
    auto qs = [&](Section s) {
        s.g = q(s.g);
        s.g2 = q(s.g2);
        s.r2 = q(s.r2);
        s.gr2g2 = q(s.gr2g2);
        s.d = q(s.d);
        s.ff = q(s.ff);
        s.g2_lp = q(s.g2_lp);
        return s;
    };
    c.c1 = qs(c.c1);
    c.c2 = qs(c.c2);
    auto qb = [&](BiquadD b) {
        b.b0 = q(b.b0);
        b.b1 = q(b.b1);
        b.b2 = q(b.b2);
        b.a1 = q(b.a1);
        b.a2 = q(b.a2);
        return b;
    };
    c.bq_raw1 = qb(c.bq_raw1);
    c.bq_raw2 = qb(c.bq_raw2);
    c.bq_norm1 = qb(c.bq_norm1);
    c.bq_norm2 = qb(c.bq_norm2);
    return c;
}

void RunCase(char const* name, qwqdsp_filter::IIRDesign::ZPK const& zpk, double w1, double w2, size_t n,
             std::vector<double> const& in) {
    Coeffs const c = MakeCoeffs(zpk, w1, w2);

    Direct direct;
    Transposed transposed;
    Cascade cascade;
    Tdf2Raw tdf2_raw;
    Tdf2Norm tdf2_norm;
    std::array<double, 4> pd{};
    std::array<double, 4> pt{};
    std::array<double, 4> pc{};
    std::array<double, 4> pr{};
    std::array<double, 4> pn{};
    auto const rd = Reference(direct, c, in, n, pd);
    auto const rt = Reference(transposed, c, in, n, pt);
    auto const rc = Reference(cascade, c, in, n, pc);
    auto const rr = Reference(tdf2_raw, c, in, n, pr);
    auto const rn = Reference(tdf2_norm, c, in, n, pn);

    double diff_dt = 0;
    double diff_dc = 0;
    double diff_dr = 0;
    double diff_dn = 0;
    for (size_t i = 0; i < n; ++i) {
        diff_dt = std::max(diff_dt, std::abs(rd[i] - rt[i]));
        diff_dc = std::max(diff_dc, std::abs(rd[i] - rc[i]));
        diff_dr = std::max(diff_dr, std::abs(rd[i] - rr[i]));
        diff_dn = std::max(diff_dn, std::abs(rd[i] - rn[i]));
    }

    std::printf("\n%s  (N=%zu, Q~%.1f)\n", name, n,
                std::sqrt(std::tan(w1 / 2) * std::tan(w2 / 2)) / (std::tan(w2 / 2) - std::tan(w1 / 2)));
    std::printf("  reference max diff: trans %.3e, cascade %.3e, tdf2 raw %.3e, tdf2 norm %.3e\n",
                diff_dt, diff_dc, diff_dr, diff_dn);
    std::printf("  state peaks direct: %.3e %.3e %.3e %.3e\n", pd[0], pd[1], pd[2], pd[3]);
    std::printf("  state peaks trans : %.3e %.3e %.3e %.3e\n", pt[0], pt[1], pt[2], pt[3]);
    std::printf("  state peaks casc  : %.3e %.3e %.3e %.3e\n", pc[0], pc[1], pc[2], pc[3]);
    std::printf("  state peaks tdf2r : %.3e %.3e %.3e %.3e\n", pr[0], pr[1], pr[2], pr[3]);
    std::printf("  state peaks tdf2n : %.3e %.3e %.3e %.3e\n", pn[0], pn[1], pn[2], pn[3]);

    std::printf("  state bits | direct | trans | cascade | tdf2 raw | tdf2 norm\n");
    for (int bits : {28, 24, 20, 18, 16}) {
        Direct dq;
        Transposed tq;
        Cascade cq;
        Tdf2Raw tr;
        Tdf2Norm tn;
        auto const ds = RunQuantized(dq, c, in, n, rd, pd, bits);
        auto const ts = RunQuantized(tq, c, in, n, rt, pt, bits);
        auto const cs = RunQuantized(cq, c, in, n, rc, pc, bits);
        auto const rs = RunQuantized(tr, c, in, n, rr, pr, bits);
        auto const ns = RunQuantized(tn, c, in, n, rn, pn, bits);
        std::printf("  %10d | %9.2f | %8.2f | %11.2f | %8.2f | %9.2f\n", bits, ds.snr_db, ts.snr_db,
                    cs.snr_db, rs.snr_db, ns.snr_db);
    }

    std::printf("  coeff bits | direct | trans | cascade | tdf2 raw | tdf2 norm\n");
    for (int bits : {24, 20, 16, 14, 12}) {
        Coeffs const cq = QuantizeCoeffs(c, bits);
        Direct dq;
        Transposed tq;
        Cascade cc;
        Tdf2Raw tr;
        Tdf2Norm tn;
        auto const ds = RunCoeffQuantized(dq, cq, in, n, rd);
        auto const ts = RunCoeffQuantized(tq, cq, in, n, rt);
        auto const cs = RunCoeffQuantized(cc, cq, in, n, rc);
        auto const rs = RunCoeffQuantized(tr, cq, in, n, rr);
        auto const ns = RunCoeffQuantized(tn, cq, in, n, rn);
        std::printf("  %10d | %9.2f | %8.2f | %11.2f | %8.2f | %9.2f\n", bits, ds.snr_db, ts.snr_db,
                    cs.snr_db, rs.snr_db, ns.snr_db);
    }
}
} // namespace

int main() {
    std::array<qwqdsp_filter::IIRDesign::ZPK, 1> proto{};
    bool const ok = qwqdsp_filter::IIRDesign::Elliptic(proto, 1, 0.5, 40.0, false);
    if (!ok || !proto[0].z) {
        std::printf("prototype design failed\n");
        return 1;
    }

    constexpr size_t kMaxN = 2'000'000;
    std::vector<double> input(kMaxN);
    std::mt19937 rng(1234567);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (auto& x : input) {
        x = dist(rng);
    }

    constexpr double pi = 3.14159265358979323846;
    RunCase("0.200-0.201 pi", proto[0], 0.200 * pi, 0.201 * pi, 1'000'000, input);
    RunCase("0.2000-0.2001 pi", proto[0], 0.2000 * pi, 0.2001 * pi, 1'000'000, input);
    RunCase("0.20000-0.20001 pi", proto[0], 0.20000 * pi, 0.20001 * pi, kMaxN, input);
    RunCase("0.1000-0.1001 pi", proto[0], 0.1000 * pi, 0.1001 * pi, 1'000'000, input);

    std::array<qwqdsp_filter::IIRDesign::ZPK, 1> proto2{};
    if (qwqdsp_filter::IIRDesign::Elliptic(proto2, 1, 0.1, 60.0, false) && proto2[0].z) {
        RunCase("0.20000-0.20001 pi, 0.1/60 dB", proto2[0], 0.20000 * pi, 0.20001 * pi, kMaxN, input);
    }
    return 0;
}
