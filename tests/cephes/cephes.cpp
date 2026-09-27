#include <qwqdsp/cephes/bessel.hpp>
#include <qwqdsp/cephes/besself.hpp>
#include <qwqdsp/cephes/elliptic.hpp>
#include <qwqdsp/cephes/ellipticf.hpp>

#include <cstdio>
#include <string>

#include "cephes_ref.inc"

static constexpr int kNum = kRefNum;

int main() {
    std::printf("Cephes 特殊函数精度验证 (scipy 参考)\n");
    std::printf("\n");
    std::printf("  %-14s  %-8s  %-16s  %-16s  %s\n",
                "名称", "n", "max_abs_err", "max_rel_err", "errs");
    std::printf("  %s\n",
               "----------------------------------------------------------------------");

    int total_errs = 0;

    // ── float ──
    float buf_f[kNum];

    auto check_f = [&](const char* name, auto fn, const float* ref,
                       float abs_tol = 1e-5f, float rel_tol = 1e-5f) {
        for (int i = 0; i < kNum; ++i)
            buf_f[i] = fn(i);
        total_errs += CheckRef(name, buf_f, ref, kNum, abs_tol, rel_tol);
    };

    auto check_f2 = [&](const char* name, const float* computed, const float* ref) {
        total_errs += CheckRef(name, computed, ref, kNum, 1e-5f, 1e-5f);
    };

    // Besself
    check_f("bf_i0",
        [](int i) { return qwqdsp_cephes::Besself::i0(10.0f * i / kNum); },
        kRef_bf_i0);
    check_f("bf_i1",
        [](int i) { return qwqdsp_cephes::Besself::i1(10.0f * i / kNum); },
        kRef_bf_i1);

    // Ellipticf
    check_f("ef_ellpk",
        [](int i) { return qwqdsp_cephes::Ellipticf::ellpk(1.0f * (i + 1) / kNum); },
        kRef_ef_ellpk);
    check_f("ef_ellpe",
        [](int i) { return qwqdsp_cephes::Ellipticf::ellpe(1.0f * i / kNum); },
        kRef_ef_ellpe);
    check_f("ef_ellik",
        [](int i) { return qwqdsp_cephes::Ellipticf::ellik(2.0f * i / kNum, 0.5f); },
        kRef_ef_ellik);
    check_f("ef_ellie",
        [](int i) { return qwqdsp_cephes::Ellipticf::ellie(2.0f * i / kNum, 0.5f); },
        kRef_ef_ellie);

    {
        float sn[kNum], cn[kNum], dn[kNum], ph[kNum];
        for (int i = 0; i < kNum; ++i)
            qwqdsp_cephes::Ellipticf::ellpj(10.0f * i / kNum, 0.5f,
                                             &sn[i], &cn[i], &dn[i], &ph[i]);
        check_f2("ef_sn", sn, kRef_ef_sn);
        check_f2("ef_cn", cn, kRef_ef_cn);
        total_errs += CheckRef("ef_dn", dn, kRef_ef_dn, kNum, 1e-4f, 1e-4f);
        check_f2("ef_ph", ph, kRef_ef_ph);
    }

    // ── double ──
    double buf_d[kNum];

    auto check_d = [&](const char* name, auto fn, const double* ref,
                       double abs_tol = 1e-12, double rel_tol = 1e-12) {
        for (int i = 0; i < kNum; ++i)
            buf_d[i] = fn(i);
        total_errs += CheckRef(name, buf_d, ref, kNum, abs_tol, rel_tol);
    };

    auto check_d2 = [&](const char* name, const double* computed, const double* ref) {
        total_errs += CheckRef(name, computed, ref, kNum, 1e-12, 1e-12);
    };

    // Bessel
    check_d("b_i0",
        [](int i) { return qwqdsp_cephes::Bessel::i0(30.0 * i / kNum); },
        kRef_b_i0);
    check_d("b_i1",
        [](int i) { return qwqdsp_cephes::Bessel::i1(30.0 * i / kNum); },
        kRef_b_i1);

    // Elliptic
    check_d("e_ellpk",
        [](int i) { return qwqdsp_cephes::Elliptic::ellpk(1.0 * (i + 1) / kNum); },
        kRef_e_ellpk);
    check_d("e_ellpe",
        [](int i) { return qwqdsp_cephes::Elliptic::ellpe(1.0 * i / kNum); },
        kRef_e_ellpe);
    check_d("e_ellik",
        [](int i) { return qwqdsp_cephes::Elliptic::ellik(2.0 * i / kNum, 0.5); },
        kRef_e_ellik);
    check_d("e_ellie",
        [](int i) { return qwqdsp_cephes::Elliptic::ellie(2.0 * i / kNum, 0.5); },
        kRef_e_ellie);

    {
        double sn[kNum], cn[kNum], dn[kNum], ph[kNum];
        for (int i = 0; i < kNum; ++i)
            qwqdsp_cephes::Elliptic::ellpj(10.0 * i / kNum, 0.5,
                                            &sn[i], &cn[i], &dn[i], &ph[i]);
        check_d2("e_sn", sn, kRef_e_sn);
        check_d2("e_cn", cn, kRef_e_cn);
        check_d2("e_dn", dn, kRef_e_dn);
        check_d2("e_ph", ph, kRef_e_ph);
    }

    std::printf("\n");
    if (total_errs == 0)
        std::printf("全部通过 ✓\n");
    else
        std::printf("发现 %d 个超差样本 ✗\n", total_errs);
    return total_errs;
}
