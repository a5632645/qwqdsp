# ------------------------------------------------------------
# labs — 探索：一次性数值探针 / 算法对照 / 论文图复现 / 脚本实验
#
# 新增：add_qwqdsp_lab(<name> <file> [EXTRA_LIBS ...]) -> qwqdsp-lab-<name>
# ------------------------------------------------------------

# ----- 一次性数值探针（只算/只打印，不产出可听可看的结果）-----
add_qwqdsp_lab(biquad filter/biquad.cpp)
add_qwqdsp_lab(fir_design filter/fir_design.cpp)
add_qwqdsp_lab(minimum_phase_fir filter/minimum_phase_fir.cpp)
add_qwqdsp_lab(residual_chebyshev filter/residual_chebyshev.cpp)
add_qwqdsp_lab(reverse_iir filter/reverse_iir.cpp)

add_qwqdsp_lab(thiran_allpass thiran_filter/allpass.cpp)
add_qwqdsp_lab(thiran_allpole thiran_filter/allpole.cpp)
add_qwqdsp_lab(thiran_zero_pole thiran_filter/zero_pole.cpp)

add_qwqdsp_lab(fft_interpolation spectral/fft_interpolation.cpp)
add_qwqdsp_lab(ifft_random_phase spectral/ifft_random_phase.cpp)
add_qwqdsp_lab(fft_fold spectral/fold.cpp)

add_qwqdsp_lab(dsf synth/dsf.cpp)
add_qwqdsp_lab(stupid_resynthsis synth/stupid_resynthsis.cpp)

add_qwqdsp_lab(convolution audiofx/convolution.cpp)
add_qwqdsp_lab(pitch_shift_transient pitch_shift_transient/transient_probe.cpp)

# ----- 参数扫描 / 算法对照 -----
add_qwqdsp_lab(lp2bp_noise filter/lp2bp/noise.cpp)
add_qwqdsp_lab(piwarp5_comb audiofx/piwarp5_comb.cpp)
add_qwqdsp_lab(pmpm-vs-pyin pitch/pmpm_vs_pyin.cpp)
add_qwqdsp_lab(autotune autotune/autotune.cpp)

# ----- raylib 实验沙盒 -----
if(QWQDSP_USE_RAYLIB)
    add_subdirectory(playing)
endif()
