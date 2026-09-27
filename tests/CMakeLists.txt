# ------------------------------------------------------------
# tests — 真正的测试：程序自己判定（断言 / 非零退出 / PASS-FAIL / 与参考逐点比对）
#
# 新增：add_qwqdsp_test(<name> <file> [EXTRA_LIBS ...]) -> qwqdsp-test-<name>
# ------------------------------------------------------------

add_qwqdsp_test(cephes cephes/cephes.cpp)

add_qwqdsp_test(parallel_allpass filter/parallel_allpass.cpp)

add_qwqdsp_test(autocorrelation_pitch_shift_rt audiofx/autocorrelation_delay_pitch_shifter_test.cpp)

add_qwqdsp_test(pyin pitch/pyin.cpp)
add_qwqdsp_test(pmpm pitch/pmpm.cpp)
add_qwqdsp_test(pitch_algorithms pitch/pitch_algorithms.cpp)

if(QWQDSP_USE_EIGEN)
    add_qwqdsp_test(swift_f0 pitch/swift_f0.cpp)
endif()

add_qwqdsp_test(minimax_sine oneset/minimax_sine_test.cpp)
