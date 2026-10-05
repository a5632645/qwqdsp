# ------------------------------------------------------------
# tests — 真正的测试：程序自己判定（断言 / 非零退出 / PASS-FAIL / 与参考逐点比对）
#
# 新增：add_qwqdsp_test(<name> <file> [EXTRA_LIBS ...]) -> qwqdsp-test-<name>
# ------------------------------------------------------------

add_qwqdsp_test(cephes cephes/cephes.cpp)

add_qwqdsp_test(parallel_allpass filter/parallel_allpass.cpp)
add_qwqdsp_test(iir_design filter/iir_design.cpp)
add_qwqdsp_test(halfband_allpass_ovs filter/halfband_allpass_ovs.cpp)

add_qwqdsp_test(autocorrelation_pitch_shift_rt audiofx/autocorrelation_delay_pitch_shifter_test.cpp)

add_qwqdsp_test(pyin pitch/pyin.cpp)
add_qwqdsp_test(pmpm pitch/pmpm.cpp)
add_qwqdsp_test(pitch_algorithms pitch/pitch_algorithms.cpp)

if(QWQDSP_USE_EIGEN)
    add_qwqdsp_test(swift_f0 pitch/swift_f0.cpp)
endif()

add_qwqdsp_test(minimax_sine oneset/minimax_sine_test.cpp)

# ------------------------------------------------------------
# 慢测试打 slow 标签
#
# 收尾默认只跑快集：`ctest -LE slow`（实测 7 例 ≈ 6 s），跑全量留给交付/大改前。
# 口径：单例 Debug 下 > 5 s 就登记进来。目前 swift_f0 ≈ 33–37 s，占全部测试时间的 ~85%
#（Debug 下跑神经网络推理，与输入长度无关——它用的是 1 秒合成正弦）。
# 改了对应子系统时仍要单独跑：`ctest -R swift_f0`。
# ------------------------------------------------------------
set(QWQDSP_SLOW_TESTS swift_f0)
foreach(_slow IN LISTS QWQDSP_SLOW_TESTS)
    if(TARGET qwqdsp-test-${_slow})
        set_tests_properties(qwqdsp-test-${_slow} PROPERTIES LABELS slow)
    endif()
endforeach()
