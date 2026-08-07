# ------------------------------------------------------------
# Pitch
# ------------------------------------------------------------
add_qwqdsp_dsp_test(pyin pyin.cpp)
add_qwqdsp_dsp_test(pmpm pmpm.cpp)
add_qwqdsp_dsp_test(pmpm-vs-pyin pmpm_vs_pyin.cpp)
add_qwqdsp_dsp_test(pitch_algorithms pitch_algorithms.cpp)
add_qwqdsp_dsp_test(pyin-vis pyin_vis.cpp EXTRA_LIBS qwqdsp_support_stb_image)
add_qwqdsp_dsp_test(vad vad.cpp)
add_qwqdsp_dsp_test(pitch_osc pitch_osc.cpp)

if(QWQDSP_USE_EIGEN)
    add_qwqdsp_dsp_test(swift_f0 swift_f0.cpp)
    add_qwqdsp_dsp_test(swift_f0-vis swift_f0_vis.cpp EXTRA_LIBS qwqdsp_support_stb_image)
endif()
