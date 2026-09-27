# ------------------------------------------------------------
# 无头 example（跑完导出 WAV/PNG 供试听试看）
#
# 新增：add_qwqdsp_example(<name> <file> [EXTRA_LIBS ...]) -> qwqdsp-example-<name>
# ------------------------------------------------------------

# ----- Audio FX -----
add_qwqdsp_example(auto_notch audiofx/auto_notch.cpp)
add_qwqdsp_example(conv_feedback audiofx/conv_feedback.cpp)
add_qwqdsp_example(conv_feedback2 audiofx/conv_feedback2.cpp)
add_qwqdsp_example(freq_shifter audiofx/freq_shifter.cpp)
add_qwqdsp_example(oversample audiofx/oversample.cpp)
add_qwqdsp_example(pitch_shifter audiofx/pitch_shifter.cpp)
add_qwqdsp_example(pitch_shift audiofx/pitch_shift.cpp)
add_qwqdsp_example(pitch_quantizers audiofx/pitch_quantizers.cpp)
add_qwqdsp_example(piwarp audiofx/piwarp.cpp)
add_qwqdsp_example(resample audiofx/resample.cpp)
add_qwqdsp_example(wiener_filter audiofx/wiener_filter.cpp)

# ----- Reverb -----
add_qwqdsp_example(stretch_reverb audiofx/reverb/stretch_reverb.cpp)
add_qwqdsp_example(spectral_reverb audiofx/reverb/spectral_reverb.cpp)
add_qwqdsp_example(reverb audiofx/reverb/reverb.cpp)
add_qwqdsp_example(conv_noise audiofx/reverb/conv_noise.cpp)
add_qwqdsp_example(velvet_noise audiofx/reverb/velvet_noise.cpp)

# ----- Filter -----
add_qwqdsp_example(any_hilbert filter/any_hilbert.cpp)

# ----- Onset -----
add_qwqdsp_example(oneset oneset/oneset.cpp)

# ----- Pitch -----
add_qwqdsp_example(vad pitch/vad.cpp)
add_qwqdsp_example(pitch_osc pitch/pitch_osc.cpp)
add_qwqdsp_example(pyin-vis pitch/pyin_vis.cpp EXTRA_LIBS qwqdsp_support_stb_image)

if(QWQDSP_USE_EIGEN)
    add_qwqdsp_example(swift_f0-vis pitch/swift_f0_vis.cpp EXTRA_LIBS qwqdsp_support_stb_image)
endif()

# ----- Spectral -----
add_qwqdsp_example(dft_filter_bank spectral/dft_filter_bank.cpp)
