# ------------------------------------------------------------
# Audio FX
# ------------------------------------------------------------
add_qwqdsp_dsp_test(convolution convolution.cpp)
add_qwqdsp_dsp_test(conv_feedback conv_feedback.cpp)
add_qwqdsp_dsp_test(conv_feedback2 conv_feedback2.cpp)
add_qwqdsp_dsp_test(freq_shifter freq_shifter.cpp)
add_qwqdsp_dsp_test(resample resample.cpp)
add_qwqdsp_dsp_test(oversample oversample.cpp)
add_qwqdsp_dsp_test(pitch_shifter pitch_shifter.cpp)
add_qwqdsp_dsp_test(pitch_quantizers pitch_quantizers.cpp)
add_qwqdsp_dsp_test(auto_notch auto_notch.cpp)
add_qwqdsp_dsp_test(wiener_filter wiener_filter.cpp)
add_qwqdsp_dsp_test(piwarp piwarp.cpp)
add_qwqdsp_dsp_test(piwarp5_comb piwarp5_comb.cpp)

add_subdirectory(reverb)
add_subdirectory(pitch_shift)
