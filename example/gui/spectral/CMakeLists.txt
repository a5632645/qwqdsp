# ------------------------------------------------------------
# Spectral (GUI)
# ------------------------------------------------------------
add_qwqdsp_gui_example(pfb_view polyphase_filter_bank_view.cpp)
target_compile_options(qwqdsp-example-pfb_view PUBLIC -mavx2 -mfma)

add_qwqdsp_gui_example(pfb2_view polyphase_filter_bank_view2.cpp)
target_compile_options(qwqdsp-example-pfb2_view PUBLIC -mavx2 -mfma)

add_qwqdsp_gui_example(resonate_bank resonate_bank.cpp)
target_compile_options(qwqdsp-example-resonate_bank PUBLIC -mavx2 -mfma)

add_qwqdsp_gui_example(octave_spectrum octave_spectrum.cpp)
add_qwqdsp_gui_example(multi_resolution_spectrum multi_resolution_spectrum.cpp)
add_qwqdsp_gui_example(portnoff_spectrum portnoff_spectrum.cpp)
add_qwqdsp_gui_example(reassignment reassignment.cpp)
