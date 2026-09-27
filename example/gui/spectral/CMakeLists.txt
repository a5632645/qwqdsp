# ------------------------------------------------------------
# Spectral (GUI)
# ------------------------------------------------------------
add_qwqdsp_gui_example(pfb_view polyphase_filter_bank_view.cpp)
target_compile_options(qwqdsp-example-pfb_view PUBLIC -mavx2 -mfma)

add_qwqdsp_gui_example(pfb2_view polyphase_filter_bank_view2.cpp)
target_compile_options(qwqdsp-example-pfb2_view PUBLIC -mavx2 -mfma)

add_qwqdsp_gui_example(resonate_bank resonate_bank.cpp)
target_compile_options(qwqdsp-example-resonate_bank PUBLIC -mavx2 -mfma)

add_qwqdsp_gui_example(spectrum spectrum/spectrum.cpp)
target_sources(qwqdsp-example-spectrum PUBLIC spectrum/spectrum_driver.cpp)

add_qwqdsp_gui_example(spectrum2 spectrum2/spectrum.cpp)
target_sources(qwqdsp-example-spectrum2 PUBLIC spectrum2/spectrum_driver.cpp)

add_qwqdsp_gui_example(spectrum3 spectrum3.cpp)
add_qwqdsp_gui_example(spectrum4 spectrum4.cpp)
add_qwqdsp_gui_example(reassignment reassignment.cpp)
