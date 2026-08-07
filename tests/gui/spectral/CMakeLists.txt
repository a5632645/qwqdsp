# ------------------------------------------------------------
# Spectral (GUI)
# ------------------------------------------------------------
add_qwqdsp_gui_test(pfb_view polyphase_filter_bank_view.cpp)
target_compile_options(qwqdsp-gui-pfb_view PUBLIC -mavx2 -mfma)

add_qwqdsp_gui_test(pfb2_view polyphase_filter_bank_view2.cpp)
target_compile_options(qwqdsp-gui-pfb2_view PUBLIC -mavx2 -mfma)

add_qwqdsp_gui_test(resonate_bank resonate_bank.cpp)
target_compile_options(qwqdsp-gui-resonate_bank PUBLIC -mavx2 -mfma)

add_qwqdsp_gui_test(spectrum spectrum/spectrum.cpp)
target_sources(qwqdsp-gui-spectrum PUBLIC spectrum/spectrum_driver.cpp)

add_qwqdsp_gui_test(spectrum2 spectrum2/spectrum.cpp)
target_sources(qwqdsp-gui-spectrum2 PUBLIC spectrum2/spectrum_driver.cpp)

add_qwqdsp_gui_test(spectrum3 spectrum3.cpp)
add_qwqdsp_gui_test(spectrum4 spectrum4.cpp)
add_qwqdsp_gui_test(reassignment reassignment.cpp)
