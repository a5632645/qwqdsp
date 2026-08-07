# ------------------------------------------------------------
# Audio FX (GUI)
# ------------------------------------------------------------
add_qwqdsp_gui_test(echo echo.cpp)
add_qwqdsp_gui_test(reverb_rt reverb.cpp)
add_qwqdsp_gui_test(pitch_quantize_rt pitch_quantize_rt.cpp)
add_qwqdsp_gui_test(psola_rt psola_rt.cpp)
add_qwqdsp_gui_test(filters filters.cpp NO_MINIAUDIO)
add_qwqdsp_gui_test(filters2 filters2.cpp NO_MINIAUDIO)
