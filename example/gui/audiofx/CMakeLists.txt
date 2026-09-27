# ------------------------------------------------------------
# Audio FX (GUI)
# ------------------------------------------------------------
add_qwqdsp_gui_example(echo echo.cpp)
add_qwqdsp_gui_example(reverb_rt reverb.cpp)
add_qwqdsp_gui_example(pitch_quantize_rt pitch_quantize_rt.cpp)
add_qwqdsp_gui_example(psola_rt psola_rt.cpp)
add_qwqdsp_gui_example(filters filters.cpp NO_MINIAUDIO)
add_qwqdsp_gui_example(filters2 filters2.cpp NO_MINIAUDIO)
add_qwqdsp_gui_example(formant_shifter formant_shifter.cpp)
add_qwqdsp_gui_example(pitch_shift_rt pitch_shift_rt/pitch_shift_rt.cpp)
