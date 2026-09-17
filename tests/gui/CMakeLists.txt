# ------------------------------------------------------------
# GUI 测试分发
# ------------------------------------------------------------
if(QWQDSP_USE_RAYLIB)
    # ----- General -----
    add_qwqdsp_gui_test(interpolations interpolations.cpp)
    add_qwqdsp_gui_test(auto_notch_rt auto_notch_rt.cpp)
    add_qwqdsp_gui_test(waveform2 waveform2.cpp)
    add_qwqdsp_gui_test(iir_design iir_design.cpp NO_MINIAUDIO)
    add_qwqdsp_gui_test(filter_designers filter_designers.cpp NO_MINIAUDIO)

    add_subdirectory(audiofx)
    add_subdirectory(spectral)
    add_subdirectory(synth)
    add_subdirectory(voicing)
    add_subdirectory(pitch)
endif()
