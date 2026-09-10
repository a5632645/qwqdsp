# ------------------------------------------------------------
# Pitch (GUI)
# ------------------------------------------------------------
if(QWQDSP_USE_EIGEN)
    add_qwqdsp_gui_test(swift_f0_rt swift_f0_rt/main.cpp)
    if(MSVC)
        target_compile_options(qwqdsp-gui-swift_f0_rt PUBLIC /arch:AVX2)
    else()
        target_compile_options(qwqdsp-gui-swift_f0_rt PUBLIC -mavx2 -mfma)
    endif()
endif()
