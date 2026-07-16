cmake_policy(SET CMP0076 NEW)

# ------------------------------------------------------------
# add a no_gui test/demo
# ------------------------------------------------------------
function(add_qwqdsp_dsp_test name code_file)
    set(multiValueArgs EXTRA_LIBS)
    cmake_parse_arguments(ARG "" "" "${multiValueArgs}" ${ARGN})
    add_executable(qwqdsp-${name}
        ${code_file}
    )
    target_link_libraries(qwqdsp-${name} PUBLIC qwqdsp qwqdsp_support ${ARG_EXTRA_LIBS})
    target_compile_definitions(qwqdsp-${name} PRIVATE
        QWQDSP_WORK_DIR="${CMAKE_CURRENT_SOURCE_DIR}/work_dir"
    )
    set_target_properties(qwqdsp-${name} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-${name} PROPERTIES FOLDER qwqdsp-dsp-tests)
endfunction()

# ------------------------------------------------------------
# add a gui test/demo
# ------------------------------------------------------------
function(add_qwqdsp_gui_test name code_file)
    set(options NO_MINIAUDIO)
    cmake_parse_arguments(ARG "${options}" "" "" ${ARGN})

    if(ARG_NO_MINIAUDIO)
        add_executable(qwqdsp-${name}
            ${code_file}
        )
    else()
        add_executable(qwqdsp-${name}
            ${code_file}
            ../support/miniaudio.c
        )
    endif()

    target_link_libraries(qwqdsp-${name} PUBLIC qwqdsp raylib qwqdsp_support)

    if(TARGET qwqdsp_support_slider)
        target_link_libraries(qwqdsp-${name} PUBLIC qwqdsp_support_slider)
    endif()

    target_compile_definitions(qwqdsp-${name} PRIVATE
        QWQDSP_WORK_DIR="${CMAKE_CURRENT_SOURCE_DIR}/work_dir"
    )
    set_target_properties(qwqdsp-${name} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-${name} PROPERTIES FOLDER qwqdsp-gui-tests)
endfunction()

# ------------------------------------------------------------
# no gui test/demo
# ------------------------------------------------------------

# ── DSP Tests (no GUI) ──
add_qwqdsp_dsp_test(cephes nogui/cephes/cephes.cpp)

# ── Filter ──
add_qwqdsp_dsp_test(biquad nogui/filter/biquad.cpp)
add_qwqdsp_dsp_test(fir_design nogui/filter/fir_design.cpp)
add_qwqdsp_dsp_test(minimum_phase_fir nogui/filter/minimum_phase_fir.cpp)
add_qwqdsp_dsp_test(paralle_allpass nogui/filter/paralle_allpass.cpp)
add_qwqdsp_dsp_test(reverse_iir nogui/filter/reverse_iir.cpp)
add_qwqdsp_dsp_test(residual_chebyshev nogui/filter/residual_chebyshev.cpp)
add_qwqdsp_dsp_test(any_hilbert nogui/filter/any_hilbert.cpp)

# ── Audio FX ──
add_qwqdsp_dsp_test(convolution nogui/audiofx/convolution.cpp)
add_qwqdsp_dsp_test(conv_feedback nogui/audiofx/conv_feedback.cpp)
add_qwqdsp_dsp_test(conv_feedback2 nogui/audiofx/conv_feedback2.cpp)
add_qwqdsp_dsp_test(reverb nogui/audiofx/reverb.cpp)
add_qwqdsp_dsp_test(freq_shifter nogui/audiofx/freq_shifter.cpp)
add_qwqdsp_dsp_test(resample nogui/audiofx/resample.cpp)
add_qwqdsp_dsp_test(oversample nogui/audiofx/oversample.cpp)
add_qwqdsp_dsp_test(pitch_shifter nogui/audiofx/pitch_shifter.cpp)
if(QWQDSP_USE_IPP)
    add_qwqdsp_dsp_test(auto_notch nogui/audiofx/auto_notch.cpp)
    add_qwqdsp_dsp_test(wiener_filter nogui/audiofx/wiener_filter.cpp)
    add_qwqdsp_dsp_test(piwarp nogui/audiofx/piwarp.cpp)
endif()

# ── Spectral ──
add_qwqdsp_dsp_test(fft_interpolation nogui/spectral/fft_interpolation.cpp)
add_qwqdsp_dsp_test(ifft_random_phase nogui/spectral/ifft_random_phase.cpp)
add_qwqdsp_dsp_test(fft_fold nogui/spectral/fold.cpp)


# ── Pitch ──
add_qwqdsp_dsp_test(pyin nogui/pitch/pyin.cpp)
add_qwqdsp_dsp_test(pmpm nogui/pitch/pmpm.cpp)
add_qwqdsp_dsp_test(pmpm-vs-pyin nogui/pitch/pmpm_vs_pyin.cpp)
add_qwqdsp_dsp_test(pitch_algorithms nogui/pitch/pitch_algorithms.cpp)
add_qwqdsp_dsp_test(pyin-vis nogui/pitch/pyin_vis.cpp EXTRA_LIBS qwqdsp_support_stb_image)
if(QWQDSP_USE_EIGEN)
    add_qwqdsp_dsp_test(swift_f0 nogui/pitch/swift_f0.cpp)
    add_qwqdsp_dsp_test(swift_f0-vis nogui/pitch/swift_f0_vis.cpp EXTRA_LIBS qwqdsp_support_stb_image)
endif()

# ── Synth ──
add_qwqdsp_dsp_test(dsf nogui/synth/dsf.cpp)
add_qwqdsp_dsp_test(stupid_resynthsis nogui/synth/stupid_resynthsis.cpp)
add_qwqdsp_dsp_test(wsola nogui/synth/wsola.cpp)

# ----- rnn noise lite -----
file(GLOB_RECURSE qwqdsp-rnnoise-lite-src "rnnoise_lite/nnom/*.c")
add_executable(qwqdsp-rnnoise-lite
    ${qwqdsp-rnnoise-lite-src}
    rnnoise_lite/main.cpp
    rnnoise_lite/mfcc.c
    rnnoise_lite/rnnoise.c
)
target_link_libraries(qwqdsp-rnnoise-lite PUBLIC qwqdsp raylib qwqdsp_support)
target_compile_definitions(qwqdsp-rnnoise-lite PRIVATE
    QWQDSP_WORK_DIR="${CMAKE_CURRENT_SOURCE_DIR}/work_dir"
)
set_target_properties(qwqdsp-rnnoise-lite PROPERTIES CXX_STANDARD 20)
set_target_properties(qwqdsp-rnnoise-lite PROPERTIES FOLDER qwqdsp-dsp-tests)

# ------------------------------------------------------------
# gui test/demo
# ------------------------------------------------------------
if(QWQDSP_USE_RAYLIB)
    # ── GUI Tests ──

    # ── General ──
    add_qwqdsp_gui_test(echo gui/audiofx/echo.cpp)
    add_qwqdsp_gui_test(reverb_rt gui/audiofx/reverb.cpp)

    add_qwqdsp_gui_test(reverb_shimmer gui/audiofx/reverb.cpp)

    add_qwqdsp_gui_test(interpolations gui/interpolations.cpp)
    if(QWQDSP_USE_IPP)
        add_qwqdsp_gui_test(auto_notch_rt gui/auto_notch_rt.cpp)
    endif()
    add_qwqdsp_gui_test(filters gui/audiofx/filters.cpp NO_MINIAUDIO)
    add_qwqdsp_gui_test(filters2 gui/audiofx/filters2.cpp NO_MINIAUDIO)

    # ── Spectral (GUI) ──
    add_qwqdsp_gui_test(pfb_view gui/spectral/polyphase_filter_bank_view.cpp)
    target_compile_options(qwqdsp-pfb_view PUBLIC -mavx2 -mfma)

    add_qwqdsp_gui_test(pfb2_view gui/spectral/polyphase_filter_bank_view2.cpp)
    target_compile_options(qwqdsp-pfb2_view PUBLIC -mavx2 -mfma)

    add_qwqdsp_gui_test(resonate_bank gui/spectral/resonate_bank.cpp)
    target_compile_options(qwqdsp-resonate_bank PUBLIC -mavx2 -mfma)

    if(QWQDSP_USE_IPP)
        add_qwqdsp_gui_test(spectrum gui/spectral/spectrum/spectrum.cpp)
        target_sources(qwqdsp-spectrum PUBLIC gui/spectral/spectrum/spectrum_driver.cpp)

        add_qwqdsp_gui_test(spectrum2 gui/spectral/spectrum2/spectrum.cpp)
        target_sources(qwqdsp-spectrum2 PUBLIC gui/spectral/spectrum2/spectrum_driver.cpp)
    endif()

    add_qwqdsp_gui_test(spectrum3 gui/spectral/spectrum3.cpp)
    add_qwqdsp_gui_test(spectrum4 gui/spectral/spectrum4.cpp)
    add_qwqdsp_gui_test(reassignment gui/spectral/reassignment.cpp)

    # ── Waveform (GUI) ──
    add_qwqdsp_gui_test(waveform2 gui/waveform2.cpp)

    # ── Synth (GUI) ──
    add_qwqdsp_gui_test(formant_filter gui/synth/formant_filter.cpp NO_MINIAUDIO)
    add_qwqdsp_gui_test(blit gui/synth/blit.cpp NO_MINIAUDIO)
    add_qwqdsp_gui_test(polyblep gui/synth/polyblep.cpp NO_MINIAUDIO)
    add_qwqdsp_gui_test(noise gui/synth/noise.cpp NO_MINIAUDIO)
endif()
