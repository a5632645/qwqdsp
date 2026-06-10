cmake_policy(SET CMP0076 NEW)

function(add_qwqdsp_dsp_test name code_file)
    add_executable(qwqdsp-${name}
        ${code_file}
    )
    target_link_libraries(qwqdsp-${name} PUBLIC qwqdsp)
    target_include_directories(qwqdsp-${name} PUBLIC .)
    set_target_properties(qwqdsp-${name} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-${name} PROPERTIES FOLDER qwqdsp-dsp-tests)
endfunction()

function(add_qwqdsp_gui_test name code_file)
    add_executable(qwqdsp-${name}
        ${code_file}
    )
    target_link_libraries(qwqdsp-${name} PUBLIC qwqdsp raylib)
    target_include_directories(qwqdsp-${name} PUBLIC .)
    set_target_properties(qwqdsp-${name} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-${name} PROPERTIES FOLDER qwqdsp-gui-tests)
endfunction()

# ── Signal Processing ──
add_qwqdsp_dsp_test(biquad signal/biquad.cpp)
add_qwqdsp_dsp_test(fir_design signal/fir_design.cpp)
add_qwqdsp_dsp_test(wiener_filter signal/wiener_filter.cpp)
add_qwqdsp_dsp_test(auto_notch signal/auto_notch.cpp)
add_qwqdsp_dsp_test(minimum_phase_fir signal/minimum_phase_fir.cpp)
add_qwqdsp_dsp_test(paralle_allpass signal/paralle_allpass.cpp)
add_qwqdsp_dsp_test(reverse_iir signal/reverse_iir.cpp)
add_qwqdsp_dsp_test(convolution signal/convolution.cpp)
add_qwqdsp_dsp_test(fft_interpolation signal/fft_interpolation.cpp)
add_qwqdsp_dsp_test(ifft_random_phase signal/ifft_random_phase.cpp)

# ── Audio FX ──
add_qwqdsp_dsp_test(freq_shifter audiofx/freq_shifter.cpp)
add_qwqdsp_dsp_test(conv_feedback audiofx/conv_feedback.cpp)
add_qwqdsp_dsp_test(dsf audiofx/dsf.cpp)

# ── Resample ──
add_qwqdsp_dsp_test(resample resample/resample.cpp)
add_qwqdsp_dsp_test(polyphase resample/polyphase.cpp)

file(GLOB_RECURSE qwqdsp-rnnoise-lite-src "rnnoise_lite/nnom/*.c")
add_executable(qwqdsp-rnnoise-lite
    ${qwqdsp-rnnoise-lite-src}
    rnnoise_lite/main.cpp
    rnnoise_lite/mfcc.c
    rnnoise_lite/rnnoise.c
)
target_link_libraries(qwqdsp-rnnoise-lite PUBLIC qwqdsp raylib)
set_target_properties(qwqdsp-rnnoise-lite PROPERTIES CXX_STANDARD 20)
set_target_properties(qwqdsp-rnnoise-lite PROPERTIES FOLDER qwqdsp-dsp-tests)

if (QWQDSP_HAVE_RAYLIB)
    # ── Realtime GUI ──
    add_qwqdsp_gui_test(interpolations realtime/interpolations.cpp)

    add_qwqdsp_gui_test(auto_notch_rt realtime/auto_notch_rt.cpp)
    target_sources(qwqdsp-auto_notch_rt PUBLIC "../playing/slider.cpp")
    target_sources(qwqdsp-auto_notch_rt PUBLIC "../playing/miniaudio.c")

    # ── Audio FX (GUI) ──
    add_qwqdsp_gui_test(formant_filter audiofx/formant_filter.cpp)
    target_sources(qwqdsp-formant_filter PUBLIC "../playing/slider.cpp")

    add_qwqdsp_gui_test(blit audiofx/blit.cpp)
    target_sources(qwqdsp-blit PUBLIC "../playing/slider.cpp")

    add_qwqdsp_gui_test(polyblep audiofx/polyblep.cpp)
    target_sources(qwqdsp-polyblep PUBLIC "../playing/slider.cpp")

    add_qwqdsp_gui_test(resonate_bank audiofx/resonate_bank.cpp)
    target_sources(qwqdsp-resonate_bank PUBLIC "../playing/miniaudio.c")
    target_compile_options(qwqdsp-resonate_bank PUBLIC -mavx2 -mfma)

    # ── Spectral ──
    add_qwqdsp_gui_test(spectrum spectral/spectrum/spectrum.cpp)
    target_sources(qwqdsp-spectrum PUBLIC spectral/spectrum/spectrum_driver.cpp ../playing/miniaudio.c)

    add_qwqdsp_gui_test(spectrum2 spectral/spectrum2/spectrum.cpp)
    target_sources(qwqdsp-spectrum2 PUBLIC spectral/spectrum2/spectrum_driver.cpp ../playing/miniaudio.c)

    add_qwqdsp_gui_test(spectrum3 spectral/spectrum3.cpp)
    target_sources(qwqdsp-spectrum3 PUBLIC ../playing/miniaudio.c)
endif()
