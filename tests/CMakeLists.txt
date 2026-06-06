cmake_policy(SET CMP0076 NEW)

function(add_qwqdsp_dsp_test name code_file)
    add_executable(qwqdsp-${name}
        ${code_file}
    )
    target_link_libraries(qwqdsp-${name} PUBLIC qwqdsp)
    set_target_properties(qwqdsp-${name} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-${name} PROPERTIES FOLDER qwqdsp-dsp-tests)
endfunction()

function(add_qwqdsp_gui_test name code_file)
    add_executable(qwqdsp-${name}
        ${code_file}
    )
    target_link_libraries(qwqdsp-${name} PUBLIC qwqdsp raylib)
    set_target_properties(qwqdsp-${name} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-${name} PROPERTIES FOLDER qwqdsp-gui-tests)
endfunction()

add_qwqdsp_dsp_test(convolution convolution.cpp)
add_qwqdsp_dsp_test(fft_interpolation fft_interpolation.cpp)
add_qwqdsp_dsp_test(minimum_phase_fir minimum_phase_fir.cpp)
add_qwqdsp_dsp_test(freq_shifter freq_shifter.cpp)
add_qwqdsp_dsp_test(fir_design fir_design.cpp)
add_qwqdsp_dsp_test(dsf dsf.cpp)
add_qwqdsp_dsp_test(resample resample.cpp)
add_qwqdsp_dsp_test(biquad biquad.cpp)
add_qwqdsp_dsp_test(paralle_allpass paralle_allpass.cpp)
add_qwqdsp_dsp_test(conv_feedback conv_feedback.cpp)
add_qwqdsp_dsp_test(reverse_iir reverse_iir.cpp)
add_qwqdsp_dsp_test(ifft_random_phase ifft_random_phase.cpp)
add_qwqdsp_dsp_test(wiener_filter wiener_filter.cpp)
add_qwqdsp_dsp_test(auto_notch auto_notch.cpp)

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
    add_qwqdsp_gui_test(interpolations interpolations.cpp)

    add_qwqdsp_gui_test(formant_filter formant_filter.cpp)
    target_sources(qwqdsp-formant_filter PUBLIC "../playing/slider.cpp")

    add_qwqdsp_gui_test(auto_notch_rt auto_notch_rt.cpp)
    target_sources(qwqdsp-auto_notch_rt PUBLIC "../playing/slider.cpp")
    target_sources(qwqdsp-auto_notch_rt PUBLIC "../playing/miniaudio.c")

    add_qwqdsp_gui_test(blit blit.cpp)
    target_sources(qwqdsp-blit PUBLIC "../playing/slider.cpp")

    add_qwqdsp_gui_test(polyblep polyblep.cpp)
    target_sources(qwqdsp-polyblep PUBLIC "../playing/slider.cpp")

    add_qwqdsp_gui_test(resonate_bank resonate_bank.cpp)
    target_sources(qwqdsp-resonate_bank PUBLIC "../playing/miniaudio.c")
    target_compile_options(qwqdsp-resonate_bank PUBLIC -mavx2 -mfma)

    add_qwqdsp_gui_test(spectrum s/spectrum.cpp)
    target_sources(qwqdsp-spectrum PUBLIC s/spectrum_driver.cpp ../playing/miniaudio.c)

    add_qwqdsp_gui_test(spectrum2 s2/spectrum.cpp)
    target_sources(qwqdsp-spectrum2 PUBLIC s2/spectrum_driver.cpp ../playing/miniaudio.c)
endif()
