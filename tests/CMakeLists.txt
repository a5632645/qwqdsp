cmake_policy(SET CMP0076 NEW)

function(add_qwqdsp_dsp_test ex_file)
    add_executable(qwqdsp-${ex_file}
        ${ex_file}.cpp
    )
    target_link_libraries(qwqdsp-${ex_file} PUBLIC qwqdsp)
    set_target_properties(qwqdsp-${ex_file} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-${ex_file} PROPERTIES FOLDER qwqdsp-dsp-tests)
endfunction()

function(add_qwqdsp_gui_test ex_file)
    add_executable(qwqdsp-${ex_file}
        ${ex_file}.cpp
    )
    target_link_libraries(qwqdsp-${ex_file} PUBLIC qwqdsp raylib)
    set_target_properties(qwqdsp-${ex_file} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-${ex_file} PROPERTIES FOLDER qwqdsp-gui-tests)
endfunction()

add_qwqdsp_dsp_test(convolution)
add_qwqdsp_dsp_test(fft_interpolation)
add_qwqdsp_dsp_test(minimum_phase_fir)
add_qwqdsp_dsp_test(freq_shifter)
add_qwqdsp_dsp_test(fir_design)
add_qwqdsp_dsp_test(dsf)
add_qwqdsp_dsp_test(resample)
add_qwqdsp_dsp_test(biquad)
add_qwqdsp_dsp_test(paralle_allpass)
add_qwqdsp_dsp_test(conv_feedback)
add_qwqdsp_dsp_test(reverse_iir)
add_qwqdsp_dsp_test(ifft_random_phase)
add_qwqdsp_dsp_test(wiener_filter)
add_qwqdsp_dsp_test(auto_notch)

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
    add_qwqdsp_gui_test(interpolations)

    add_qwqdsp_gui_test(formant_filter)
    target_sources(qwqdsp-formant_filter PUBLIC "../playing/slider.cpp")

    add_qwqdsp_gui_test(auto_notch_rt)
    target_sources(qwqdsp-auto_notch_rt PUBLIC "../playing/slider.cpp")
    target_sources(qwqdsp-auto_notch_rt PUBLIC "../playing/miniaudio.c")

    add_qwqdsp_gui_test(blit)
    target_sources(qwqdsp-blit PUBLIC "../playing/slider.cpp")

    add_qwqdsp_gui_test(polyblep)
    target_sources(qwqdsp-polyblep PUBLIC "../playing/slider.cpp")

    add_qwqdsp_gui_test(resonate_bank)
    target_sources(qwqdsp-resonate_bank PUBLIC "../playing/miniaudio.c")
    target_compile_options(qwqdsp-resonate_bank PUBLIC -mavx2 -mfma)

    add_qwqdsp_gui_test(spectrum)
    target_sources(qwqdsp-spectrum PUBLIC s/spectrum_driver.cpp ../playing/miniaudio.c)
endif()
