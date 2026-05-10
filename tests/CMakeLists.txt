cmake_policy(SET CMP0076 NEW)

enable_testing()

function(add_qwqdsp_test ex_file)
    add_executable(qwqdsp-${ex_file}
        ${ex_file}.cpp
    )
    target_link_libraries(qwqdsp-${ex_file} PUBLIC qwqdsp raylib)
    set_target_properties(qwqdsp-${ex_file} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-${ex_file} PROPERTIES FOLDER qwqdsp-tests)
endfunction()

add_qwqdsp_test(convolution)
add_qwqdsp_test(fft_interpolation)
add_qwqdsp_test(interpolations)
add_qwqdsp_test(minimum_phase_fir)
add_qwqdsp_test(freq_shifter)
add_qwqdsp_test(fir_design)
add_qwqdsp_test(formant_filter)
target_sources(qwqdsp-formant_filter PUBLIC "../playing/slider.cpp")
add_qwqdsp_test(dsf)
add_qwqdsp_test(blit)
target_sources(qwqdsp-blit PUBLIC "../playing/slider.cpp")
add_qwqdsp_test(polyblep)
target_sources(qwqdsp-polyblep PUBLIC "../playing/slider.cpp")
add_qwqdsp_test(resample)
add_qwqdsp_test(biquad)
add_qwqdsp_test(paralle_allpass)

file(GLOB_RECURSE qwqdsp-rnnoise-lite-src "rnnoise_lite/nnom/*.c")
add_executable(qwqdsp-rnnoise-lite
    ${qwqdsp-rnnoise-lite-src}
    rnnoise_lite/main.cpp
    rnnoise_lite/mfcc.c
    rnnoise_lite/rnnoise.c
)
target_link_libraries(qwqdsp-rnnoise-lite PUBLIC qwqdsp raylib)
set_target_properties(qwqdsp-rnnoise-lite PROPERTIES CXX_STANDARD 20)
set_target_properties(qwqdsp-rnnoise-lite PROPERTIES FOLDER qwqdsp-tests)

add_qwqdsp_test(conv_feedback)
add_qwqdsp_test(reverse_iir)
add_qwqdsp_test(ifft_random_phase)
