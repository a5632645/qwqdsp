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

function(add_qwqdsp_test2 ex_file)
    add_executable(qwqdsp-${ex_file}
        ${ex_file}.cpp
    )
    target_link_libraries(qwqdsp-${ex_file} PUBLIC qwqdsp)
    set_target_properties(qwqdsp-${ex_file} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-${ex_file} PROPERTIES FOLDER qwqdsp-tests)
endfunction()

add_qwqdsp_test2(convolution)
add_qwqdsp_test2(fft_interpolation)
add_qwqdsp_test2(minimum_phase_fir)
add_qwqdsp_test2(freq_shifter)
add_qwqdsp_test2(fir_design)
add_qwqdsp_test2(dsf)
add_qwqdsp_test2(resample)
add_qwqdsp_test2(biquad)
add_qwqdsp_test2(paralle_allpass)
add_qwqdsp_test2(residual_chebyshev)

if (qwqdsp_have_raylib)
    add_qwqdsp_test(interpolations)
    add_qwqdsp_test(formant_filter)
    target_sources(qwqdsp-formant_filter PUBLIC "../playing/slider.cpp")
    add_qwqdsp_test(blit)
    target_sources(qwqdsp-blit PUBLIC "../playing/slider.cpp")
    add_qwqdsp_test(polyblep)
    target_sources(qwqdsp-polyblep PUBLIC "../playing/slider.cpp")
endif()
