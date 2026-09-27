add_executable(qwqdsp-lab-playing
    "playing.cpp"
    "${QWQDSP_SRC_DIR}/support/miniaudio.c"
)

target_link_libraries(qwqdsp-lab-playing PUBLIC qwqdsp raylib qwqdsp_support)
target_link_libraries(qwqdsp-lab-playing PRIVATE qwqdsp_support_slider qwqdsp_support)
set_target_properties(qwqdsp-lab-playing PROPERTIES CXX_STANDARD 20)
set_target_properties(qwqdsp-lab-playing PROPERTIES FOLDER qwqdsp-labs)

if(MSVC)
    target_compile_options(qwqdsp-lab-playing PUBLIC /arch:AVX2)
else()
    target_compile_options(qwqdsp-lab-playing PUBLIC -mavx2)
endif()
