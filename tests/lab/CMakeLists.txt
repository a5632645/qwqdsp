add_executable(playing
    "playing.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../support/miniaudio.c"
)

target_link_libraries(playing PUBLIC qwqdsp raylib qwqdsp_support)
target_link_libraries(playing PRIVATE qwqdsp_support_slider qwqdsp_support)
set_target_properties(playing PROPERTIES CXX_STANDARD 20)

if(MSVC)
    target_compile_options(playing PUBLIC /arch:AVX2)
else()
    target_compile_options(playing PUBLIC -mavx2)
endif()
