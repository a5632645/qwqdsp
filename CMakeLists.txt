cmake_minimum_required(VERSION 3.11)
project(qwqdsp)

option(QWQDSP_BUILD_TESTS "Disable building tests" ON)

# check for intel ipp
set(IPP_SHARED OFF)
find_package(IPP)
set(qwqdsp_have_ipp ${IPP_FOUND})

# main library
add_library(qwqdsp
    STATIC
        "source/oouras.cpp"
        "source/real_fft.cpp"
        "source/complex_fft.cpp"
)
target_include_directories(qwqdsp PUBLIC "./include")
set_target_properties(qwqdsp PROPERTIES CXX_STANDARD 20)

# external library
if (TARGET eigen)
    target_link_libraries(qwqdsp PUBLIC Eigen3::Eigen)
    target_compile_definitions(qwqdsp PUBLIC QWQDSP_HAVE_EIGEN=1)
endif()
if (qwqdsp_have_ipp)
    target_link_libraries(qwqdsp PUBLIC IPP::ipps)
    target_compile_definitions(qwqdsp PUBLIC QWQDSP_HAVE_IPP=1)
endif()
if (TARGET simde)
    target_link_libraries(qwqdsp PUBLIC simde)
    target_compile_definitions(qwqdsp PUBLIC QWQDSP_HAVE_SIMDE)
endif()

add_subdirectory(support)

# enable tests
if (QWQDSP_BUILD_TESTS)
    if (QWQDSP_HAVE_RAYLIB)
        add_executable(playing
            "playing/playing.cpp"
        )
        target_link_libraries(playing PUBLIC qwqdsp raylib qwqdsp_support)
        target_link_libraries(playing PRIVATE qwqdsp_support_slider qwqdsp_support)
        set_target_properties(playing PROPERTIES CXX_STANDARD 20)

        if (MSVC)
            target_compile_options(playing PUBLIC /arch:AVX2)
        else()
            target_compile_options(playing PUBLIC -mavx2)
        endif()
    endif()
    
    add_subdirectory(tests)
endif()
