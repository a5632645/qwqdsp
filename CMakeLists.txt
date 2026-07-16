cmake_minimum_required(VERSION 3.22)
project(qwqdsp)

option(QWQDSP_BUILD_TESTS "building tests" ON)

# 依赖由外部 CMake 控制:
# QWQDSP_USE_EIGEN  — 定义 QWQDSP_HAVE_EIGEN  (需外部链接 qwqdsp_eigen)
# QWQDSP_USE_IPP    — 定义 QWQDSP_HAVE_IPP    (需外部链接 IPP::ipps)
# QWQDSP_USE_SIMDE  — 定义 QWQDSP_HAVE_SIMDE  (需外部链接 simde)
# QWQDSP_USE_RAYLIB — 构建 playing 和 GUI 测试 (需外部提供 raylib)

# main library
add_library(qwqdsp
    STATIC
    "source/oouras.cpp"
)
target_include_directories(qwqdsp PUBLIC "./include")
set_target_properties(qwqdsp PROPERTIES CXX_STANDARD 20)

# 仅定义编译宏 — 链接由外部 CMake 负责
if(QWQDSP_USE_EIGEN)
    target_compile_definitions(qwqdsp PUBLIC QWQDSP_HAVE_EIGEN)
endif()

if(QWQDSP_USE_IPP)
    target_compile_definitions(qwqdsp PUBLIC QWQDSP_HAVE_IPP)
    message(STATUS "[qwqdsp] using intel ipp")
elseif(QWQDSP_USE_ACCELERATE)
    target_compile_definitions(qwqdsp PUBLIC QWQDSP_HAVE_ACCELERATE)
    message(STATUS "[qwqdsp] using Apple Accelerate")
else()
    message(STATUS "[qwqdsp] no intel ipp/accelerate, fallback to oouras")
endif()

if(QWQDSP_USE_SIMDE)
    target_compile_definitions(qwqdsp PUBLIC QWQDSP_HAVE_SIMDE)
endif()

# enable tests
if(QWQDSP_BUILD_TESTS)
    add_subdirectory(support)
    add_subdirectory(tests)

    if(QWQDSP_USE_RAYLIB)
        add_executable(playing
            "playing/playing.cpp"
        )
        target_link_libraries(playing PUBLIC qwqdsp raylib qwqdsp_support)
        target_link_libraries(playing PRIVATE qwqdsp_support_slider qwqdsp_support)
        set_target_properties(playing PROPERTIES CXX_STANDARD 20)

        if(MSVC)
            target_compile_options(playing PUBLIC /arch:AVX2)
        else()
            target_compile_options(playing PUBLIC -mavx2)
        endif()
    endif()
endif()
