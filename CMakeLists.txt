cmake_minimum_required(VERSION 3.22)
project(qwqdsp)

# 依赖由宏控制，外部链接
# QWQDSP_HAVE_IPP (需外部链接ipp)
# QWQDSP_HAVE_ACCELERATE (-framework Accelerate)

# 依赖由外部 CMake 控制:
# QWQDSP_USE_EIGEN  — 定义QWQDSP_HAVE_EIGEN  (需外部链接 qwqdsp_eigen)
# QWQDSP_USE_SIMDE  — 定义QWQDSP_HAVE_SIMDE  (需外部链接 simde)
# QWQDSP_USE_RAYLIB — playing 和 GUI测试 依赖 (需外部提供 raylib)
# QWQDSP_BUILD_TESTS

# main library
add_library(qwqdsp INTERFACE)
target_include_directories(qwqdsp INTERFACE "./include")
set_target_properties(qwqdsp PROPERTIES CXX_STANDARD 20)

# 仅定义编译宏 — 链接由外部 CMake 负责
if(QWQDSP_USE_EIGEN)
    target_compile_definitions(qwqdsp INTERFACE QWQDSP_HAVE_EIGEN)
endif()

if(QWQDSP_USE_SIMDE)
    target_compile_definitions(qwqdsp PUBLIC QWQDSP_HAVE_SIMDE)
endif()

# enable tests
if(QWQDSP_BUILD_TESTS)
    add_subdirectory(support)
    add_subdirectory(tests)
endif()
