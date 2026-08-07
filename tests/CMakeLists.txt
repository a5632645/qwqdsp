cmake_policy(SET CMP0076 NEW)

# ------------------------------------------------------------
# 测试源码根目录（供子目录调用的 add_qwqdsp_* 定位 work_dir / support）
# ------------------------------------------------------------
set(QWQDSP_TESTS_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# ------------------------------------------------------------
# add a no_gui test/demo
# ------------------------------------------------------------
function(add_qwqdsp_dsp_test name code_file)
    set(multiValueArgs EXTRA_LIBS)
    cmake_parse_arguments(ARG "" "" "${multiValueArgs}" ${ARGN})
    add_executable(qwqdsp-nogui-${name}
        ${code_file}
    )
    target_link_libraries(qwqdsp-nogui-${name} PUBLIC qwqdsp qwqdsp_support ${ARG_EXTRA_LIBS})
    target_compile_definitions(qwqdsp-nogui-${name} PRIVATE
        QWQDSP_WORK_DIR="${QWQDSP_TESTS_SRC_DIR}/work_dir"
    )
    set_target_properties(qwqdsp-nogui-${name} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-nogui-${name} PROPERTIES FOLDER qwqdsp-dsp-tests)
endfunction()

# ------------------------------------------------------------
# add a gui test/demo
# ------------------------------------------------------------
function(add_qwqdsp_gui_test name code_file)
    set(options NO_MINIAUDIO)
    cmake_parse_arguments(ARG "${options}" "" "" ${ARGN})

    if(ARG_NO_MINIAUDIO)
        add_executable(qwqdsp-gui-${name}
            ${code_file}
        )
    else()
        add_executable(qwqdsp-gui-${name}
            ${code_file}
            ${QWQDSP_TESTS_SRC_DIR}/../support/miniaudio.c
        )
    endif()

    target_link_libraries(qwqdsp-gui-${name} PUBLIC qwqdsp raylib qwqdsp_support)

    if(TARGET qwqdsp_support_slider)
        target_link_libraries(qwqdsp-gui-${name} PUBLIC qwqdsp_support_slider)
    endif()

    target_compile_definitions(qwqdsp-gui-${name} PRIVATE
        QWQDSP_WORK_DIR="${QWQDSP_TESTS_SRC_DIR}/work_dir"
    )
    set_target_properties(qwqdsp-gui-${name} PROPERTIES CXX_STANDARD 20)
    set_target_properties(qwqdsp-gui-${name} PROPERTIES FOLDER qwqdsp-gui-tests)
endfunction()

# ------------------------------------------------------------
# 分发到子文件夹
# ------------------------------------------------------------

# ----- DSP 测试（无 GUI） -----
add_subdirectory(nogui)

# ----- GUI 测试 -----
add_subdirectory(gui)
