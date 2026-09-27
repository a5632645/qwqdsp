# ------------------------------------------------------------
# 目标构造器：tests / example / labs 三个桶共用
#
# 需要调用方先设置 QWQDSP_SRC_DIR 为仓库根（定位 work_dir / support）。
# 命名规则：目标名 = 桶前缀 + 文件名，产物落在 build 里与源码同构的目录。
#   qwqdsp-test-<name>     无头测试（自己判定是否通过）
#   qwqdsp-example-<name>  无头 example（导出 WAV/PNG 给人试听试看）
#   qwqdsp-example-<name>  GUI example（开 raylib 窗口，见 add_qwqdsp_gui_example）
#   qwqdsp-lab-<name>      labs 探索程序
# ------------------------------------------------------------

# ----- 无头：test / example / lab 共用的骨架 -----
function(_qwqdsp_headless prefix folder name code_file)
    set(multiValueArgs EXTRA_LIBS)
    cmake_parse_arguments(ARG "" "" "${multiValueArgs}" ${ARGN})

    set(target ${prefix}-${name})
    add_executable(${target}
        ${code_file}
    )
    target_link_libraries(${target} PUBLIC qwqdsp qwqdsp_support ${ARG_EXTRA_LIBS})
    target_compile_definitions(${target} PRIVATE
        QWQDSP_WORK_DIR="${QWQDSP_SRC_DIR}/work_dir"
    )
    set_target_properties(${target} PROPERTIES CXX_STANDARD 20)
    set_target_properties(${target} PROPERTIES FOLDER ${folder})
endfunction()

# ------------------------------------------------------------
# 无头测试
# ------------------------------------------------------------
function(add_qwqdsp_test name code_file)
    _qwqdsp_headless(qwqdsp-test qwqdsp-tests ${name} ${code_file} ${ARGN})
endfunction()

# ------------------------------------------------------------
# 无头 example（导出音频/图片供人查看）
# ------------------------------------------------------------
function(add_qwqdsp_example name code_file)
    _qwqdsp_headless(qwqdsp-example qwqdsp-examples ${name} ${code_file} ${ARGN})
endfunction()

# ------------------------------------------------------------
# labs 探索程序
# ------------------------------------------------------------
function(add_qwqdsp_lab name code_file)
    _qwqdsp_headless(qwqdsp-lab qwqdsp-labs ${name} ${code_file} ${ARGN})
endfunction()

# ------------------------------------------------------------
# GUI example（raylib；默认再编译一份 miniaudio）
# ------------------------------------------------------------
function(add_qwqdsp_gui_example name code_file)
    set(options NO_MINIAUDIO)
    cmake_parse_arguments(ARG "${options}" "" "" ${ARGN})

    set(target qwqdsp-example-${name})
    if(ARG_NO_MINIAUDIO)
        add_executable(${target}
            ${code_file}
        )
    else()
        add_executable(${target}
            ${code_file}
            ${QWQDSP_SRC_DIR}/support/miniaudio.c
        )
    endif()

    target_link_libraries(${target} PUBLIC qwqdsp raylib qwqdsp_support)

    if(TARGET qwqdsp_support_slider)
        target_link_libraries(${target} PUBLIC qwqdsp_support_slider)
    endif()

    target_compile_definitions(${target} PRIVATE
        QWQDSP_WORK_DIR="${QWQDSP_SRC_DIR}/work_dir"
    )
    set_target_properties(${target} PROPERTIES CXX_STANDARD 20)
    set_target_properties(${target} PROPERTIES FOLDER qwqdsp-examples-gui)
endfunction()
