add_library(qwqdsp_support INTERFACE)
target_include_directories(qwqdsp_support INTERFACE .)

if (TARGET raylib)
    add_library(qwqdsp_support_slider STATIC
        slider.cpp
    )
    target_include_directories(qwqdsp_support_slider PUBLIC .)
    target_link_libraries(qwqdsp_support_slider PUBLIC raylib)
    set_target_properties(qwqdsp_support_slider PROPERTIES CXX_STANDARD 20)
endif()
