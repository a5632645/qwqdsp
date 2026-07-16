add_library(qwqdsp_support INTERFACE)
target_include_directories(qwqdsp_support INTERFACE .)

# stb_image_write 图像输出库
add_library(qwqdsp_support_stb_image
    OBJECT
    "stb_image_write.c"
)
target_include_directories(qwqdsp_support_stb_image PUBLIC .)
target_link_libraries(qwqdsp_support_stb_image PUBLIC qwqdsp_support)

if(QWQDSP_USE_RAYLIB)
    add_library(qwqdsp_support_slider OBJECT
        slider.cpp
    )
    target_include_directories(qwqdsp_support_slider PUBLIC .)
    target_link_libraries(qwqdsp_support_slider PUBLIC raylib)
    set_target_properties(qwqdsp_support_slider PROPERTIES CXX_STANDARD 20)
endif()
