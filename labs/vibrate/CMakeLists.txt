add_executable(qwqdsp-lab-vibrate-chain1d "chain_1d.cpp")
target_link_libraries(qwqdsp-lab-vibrate-chain1d PRIVATE raylib)
set_target_properties(qwqdsp-lab-vibrate-chain1d PROPERTIES CXX_STANDARD 20 FOLDER qwqdsp-labs)

add_executable(qwqdsp-lab-vibrate-cross "cross_strings.cpp")
target_link_libraries(qwqdsp-lab-vibrate-cross PRIVATE raylib)
set_target_properties(qwqdsp-lab-vibrate-cross PROPERTIES CXX_STANDARD 20 FOLDER qwqdsp-labs)

add_executable(qwqdsp-lab-vibrate-membrane "membrane.cpp")
target_link_libraries(qwqdsp-lab-vibrate-membrane PRIVATE raylib)
set_target_properties(qwqdsp-lab-vibrate-membrane PROPERTIES CXX_STANDARD 20 FOLDER qwqdsp-labs)

add_executable(qwqdsp-lab-vibrate-disk "disk.cpp")
target_link_libraries(qwqdsp-lab-vibrate-disk PRIVATE raylib)
set_target_properties(qwqdsp-lab-vibrate-disk PROPERTIES CXX_STANDARD 20 FOLDER qwqdsp-labs)

add_executable(qwqdsp-lab-vibrate-ring "ring_strings.cpp")
target_link_libraries(qwqdsp-lab-vibrate-ring PRIVATE raylib)
set_target_properties(qwqdsp-lab-vibrate-ring PROPERTIES CXX_STANDARD 20 FOLDER qwqdsp-labs)

add_executable(qwqdsp-lab-vibrate-waterphone "waterphone.cpp")
target_link_libraries(qwqdsp-lab-vibrate-waterphone PRIVATE raylib)
set_target_properties(qwqdsp-lab-vibrate-waterphone PROPERTIES CXX_STANDARD 20 FOLDER qwqdsp-labs)
