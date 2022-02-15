
FetchContent_GetProperties(r-index)
FetchContent_GetProperties(ksw2)

add_executable(split_fa split_fa.cpp)
target_link_libraries(split_fa klib z)

add_executable(visualizer blast-like-visualizer.cpp)
target_link_libraries(visualizer common ksw2 sdsl)
target_include_directories(visualizer PUBLIC    "../include/common"
                                                "${ksw2_SOURCE_DIR}")
target_compile_options(visualizer PUBLIC "-std=c++17")