
FetchContent_GetProperties(r-index)
FetchContent_GetProperties(ksw2)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(klib)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)

add_executable(split_fa split_fa.cpp)
target_link_libraries(split_fa klib z)

add_executable(visualizer blast-like-visualizer.cpp)
target_link_libraries(visualizer common ksw2 sdsl malloc_count)
target_include_directories(visualizer PUBLIC    "../include/common"
                                                "${ksw2_SOURCE_DIR}")
target_compile_options(visualizer PUBLIC "-std=c++17")

add_executable(ksw_ksw2_global ksw_ksw2_global.cpp ${klib_SOURCE_DIR}/ksw.c)
target_link_libraries(ksw_ksw2_global common ksw2 sdsl klib z malloc_count)
target_include_directories(ksw_ksw2_global PUBLIC    "../include/common"
                                                "${ksw2_SOURCE_DIR}"
                                                "${klib_SOURCE_DIR}"
                                                )
target_compile_options(ksw_ksw2_global PUBLIC "-std=c++17")

add_executable(ra random_access.cpp)
target_link_libraries(ra common sdsl lvsam klib z malloc_count)
target_include_directories(ra PUBLIC    "../include/common"
                                                "../include/aligner"
                                                "${shaped_slp_SOURCE_DIR}" 
                                                "${FOLCA_SOURCE_DIR}" 
                                                "${SUX_SOURCE_DIR}/function" 
                                                "${SUX_SOURCE_DIR}/support"
                                                )
target_compile_options(ra PUBLIC "-std=c++17")

add_executable(extract extract.cpp)
target_link_libraries(extract common sdsl lvsam klib z malloc_count)
target_include_directories(extract PUBLIC    "../include/common"
                                                "../include/aligner"
                                                "${shaped_slp_SOURCE_DIR}" 
                                                "${FOLCA_SOURCE_DIR}" 
                                                "${SUX_SOURCE_DIR}/function" 
                                                "${SUX_SOURCE_DIR}/support"
                                                )
target_compile_options(extract PUBLIC "-std=c++17")