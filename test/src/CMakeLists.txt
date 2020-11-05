
add_executable(compress_dictionary compress_dictionary.cpp)
target_link_libraries(compress_dictionary common sdsl malloc_count)

FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(ssw)
FetchContent_GetProperties(ksw2)
FetchContent_GetProperties(klib)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)


add_executable(ms matching_statistics.cpp)
target_link_libraries(ms common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(ms PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(ms PUBLIC "-std=c++17")

add_executable(rlebwt_ms_build rlebwt_ms_build.cpp)
target_link_libraries(rlebwt_ms_build common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(rlebwt_ms_build PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(rlebwt_ms_build PUBLIC "-std=c++17")

add_executable(shapedslp_test shapedslp_test.cpp)
target_link_libraries(shapedslp_test common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(shapedslp_test PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(shapedslp_test PUBLIC "-std=c++17")

add_executable(align align.cpp)
target_link_libraries(align common sdsl divsufsort divsufsort64 malloc_count ri ssw)
target_include_directories(align PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ssw_SOURCE_DIR}/src"
                                        )
target_compile_options(align PUBLIC "-std=c++17")


add_executable(align_omp align_omp.cpp)
target_link_libraries(align_omp common sdsl divsufsort divsufsort64 malloc_count ri ssw gomp)
target_include_directories(align_omp PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ssw_SOURCE_DIR}/src"
                                        )
target_link_directories(align_omp PUBLIC /s/apu/a/homes/rossi.m/miniconda3/envs/moni/lib)
target_compile_options(align_omp PUBLIC "-std=c++17")
target_compile_options(align_omp PUBLIC "-fopenmp")

add_executable(align_omp_mems align_omp_mems.cpp)
target_link_libraries(align_omp_mems common sdsl divsufsort divsufsort64 malloc_count ri ssw gomp)
target_include_directories(align_omp_mems PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ssw_SOURCE_DIR}/src"
                                        )
target_link_directories(align_omp_mems PUBLIC /s/apu/a/homes/rossi.m/miniconda3/envs/moni/lib)
target_compile_options(align_omp_mems PUBLIC "-std=c++17")
target_compile_options(align_omp_mems PUBLIC "-fopenmp")

add_executable(align_ksw2_omp align_ksw2_omp.cpp)
target_link_libraries(align_ksw2_omp common sdsl malloc_count divsufsort divsufsort64 ri ksw2 gomp)
target_include_directories(align_ksw2_omp PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ksw2_SOURCE_DIR}"
                                        )
target_compile_options(align_ksw2_omp PUBLIC "-std=c++17")
target_compile_options(align_ksw2_omp PUBLIC "-fopenmp")

add_executable(align_klib_omp align_klib_omp.cpp ${klib_SOURCE_DIR}/ksw.c)
target_link_libraries(align_klib_omp common malloc_count sdsl divsufsort divsufsort64 ri klib gomp)
target_include_directories(align_klib_omp PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ssw_SOURCE_DIR}"
                                            "${klib_SOURCE_DIR}"
                                        )
target_compile_options(align_klib_omp PUBLIC "-std=c++17")
target_compile_options(align_klib_omp PUBLIC "-fopenmp")
