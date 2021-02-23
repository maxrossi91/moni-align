
add_executable(compress_dictionary compress_dictionary.cpp)
target_link_libraries(compress_dictionary common sdsl malloc_count)

FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(ssw)
FetchContent_GetProperties(ksw2)
FetchContent_GetProperties(klib)
FetchContent_GetProperties(bigbwt)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)


add_executable(ms matching_statistics.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(ms common sdsl divsufsort divsufsort64 malloc_count ri pthread)
target_include_directories(ms PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        "${bigbwt_SOURCE_DIR}"
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

add_executable(rlebwt_ms_full_build rlebwt_ms_full_build.cpp)
target_link_libraries(rlebwt_ms_full_build common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(rlebwt_ms_full_build PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(rlebwt_ms_full_build PUBLIC "-std=c++17")

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
target_link_libraries(align_klib_omp common malloc_count sdsl divsufsort divsufsort64 ri klib ssw gomp)
target_include_directories(align_klib_omp PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ssw_SOURCE_DIR}/src"
                                            "${klib_SOURCE_DIR}"
                                        )
target_compile_options(align_klib_omp PUBLIC "-std=c++17")
target_compile_options(align_klib_omp PUBLIC "-fopenmp")

add_executable(align_klib align_klib.cpp ${klib_SOURCE_DIR}/ksw.c ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(align_klib common malloc_count sdsl divsufsort divsufsort64 ri klib ssw pthread)
target_include_directories(align_klib PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ssw_SOURCE_DIR}/src"
                                            "${klib_SOURCE_DIR}"
                                            "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(align_klib PUBLIC "-std=c++17")

add_executable(align_full_klib align_full_klib.cpp ${klib_SOURCE_DIR}/ksw.c ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(align_full_klib common malloc_count sdsl divsufsort divsufsort64 ri klib ssw pthread)
target_include_directories(align_full_klib PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ssw_SOURCE_DIR}/src"
                                            "${klib_SOURCE_DIR}"
                                            "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(align_full_klib PUBLIC "-std=c++17")

add_executable(align_full_ksw2 align_full_ksw2.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(align_full_ksw2 common sdsl malloc_count divsufsort divsufsort64 ri ksw2 pthread)
target_include_directories(align_full_ksw2 PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ksw2_SOURCE_DIR}"
                                            "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(align_full_ksw2 PUBLIC "-std=c++17")
                                            
add_executable(ms_stats ms_stats.cpp)
target_link_libraries(ms_stats common sdsl divsufsort divsufsort64 malloc_count ri ssw)
target_include_directories(ms_stats PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(ms_stats PUBLIC "-std=c++17")
                                            
add_executable(ms_full_stats ms_full_stats.cpp)
target_link_libraries(ms_full_stats common sdsl divsufsort divsufsort64 malloc_count ri ssw)
target_include_directories(ms_full_stats PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(ms_full_stats PUBLIC "-std=c++17")

add_executable(ms_full_test ms_full_test.cpp)
target_link_libraries(ms_full_test common sdsl divsufsort divsufsort64 malloc_count ri ssw)
target_include_directories(ms_full_test PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(ms_full_test PUBLIC "-std=c++17")

add_executable(thr_stats thr_stats.cpp)
target_link_libraries(thr_stats common sdsl divsufsort divsufsort64 malloc_count ri ssw)
target_include_directories(thr_stats PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(thr_stats PUBLIC "-std=c++17")

add_executable(thr_test thr_test.cpp)
target_link_libraries(thr_test common sdsl divsufsort divsufsort64 malloc_count ri ssw)
target_include_directories(thr_test PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(thr_test PUBLIC "-std=c++17")
