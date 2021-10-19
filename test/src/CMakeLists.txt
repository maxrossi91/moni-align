
add_executable(compress_dictionary compress_dictionary.cpp)
target_link_libraries(compress_dictionary sdsl malloc_count)
target_include_directories(compress_dictionary PUBLIC "../../include/common")

FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(ssw)
FetchContent_GetProperties(ksw2)
FetchContent_GetProperties(klib)
FetchContent_GetProperties(bigbwt)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)

add_executable(run_spumoni run_spumoni.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(run_spumoni common sdsl divsufsort divsufsort64 malloc_count ri pthread)
target_include_directories(run_spumoni PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(run_spumoni PUBLIC "-std=c++17")

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

add_executable(build_spumoni build_spumoni.cpp)
target_link_libraries(build_spumoni common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(build_spumoni PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        )
target_compile_options(build_spumoni PUBLIC "-std=c++17")

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