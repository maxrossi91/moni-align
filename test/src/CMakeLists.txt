
FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(ssw)
FetchContent_GetProperties(ksw2)
FetchContent_GetProperties(klib)
FetchContent_GetProperties(bigbwt)
FetchContent_GetProperties(leviosam)
FetchContent_MakeAvailable(Catch2)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)

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
                                            
# add_executable(align_full_klib align_full_klib.cpp ${klib_SOURCE_DIR}/ksw.c ${bigbwt_SOURCE_DIR}/xerrors.c)
# target_link_libraries(align_full_klib common malloc_count sdsl divsufsort divsufsort64 ri klib ssw pthread)
# target_include_directories(align_full_klib PUBLIC    "../../include/ms" 
#                                             "../../include/common" 
#                                             "${shaped_slp_SOURCE_DIR}" 
#                                             "${FOLCA_SOURCE_DIR}" 
#                                             "${SUX_SOURCE_DIR}/function" 
#                                             "${SUX_SOURCE_DIR}/support"
#                                             "${ssw_SOURCE_DIR}/src"
#                                             "${klib_SOURCE_DIR}"
#                                             "${bigbwt_SOURCE_DIR}"
#                                         )
# target_compile_options(align_full_klib PUBLIC "-std=c++17")

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

add_executable(lifting_test lifting_test.cpp)
target_link_libraries(lifting_test common sdsl divsufsort divsufsort64 malloc_count klib z lvsam Catch2::Catch2)
target_include_directories(lifting_test PUBLIC
                                        "../../include/common" 
                                        "../../include/aligner"
                                        )
target_compile_options(lifting_test PUBLIC "-std=c++17")

add_executable(liftidx_test liftidx_test.cpp)
target_link_libraries(liftidx_test common sdsl divsufsort divsufsort64 malloc_count klib z lvsam)
target_include_directories(liftidx_test PUBLIC
                                            "../../include/common" 
                                            "../../include/aligner"
                                        )
target_compile_options(liftidx_test PUBLIC "-std=c++17")

add_executable(r-index_test r-index_test.cpp)
target_link_libraries(r-index_test common sdsl divsufsort divsufsort64 malloc_count ri klib z lvsam Catch2::Catch2)
target_include_directories(r-index_test PUBLIC
                                        "../../include/common" 
                                        "../../include/aligner"
                                        "../../include/ms"
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(r-index_test PUBLIC "-std=c++17")


add_executable(ldx_slp_test ldx_slp_test.cpp)
target_link_libraries(ldx_slp_test common common sdsl lvsam klib z Catch2::Catch2)
target_include_directories(ldx_slp_test PUBLIC
                                        "../../include/common" 
                                        "../../include/aligner"
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(ldx_slp_test PUBLIC "-std=c++17")

add_executable(sigsegv_bt_handler_test sigsegv_bt_handler_test.cpp)
target_link_libraries(sigsegv_bt_handler_test common sdsl klib backward)
target_include_directories(sigsegv_bt_handler_test PUBLIC
                                        "../../include/common" 
                                        )
target_compile_options(sigsegv_bt_handler_test PUBLIC "-std=c++17")

add_executable(backward_test backward_test.cpp)
target_link_libraries(backward_test backward)
target_compile_options(backward_test PUBLIC "-std=c++17")


add_executable(rle_string_test rle_string_test.cpp)
target_link_libraries(rle_string_test common sdsl divsufsort divsufsort64 malloc_count ri klib z lvsam Catch2::Catch2)
target_include_directories(rle_string_test PUBLIC
                                        "../../include/common" 
                                        "../../include/ms"
                                        )
target_compile_options(rle_string_test PUBLIC "-std=c++17")


# list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
# include(CTest)
# include(Catch)
# catch_discover_tests(lifting_test)

