FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(ssw)
FetchContent_GetProperties(ksw2)
FetchContent_GetProperties(klib)
FetchContent_GetProperties(bigbwt)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)


add_executable(rlebwt_ms_full_build rlebwt_ms_full_build.cpp)
target_link_libraries(rlebwt_ms_full_build common sdsl divsufsort divsufsort64 malloc_count ri lvsam)
target_include_directories(rlebwt_ms_full_build PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(rlebwt_ms_full_build PUBLIC "-std=c++17")


add_executable(align_full_ksw2 align_full_ksw2.cpp ${bigbwt_SOURCE_DIR}/xerrors.c ${klib_SOURCE_DIR}/ksw.c )
target_link_libraries(align_full_ksw2 common sdsl malloc_count divsufsort divsufsort64 ri ksw2 klib pthread lvsam)
target_include_directories(align_full_ksw2 PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "../../include/aligner" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ksw2_SOURCE_DIR}"
                                            "${klib_SOURCE_DIR}"
                                            "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(align_full_ksw2 PUBLIC "-std=c++17")
