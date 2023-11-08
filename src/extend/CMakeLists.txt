FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(ssw)
FetchContent_GetProperties(ksw2)
FetchContent_GetProperties(klib)
FetchContent_GetProperties(bigbwt)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)

add_executable(extend_klib extend_klib.cpp ${klib_SOURCE_DIR}/ksw.c ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(extend_klib common malloc_count sdsl divsufsort divsufsort64 ri klib ssw pthread)
target_include_directories(extend_klib PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "../../include/extender" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ssw_SOURCE_DIR}/src"
                                            "${klib_SOURCE_DIR}"
                                            "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(extend_klib PUBLIC "-std=c++17")

add_executable(extend_ksw2 extend_ksw2.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(extend_ksw2 common sdsl malloc_count divsufsort divsufsort64 ri ksw2 pthread)
target_include_directories(extend_ksw2 PUBLIC    "../../include/ms" 
                                            "../../include/common"
                                            "${ksw2_SOURCE_DIR}"
                                            "../../include/extender" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(extend_ksw2 PUBLIC "-std=c++17")