
FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(ssw)
FetchContent_GetProperties(ksw2)
FetchContent_GetProperties(klib)
FetchContent_GetProperties(bigbwt)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)


add_executable(build_spumoni build_spumoni.cpp)
target_link_libraries(build_spumoni common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(build_spumoni PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        )
target_compile_options(build_spumoni PUBLIC "-std=c++17")

add_executable(run_spumoni run_spumoni.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(run_spumoni common sdsl divsufsort divsufsort64 malloc_count ri pthread)
target_include_directories(run_spumoni PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(run_spumoni PUBLIC "-std=c++17")
