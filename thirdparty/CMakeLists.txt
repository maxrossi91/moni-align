# Thirdparty libraries
include(FetchContent)
include(FetchContent_Install)

# Log directory for external content builds
make_directory(${CMAKE_CURRENT_BINARY_DIR}/logs)
set(INSTALL_FILE  ${CMAKE_CURRENT_BINARY_DIR}/install.cmake)

# Set environment

list(APPEND CMAKE_LIBRARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/lib ")
list(APPEND CMAKE_INCLUDE_PATH "${CMAKE_CURRENT_BINARY_DIR}/include ")
list(APPEND CMAKE_PREFIX_PATH  "${CMAKE_CURRENT_BINARY_DIR}")

# External Project resources
# ------------------------------------------------------------------------------

# DIVSUFSORT
find_library(DIVSUFSORT_LIB divsufsort HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
if(NOT DIVSUFSORT_LIB)
    message(STATUS "divsufsort library not found. Building as an external content")
else()
    message(STATUS "divsufsort library found at ${DIVSUFSORT_LIB}.")
endif()

# DIVSUFSORT64
find_library(DIVSUFSORT64_LIB divsufsort64 HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
if(NOT DIVSUFSORT64_LIB)
    message(STATUS "divsufsort64 library not found. Building as an external content")
else()
    message(STATUS "divsufsort64 library found at ${DIVSUFSORT64_LIB}.")
endif()

# SDSL
find_library(SDSL_LIB sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
find_path(SDSL_SRC sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/include PATHS ${CONDA_PREFIX}/include)
if(NOT SDSL_LIB)
    message(STATUS "sdsl library not found. Building as an external content")
    FetchContent_Declare(sdsl
            GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
            # GIT_TAG "v2.1.1"
    )

    if(NOT sdsl_POPULATED)
        FetchContent_Populate(sdsl)
        # Install in the current binary dir
        execute_process(COMMAND ./install.sh ${CMAKE_CURRENT_BINARY_DIR}
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${sdsl_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log)
        file(APPEND ${INSTALL_FILE} "MESSAGE(STATUS \"Installing sdsl.\")\n")
        file(APPEND ${INSTALL_FILE} "execute_process(COMMAND ./install.sh \${CMAKE_INSTALL_PREFIX}
                WORKING_DIRECTORY ${sdsl_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log)\n")
        if(result)
            message(FATAL_ERROR "CMake step for sdsl-lite failed: ${result}")
        endif()
    endif()

    set(SDSL_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(DIVSUFSORT_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(SDSL_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a)
    set(DIVSUFSORT_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort.a)
    set(DIVSUFSORT64_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort64.a)

else()
    message(STATUS "sdsl library found at ${SDSL_LIB}.")
    message(STATUS "sdsl sources found at ${SDSL_SRC}.")
endif()

add_library(sdsl STATIC IMPORTED GLOBAL)
set_property(TARGET sdsl PROPERTY IMPORTED_LOCATION ${SDSL_LIB})
target_include_directories(sdsl INTERFACE ${SDSL_SRC})

add_library(divsufsort STATIC IMPORTED GLOBAL)
set_property(TARGET divsufsort PROPERTY IMPORTED_LOCATION ${DIVSUFSORT_LIB})
    
add_library(divsufsort64 STATIC IMPORTED GLOBAL)
set_property(TARGET divsufsort64 PROPERTY IMPORTED_LOCATION ${DIVSUFSORT64_LIB})

# HTSLIB
find_library(HTS_LIB hts)
find_path(HTS_SRC htslib)
if(NOT HTS_LIB)
    message(STATUS "htslib library not found. Building as an external content.")
    FetchContent_Declare(htslib
            GIT_REPOSITORY "https://github.com/samtools/htslib"
            GIT_TAG "1.14"
    )

    if(NOT htslib_POPULATED)
        FetchContent_Populate(htslib)
        # Install in the current binary dir
        execute_process(COMMAND autoreconf -i
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${htslib_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-autoreconf.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-autoreconf.log)
                
        file(APPEND ${INSTALL_FILE} "MESSAGE(STATUS \"Installing htslib.\")\n")
        file(APPEND ${INSTALL_FILE} "execute_process(COMMAND autoreconf -i
                WORKING_DIRECTORY ${htslib_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-autoreconf.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-autoreconf.log)\n")
        if(result)
            message(FATAL_ERROR "CMake step for htslib failed: ${result}")
        endif()

        execute_process(COMMAND ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${htslib_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-configure.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-configure.log)
        file(APPEND ${INSTALL_FILE} "execute_process(COMMAND ./configure --prefix=\${CMAKE_INSTALL_PREFIX}
                WORKING_DIRECTORY ${htslib_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-configure.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-configure.log)\n")

        if(result)
            message(FATAL_ERROR "CMake step for htslib failed: ${result}")
        endif()

        execute_process(COMMAND make
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${htslib_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-make.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-make.log)
        file(APPEND ${INSTALL_FILE} "execute_process(COMMAND make
                WORKING_DIRECTORY ${htslib_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-make.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-make.log)\n")
        if(result)
            message(FATAL_ERROR "CMake step for htslib failed: ${result}")
        endif()

        execute_process(COMMAND make install
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${htslib_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-install.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-install.log)
        file(APPEND ${INSTALL_FILE} "execute_process(COMMAND make install
                WORKING_DIRECTORY ${htslib_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-install.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/htslib-install.log)\n")
        if(result)
            message(FATAL_ERROR "CMake step for htslib failed: ${result}")
        endif()
    endif()

    set(HTS_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(HTS_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a)

else()
    message(STATUS "htslib library found at ${HTS_LIB}.")
    message(STATUS "htslib sources found at ${HTS_SRC}.")
endif()

add_library(htslib STATIC IMPORTED GLOBAL)
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${HTS_LIB})
target_include_directories(htslib INTERFACE ${HTS_SRC})

# PFP++
find_program(PFP_CHECK check)
if(NOT PFP_CHECK)
    message(STATUS "pfp check not found. Building as an external content")
else()
    message(STATUS "pfp check found at ${PFP_CHECK}.")
endif()

find_program(PFP pfp++)
if(NOT PFP)
    message(STATUS "pfp not found. Building as an external content.")
    FetchContent_Declare(pfp
            GIT_REPOSITORY "https://github.com/maxrossi91/pfp"
            GIT_TAG leviosam_integration
    )

    if(NOT pfp_POPULATED)
        FetchContent_Populate(pfp)
        FetchContent_Install(   pfp  
                                SOURCE ${pfp_SOURCE_DIR} 
                                BINARY ${pfp_BINARY_DIR} 
                                DESTINATION ${CMAKE_CURRENT_BINARY_DIR} 
                                LOG ${CMAKE_CURRENT_BINARY_DIR}/logs
                                FILE ${INSTALL_FILE})

    endif()
    set(PFP ${CMAKE_CURRENT_BINARY_DIR}/bin/pfp++)
    set(PFP_CHECK ${CMAKE_CURRENT_BINARY_DIR}/bin/check)

else()
    message(STATUS "pfp found at ${PFP}.")
endif()


## BigRePair
find_program(largerepair largeb_irepair)
find_program(preprocess procdic)
find_program(postproc postproc)
if(NOT largerepair OR NOT preprocess OR NOT postproc)
    message(STATUS "BigRePair not found. Building as an external content.")
    FetchContent_Declare(
        bigrepair
        GIT_REPOSITORY https://gitlab.com/maxrossi91/bigrepair.git
    )
    
    FetchContent_GetProperties(bigrepair)
    if(NOT bigrepair_POPULATED)

        FetchContent_Populate(bigrepair)
        FetchContent_Install(   bigrepair  
                                SOURCE ${bigrepair_SOURCE_DIR} 
                                BINARY ${bigrepair_BINARY_DIR} 
                                DESTINATION ${CMAKE_CURRENT_BINARY_DIR} 
                                LOG ${CMAKE_CURRENT_BINARY_DIR}/logs
                                FILE ${INSTALL_FILE})

    endif()

    set(largerepair ${CMAKE_CURRENT_BINARY_DIR}/bin/largeb_irepair)
    set(preprocess ${CMAKE_CURRENT_BINARY_DIR}/bin/procdic)
    set(postproc ${CMAKE_CURRENT_BINARY_DIR}/bin/irepair)


else()
    message(STATUS "BigRePair components found at ${largerepair} ${preprocess} ${prostproc}.")
endif()

## pfp-thresholds
find_program(pfp_thresholds pfp_thresholds)
if(NOT pfp_thresholds)
    message(STATUS "pfp_thresholds not found. Building as an external content.")
    FetchContent_Declare(
        pfp_thresholds
        GIT_REPOSITORY https://github.com/maxrossi91/pfp-thresholds.git
        GIT_TAG reorganization
    )
    
    FetchContent_GetProperties(pfp_thresholds)
    if(NOT pfp_thresholds_POPULATED)

        FetchContent_Populate(pfp_thresholds)
        FetchContent_Install(   pfp_thresholds  
                                SOURCE ${pfp_thresholds_SOURCE_DIR} 
                                BINARY ${pfp_thresholds_BINARY_DIR} 
                                DESTINATION ${CMAKE_CURRENT_BINARY_DIR} 
                                LOG ${CMAKE_CURRENT_BINARY_DIR}/logs
                                FILE ${INSTALL_FILE})

    endif()

    set(pfp_thresholds ${CMAKE_CURRENT_BINARY_DIR}/bin/pfp_thresholds)

else()
    message(STATUS "pfp_thresholds found at ${pfp_thresholds}.")
endif()


install(SCRIPT ${INSTALL_FILE})


# In project libraries
# ------------------------------------------------------------------------------

## Add Big-BWT
FetchContent_Declare(
  bigbwt
  GIT_REPOSITORY https://github.com/alshai/Big-BWT.git
  )
  
FetchContent_GetProperties(bigbwt)
if(NOT bigbwt_POPULATED)
  FetchContent_Populate(bigbwt)
  add_subdirectory(${bigbwt_SOURCE_DIR} ${bigbwt_BINARY_DIR})
  
endif()

## Add malloc_count
FetchContent_Declare(
  malloc_count
  GIT_REPOSITORY https://github.com/bingmann/malloc_count
  )
  
FetchContent_GetProperties(malloc_count)
if(NOT malloc_count_POPULATED)
  FetchContent_Populate(malloc_count)

  add_library(malloc_count OBJECT ${malloc_count_SOURCE_DIR}/malloc_count.c ${malloc_count_SOURCE_DIR}/malloc_count.h)
  target_link_libraries(malloc_count dl)
  target_include_directories(malloc_count PUBLIC "${malloc_count_SOURCE_DIR}")

  add_library(memprofile OBJECT ${malloc_count_SOURCE_DIR}/memprofile.h)
  target_include_directories(memprofile PUBLIC "${malloc_count_SOURCE_DIR}")
endif()

  ## Add klib
FetchContent_Declare(
  klib
  GIT_REPOSITORY https://github.com/attractivechaos/klib
)

FetchContent_GetProperties(klib)
if(NOT klib_POPULATED)
  FetchContent_Populate(klib)
  
  add_library(klib INTERFACE)

  target_include_directories(klib INTERFACE ${klib_SOURCE_DIR})
endif()

## Add r-index
FetchContent_Declare(
  r-index
  GIT_REPOSITORY https://github.com/maxrossi91/r-index.git
)

FetchContent_GetProperties(r-index)
if(NOT r-index_POPULATED)
  FetchContent_Populate(r-index)

  add_subdirectory(${r-index_SOURCE_DIR} ${r-index_BINARY_DIR} )#EXCLUDE_FROM_ALL)
  add_library(ri INTERFACE)
  target_link_libraries(ri INTERFACE klib z)
  target_include_directories(ri INTERFACE ${r-index_SOURCE_DIR}/internal)
endif()

## Add ShapedSlp
FetchContent_Declare(
  shaped_slp
  GIT_REPOSITORY https://github.com/maxrossi91/ShapedSlp.git
  GIT_TAG master
  )
  
FetchContent_GetProperties(shaped_slp)
if(NOT shaped_slp_POPULATED)
  FetchContent_Populate(shaped_slp)
  add_subdirectory(${shaped_slp_SOURCE_DIR} ${shaped_slp_BINARY_DIR})
  set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
  set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)
endif()

## Add SSW
FetchContent_Declare(
    ssw
    GIT_REPOSITORY https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
    GIT_TAG master
  )
  
  FetchContent_GetProperties(ssw)
  if(NOT ssw_POPULATED)
    FetchContent_Populate(ssw)
    
    add_library(ssw OBJECT ${ssw_SOURCE_DIR}/src/ssw_cpp.cpp ${ssw_SOURCE_DIR}/src/ssw.c)
    target_include_directories(ssw PUBLIC ${ssw_SOURCE_DIR}/src)
  endif()

## Add Ksw2
FetchContent_Declare(
    ksw2
    GIT_REPOSITORY https://github.com/lh3/ksw2
    GIT_TAG master
  )
  
  FetchContent_GetProperties(ksw2)
  if(NOT ksw2_POPULATED)
    FetchContent_Populate(ksw2)
    
    add_library(ksw2 OBJECT ${ksw2_SOURCE_DIR}/kalloc.c 
                            ${ksw2_SOURCE_DIR}/ksw2_gg.c 
                            ${ksw2_SOURCE_DIR}/ksw2_gg2.c 
                            ${ksw2_SOURCE_DIR}/ksw2_gg2_sse.c 
                            ${ksw2_SOURCE_DIR}/ksw2_extz.c 
                            ${ksw2_SOURCE_DIR}/ksw2_extz2_sse.c
                            ${ksw2_SOURCE_DIR}/ksw2_extd.c 
                            ${ksw2_SOURCE_DIR}/ksw2_extd2_sse.c 
                            ${ksw2_SOURCE_DIR}/ksw2_extf2_sse.c 
                            ${ksw2_SOURCE_DIR}/ksw2_exts2_sse.c)
    target_include_directories(ksw2 PUBLIC ${ksw2_SOURCE_DIR}/src)
  endif()

  ## Add leviosam
  FetchContent_Declare(
    leviosam
    GIT_REPOSITORY https://github.com/maxrossi91/levioSAM.git
    GIT_TAG master
    )
    
  FetchContent_GetProperties(leviosam)
  if(NOT leviosam_POPULATED)
    FetchContent_Populate(leviosam)
    set(BUILD_TESTS OFF CACHE BOOL "Avoid building LevioSAM tests")

    add_subdirectory(${leviosam_SOURCE_DIR} ${leviosam_BINARY_DIR} EXCLUDE_FROM_ALL)
    target_include_directories(lvsam PUBLIC ${leviosam_SOURCE_DIR}/src ${CMAKE_INCLUDE_PATH})
    target_link_libraries(lvsam htslib)
  endif()
  

  ## Add Catch2
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.0.0-preview3
  )

  FetchContent_MakeAvailable(Catch2)