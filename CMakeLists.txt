cmake_minimum_required(VERSION 3.15)
include_guard(GLOBAL)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/CMakeModules")

# Set a default build type if none was specified
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING 
      "Choose the type of build, options are: Release|Debug|RelWithDebInfo (for distros)." FORCE)
endif()

message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")

# About this project
# ------------------------------------------------------------------------------
project(moni-align)
set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "0")
set(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Set environment
# ------------------------------------------------------------------------------
find_package(Git REQUIRED)

# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include") # This is an hack because include(GUIInstallDirs) doesn't work

# Set the output directories for executables and libraries
if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY not defined. Setting variable to ${CMAKE_CURRENT_BINARY_DIR}/bin")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
else()
  message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY user defined as ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
endif()
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY not defined. Setting variable to ${CMAKE_CURRENT_BINARY_DIR}/lib")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
else()
  message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY user defined as ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
endif()
if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  message(STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY not defined. Setting variable to ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
else()
  message(STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY user defined as ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
endif()

include(AddGitSubmodule)
include(FetchContent)
include(ExternalProject)

# Set up CMAKE variables that control where to look 
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_BINARY_DIR}") #find_package
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}") #find_library

# ZLIB
find_package(ZLIB REQUIRED)

# SDSL
find_library(SDSL_LIB sdsl PATHS ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} NO_DEFAULT_PATH) # I am almost positive it is important that do not use local sdsl or htslib
find_library(DIVSUFSORT_LIB divsufsort PATHS ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} NO_DEFAULT_PATH)
find_library(DIVSUFSORT64_LIB divsufsort64 PATHS ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} NO_DEFAULT_PATH)
if (NOT SDSL_LIB OR NOT DIVSUFSORT_LIB OR NOT DIVSUFSORT64_LIB)
  message(STATUS "sdsl not found. Building from submodule.")
  add_git_submodule(${CMAKE_CURRENT_LIST_DIR}/thirdparty/sdsl-lite)
  set(SDSL_SRC ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/sdsl-lite/include/ CACHE PATH "Path to sdsl headers" FORCE)
  set(DIVSUFSORT_SRC ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/sdsl-lite/external/libdivsufsort/include CACHE PATH "Path to libdivsufsort" FORCE)
  set(SDSL_LIB ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/libsdsl.a CACHE PATH "Path to libsdsl.a" FORCE)
  set(DIVSUFSORT_LIB ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/libdivsufsort.a CACHE PATH "Path to libdivsufsort.a" FORCE)
  set(DIVSUFSORT64_LIB ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/libdivsufsort64.a CACHE PATH "Path to libdivsufsort64.a" FORCE)
else()
  message(STATUS "sdsl library found at ${SDSL_LIB}.")
  message(STATUS "sdsl sources found at ${SDSL_SRC}.")
  message(STATUS "divsufsort library found at ${DIVSUFSORT_LIB}.")
  message(STATUS "divsufsort64 library found at ${DIVSUFSORT64_LIB}.")
  message(STATUS "divsufsort sources found at ${DIVSUFSORT_SRC}.")
endif()

if (TARGET sdsl AND TARGET divsufsort AND TARGET divsufsort64)
  # Need to add this for PFP-thresholds
  target_include_directories(divsufsort INTERFACE ${DIVSUFSORT_SRC})
  target_include_directories(divsufsort64 INTERFACE ${DIVSUFSORT_SRC})
else()
  message(STATUS "Creating sdsl global library accessible to all CMakeLists.txt")
  add_library(sdsl STATIC IMPORTED GLOBAL)
  set_property(TARGET sdsl PROPERTY IMPORTED_LOCATION ${SDSL_LIB})
  target_include_directories(sdsl INTERFACE ${SDSL_SRC})

  message(STATUS "Creating divsufsort global library accessible to all CMakeLists.txt")
  add_library(divsufsort STATIC IMPORTED GLOBAL)
  set_property(TARGET divsufsort PROPERTY IMPORTED_LOCATION ${DIVSUFSORT_LIB})
  target_include_directories(divsufsort INTERFACE ${DIVSUFSORT_SRC})

  message(STATUS "Creating divsufsort64 global library accessible to all CMakeLists.txt")
  add_library(divsufsort64 STATIC IMPORTED GLOBAL)
  set_property(TARGET divsufsort64 PROPERTY IMPORTED_LOCATION ${DIVSUFSORT64_LIB})
  target_include_directories(divsufsort64 INTERFACE ${DIVSUFSORT_SRC})
endif()

# HTSLIB
get_filename_component(CMAKE_RUNTIME_OUTPUT_DIRECTORY_PARENT ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} DIRECTORY) # Get directory above bin directory
find_library(HTS_LIB hts PATHS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} NO_DEFAULT_PATH)
if (NOT HTS_LIB)
  message(STATUS "htslib library not found. Building from submodule.")
  make_directory(${CMAKE_CURRENT_BINARY_DIR}/thirdparty/htslib)
  ExternalProject_Add(
    htslib_proj # Seemingly creates a target from the name which is why cannot name as htslib
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/htslib
    BINARY_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/htslib
    UPDATE_COMMAND autoreconf -i
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_RUNTIME_OUTPUT_DIRECTORY_PARENT} # Might have to modify
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND $(MAKE) install)
  
  set(HTS_SRC ${CMAKE_CURRENT_BINARY_DIR}/include CACHE PATH "Path to htslib headers" FORCE)
  set(HTS_LIB ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhts.so CACHE FILEPATH "Path to libhts.so" FORCE)
else()
  message(STATUS "htslib library found at ${HTS_LIB}.")
  message(STATUS "htslib sources found at ${HTS_SRC}.")
endif()

# This should work since the find_library command should only look in specified lib directory
if (NOT TARGET htslib)
  message(STATUS "Creating htslib global library accessible to all CMakeLists.txt")
  add_library(htslib STATIC IMPORTED GLOBAL)
  set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${HTS_LIB})
  target_include_directories(htslib INTERFACE ${HTS_SRC})
endif() 

# Create FindHTSlib.cmake for PFP
configure_file(${CMAKE_CURRENT_LIST_DIR}/pipeline/FindHTSlib.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/FindHTSlib.cmake @ONLY)

# PFP
list(APPEND CMAKE_INCLUDE_PATH "${SDSL_SRC}" "${HTS_SRC}/htslib")
find_program(PFP pfp++ PATHS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
if (NOT PFP)
  message(STATUS "pfp++ not found. Building from submodule")
  add_git_submodule(${CMAKE_CURRENT_LIST_DIR}/thirdparty/pfp)
  set(PFP ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pfp++ CACHE PATH "Path to pfp++ executable" FORCE)
  set(PFP_CHECK ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/check CACHE PATH "Path to check executable" FORCE)
else()
  message(STATUS "pfp++ executable found at ${PFP}")
endif()

# BigRePair
find_program(largerepair largeb_irepair PATHS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
find_program(preprocess procdic PATHS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
find_program(postproc postproc PATHS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
if(NOT largerepair OR NOT preprocess OR NOT postproc)
  message(STATUS "bigrepair not found. Building from submodule")
  set(CMAKE_INSTALL_PREFIX ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  ExternalProject_Add(
    bigrepair_proj 
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/bigrepair
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/bigrepair
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND "") # make install DESTDIR=${} installs the executables in ${} but ${}/usr/local/bin
                        # make install prefix=${} doesn't seem to work

  # Copy the executables to the bin directory 
  add_custom_command(
    TARGET bigrepair_proj
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/bigrepair/procdic ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/bigrepair/postproc ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/bigrepair/largeb_repair/largeb_irepair ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
    DEPENDS bigrepair_proj  # Ensure this runs after the project is built
  )
  set(largerepair ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/largeb_irepair CACHE PATH "Path to largeb_irepair" FORCE)
  set(preprocess ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/procdic CACHE PATH "Path to procdic" FORCE)
  set(postproc ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/postproc CACHE PATH "Path to postproc" FORCE)

else()
  message(STATUS "BigRePair components found at ${largerepair} ${preprocess} ${postproc}.")
endif()

# PFP-Thresholds
find_program(pfp_thresholds pfp_thresholds PATHS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
if (NOT pfp_thresholds)
  message(STATUS "pfp_thresholds not found. Building as submodule.")
  add_git_submodule(${CMAKE_CURRENT_LIST_DIR}/thirdparty/pfp-thresholds)
  set(pfp_thresholds ${CMAKE_CURRENT_BINARY_DIR}/bin/pfp_thresholds CACHE PATH "Path to pfp_thresholds" FORCE)
else()
  message(STATUS "pfp_thresholds found at ${pfp_thresholds}.")
endif()

# In project libraries
#--------------------------------------------------------------------------
# BigBWT
if (NOT bigbwt_POPULATED)
  message(STATUS "bigbwt not found. Building as submodule.")
  add_git_submodule(${CMAKE_CURRENT_LIST_DIR}/thirdparty/Big-BWT)
  set(bigbwt_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/Big-BWT CACHE PATH "Path to BigBWT root directory")
  set(bigbwt_POPULATED TRUE CACHE BOOL "Set bigbwt_POPULATED to be TRUE")
else()
  message(STATUS "bigbwt found. Do not build from submodule.")
endif()

# Malloc Count
if (NOT TARGET malloc_count)
  message(STATUS "malloc_count not found. Creating library.")
  set(MALLOC_COUNT_SRC ${CMAKE_CURRENT_LIST_DIR}/thirdparty/malloc_count)
  add_library(malloc_count OBJECT ${MALLOC_COUNT_SRC}/malloc_count.c ${MALLOC_COUNT_SRC}/malloc_count.h)
  target_link_libraries(malloc_count dl)
  target_include_directories(malloc_count PUBLIC "${MALLOC_COUNT_SRC}")
else()
  message(STATUS "malloc_count found. Do not create the library.")
endif()

if (NOT TARGET memprofile)
  message(STATUS "memprofile not found. Creating library")
  add_library(memprofile OBJECT ${MALLOC_COUNT_SRC}/memprofile.h)
  target_include_directories(memprofile PUBLIC "${MALLOC_COUNT_SRC}")
else()
  message(STATUS "memprofile found. Do not create the library.")
endif()

# Klib
if (NOT TARGET klib)
  message(STATUS "klib not found. Creating library.")
  set(klib_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/klib)
  add_library(klib INTERFACE)
  target_include_directories(klib INTERFACE ${klib_SOURCE_DIR})
else()
  message(STATUS "klib found. Do not create the library.")
endif()

# r-index
if (NOT r-index_POPULATED)
  message(STATUS "r-index not found. Building as submodule.")
  add_git_submodule(${CMAKE_CURRENT_LIST_DIR}/thirdparty/r-index)
  set(r-index_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/r-index CACHE PATH "Path to r-index root directory")
  add_library(ri INTERFACE)
  target_link_libraries(ri INTERFACE klib z)
  target_include_directories(ri INTERFACE ${r-index_SOURCE_DIR}/internal)
  set(r-index_POPULATED TRUE CACHE BOOL "Set r-index_POPULATED to be TRUE")
else()
  message(STATUS "r-index found. Do not build from submodule.")
endif()

# ShapedSLP
if (NOT shaped_slp_POPULATED)
  message(STATUS "ShapedSLP not found. Building as submodule.")
  add_git_submodule(${CMAKE_CURRENT_LIST_DIR}/thirdparty/ShapedSlp)
  set(shaped_slp_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/ShapedSlp CACHE PATH "Path to ShapedSLP root directory" FORCE)
  set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
  set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)
  set(shaped_slp_POPULATED TRUE CACHE BOOL "Set shaped_slp_POPULATED to be TRUE")
else()
  message(STATUS "ShapedSLP found. Do not build as submodule.")
endif()

# SSW
if (NOT ssw_POPULATED)
  message(STATUS "Complete Striped Smith Waterman library not found. Creating library")
  set(ssw_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/Complete-Striped-Smith-Waterman-Library CACHE PATH "Path to ssw root directory" FORCE)
  add_library(ssw OBJECT ${ssw_SOURCE_DIR}/src/ssw_cpp.cpp ${ssw_SOURCE_DIR}/src/ssw.c)
  target_include_directories(ssw PUBLIC ${ssw_SOURCE_DIR}/src)
  set(ssw_POPULATED TRUE CACHE BOOL "Set ssw_POPULATED to be TRUE")
else()
  message(STATUS "Complete Striped Smith Waterman found. Do not create library.")
endif()

# KSW2
if (NOT ksw2_POPULATED)
  message(STATUS "ksw2 library not found. Creating library.")
  set(ksw2_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/ksw2 CACHE PATH "Path to ksw2 src directory")
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
  set(ksw2_POPULATED TRUE CACHE BOOL "Set ksw2_POPULATED to be TRUE")
else()
  message(STATUS "ksw2 library found. Do not create the library.")
endif()

# LevioSAM
if (NOT TARGET lvsam)
  message(STATUS "leviosam library not found. Building as submodule.")
  add_git_submodule(${CMAKE_CURRENT_LIST_DIR}/thirdparty/levioSAM)
  set(levioSAM_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/levioSAM/src CACHE PATH "Path to leviosam src directory")
  target_include_directories(lvsam PUBLIC ${leviosam_SOURCE_DIR} ${CMAKE_INCLUDE_PATH})
  target_link_libraries(lvsam htslib)
else()
  message(STATUS "leviosam library found. Do not create library")
endif()

# backward
if (NOT backward_POPULATED)
  message(STATUS "backward library not found. Building as submodule.")
  set(STACK_DETAILS_AUTO_DETECT FALSE CACHE BOOL "Use default stack trace since cannot guarantee that bfd dependencies exist on system" FORCE)
  set(BACKWARD_TESTS OFF CACHE BOOL "Avoid building backward-cpp tests")
  add_git_submodule(${CMAKE_CURRENT_LIST_DIR}/thirdparty/backward-cpp)
  set(backward_POPULATED TRUE CACHE BOOL "Set backward_POPULATED to be TRUE")
  target_link_libraries(backward INTERFACE ${BACKWARD_LIBRARIES})
else()
   message(STATUS "backward library found. Do not build from submodule.")
endif()

# Catch2
if (NOT catch2_POPULATED)
  message(STATUS "catch2 not found. Building as submodule.")
  add_git_submodule(${CMAKE_CURRENT_LIST_DIR}/thirdparty/Catch2)
  set(catch2_POPULATED TRUE CACHE BOOL "Set catch2_POPULATED to be TRUE")
else()
  message(STATUS "catch2 found. Do not build from submodule.")
endif()


#add_subdirectory(thirdparty)

# Configure the compiler with the appropriate flags
# ------------------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  include(ConfigureCompilerClang)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  include(ConfigureCompilerGcc)
else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message(FATAL_ERROR "Only the compiler gcc and clang are supported")
endif()

message(STATUS "Adding include CMakeLists")
add_subdirectory(include)
message(STATUS "Adding src CMakeLists")
add_subdirectory(src)
# message(STATUS "Adding test CMakeLists")
# add_subdirectory(test/src)
message(STATUS "Adding utils CMakeLists")
add_subdirectory(utils)

set(PFP_EXE ${PFP})
set(PFP_CHECK_EXE ${PFP_CHECK})
set(PFP_THRESHOLDS_EXE ${pfp_thresholds})
set(LARGEREPAIR_EXE ${largerepair})
set(PREPROCESS_EXE ${preprocess})
set(POSTPROCESS_EXE ${postproc})
set(MONI_VERSION ${VERSION})

# Configure pipeline for build folder
set(USE_INSTALL_PATH False)
configure_file(${PROJECT_SOURCE_DIR}/pipeline/moni.in ${PROJECT_BINARY_DIR}/moni @ONLY)

# Configure pipeline for install folder
set(USE_INSTALL_PATH True)
configure_file(${PROJECT_SOURCE_DIR}/pipeline/moni.in ${PROJECT_BINARY_DIR}/moni.install @ONLY)


install(TARGETS ms mems rlebwt_ms_build extend_ksw2 compress_dictionary build_seqidx TYPE RUNTIME)
install(TARGETS rlebwt_ms_full_build align_full_ksw2 TYPE RUNTIME)
install(TARGETS SlpEncBuild TYPE RUNTIME)
# install(TARGETS pfp_thresholds pfp_thresholds64 TYPE RUNTIME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/moni.install RENAME moni TYPE BIN)

# CPack
# ------------------------------------------------------------------------------

set(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set(CPACK_PACKAGE_VERSION "${VERSION}")

include(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "STGZ;TGZ;DEB")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "University of Florida")
set(CPACK_PACKAGE_CONTACT "rossi.m@ufl.edu")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MONI-align - A Read Aligner with Multi-Genome References")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-sources")

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Massimiliano Rossi")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT) 
set(CPACK_COMPONENTS_GROUPING ALL_COMPONENTS_IN_ONE) # Group all components
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEB_COMPONENT_INSTALL YES)
include(CPack)