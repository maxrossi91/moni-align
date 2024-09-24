cmake_minimum_required(VERSION 3.15)
include_guard(GLOBAL)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

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
endif()
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY not defined. Setting variable to ${CMAKE_CURRENT_BINARY_DIR}/lib")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
endif()
if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  message(STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY not defined. Setting variable to ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif()

include(AddGitSubmodule)
include(FetchContent)
include(ExternalProject)

# SDSL
find_library(SDSL_LIB sdsl PATHS ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
if (NOT SDSL_LIB)
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
endif()

# HTSLIB
find_library(HTS_LIB hts PATHS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
if (NOT HTS_LIB)
  message(STATUS "htslib library not found. Building from submodule.")
  make_directory(${CMAKE_CURRENT_BINARY_DIR}/thirdparty/htslib)
  ExternalProject_Add(
    htslib_proj # Seemingly creates a target from the name which is why cannot name as htslib
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/htslib
    BINARY_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty/htslib
    UPDATE_COMMAND autoreconf -i
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND $(MAKE) install)
  
  set(HTS_SRC ${CMAKE_CURRENT_BINARY_DIR}/include CACHE PATH "Path to htslib headers" FORCE)
  set(HTS_LIB ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhts.so CACHE FILEPATH "Path to libhts.so" FORCE)
else()
  message(STATUS "htslib library found at ${HTS_LIB}.")
  message(STATUS "htslib sources found at ${HTS_SRC}.")
endif()

# This should work since the find_library command should only look in specified lib directory
add_library(htslib STATIC IMPORTED GLOBAL)
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${HTS_LIB})
target_include_directories(htslib INTERFACE ${HTS_SRC})    

# PFP
find_program(PFP pfp++ PATHS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
if (NOT PFP)
  message(STATUS "pfp++ not found. Building from submodule")
  add_git_submodule(${CMAKE_CURRENT_LIST_DIR}/thirdparty/pfp)
  set(PFP ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pfp++ CACHE PATH "Path to pfp++ executable" FORCE)
  set(PFP_CHECK ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/check CACHE PATH "Path to check executable" FORCE)
else()
  message(STATUS "pfp++ executable found at ${PFP}")
endif()

#add_subdirectory(thirdparty)


# Configure the compiler with the appropriate flags
# ------------------------------------------------------------------------------
# if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
#   # using Clang
#   include(ConfigureCompilerClang)
# elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#   # using GCC
#   include(ConfigureCompilerGcc)
# else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
# 	message(FATAL_ERROR "Only the compiler gcc and clang are supported")
# endif()


# add_subdirectory(include)
# add_subdirectory(src)
# add_subdirectory(test/src)
# add_subdirectory(utils)
# # add_subdirectory(benchmarks/src)

# set(PFP_EXE ${PFP})
# set(PFP_CHECK_EXE ${PFP_CHECK})
# set(PFP_THRESHOLDS_EXE ${pfp_thresholds})
# set(LARGEREPAIR_EXE ${largerepair})
# set(PREPROCESS_EXE ${preprocess})
# set(POSTPROCESS_EXE ${postproc})
# set(MONI_VERSION ${VERSION})

# # Configure pipeline for build folder
# set(USE_INSTALL_PATH False)
# configure_file(${PROJECT_SOURCE_DIR}/pipeline/moni.in ${PROJECT_BINARY_DIR}/moni @ONLY)

# # Configure pipeline for install folder
# set(USE_INSTALL_PATH True)
# configure_file(${PROJECT_SOURCE_DIR}/pipeline/moni.in ${PROJECT_BINARY_DIR}/moni.install @ONLY)


# install(TARGETS ms mems rlebwt_ms_build extend_ksw2 compress_dictionary build_seqidx TYPE RUNTIME)
# install(TARGETS rlebwt_ms_full_build align_full_ksw2 TYPE RUNTIME)
# install(TARGETS SlpEncBuild TYPE RUNTIME)
# # isntall(TARGETS pfp_thresholds pfp_thresholds64 TYPE RUNTIME)
# install(PROGRAMS ${PROJECT_BINARY_DIR}/moni.install RENAME moni TYPE BIN)

# # CPack
# # ------------------------------------------------------------------------------

# set(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
# set(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
# set(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
# set(CPACK_PACKAGE_VERSION "${VERSION}")

# include(InstallRequiredSystemLibraries)
# set(CPACK_GENERATOR "STGZ;TGZ;DEB")
# set(CPACK_SOURCE_GENERATOR "TGZ")
# set(CPACK_PACKAGE_VENDOR "University of Florida")
# set(CPACK_PACKAGE_CONTACT "rossi.m@ufl.edu")
# set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MONI-align - A Read Aligner with Multi-Genome References")
# set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
# set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
# set(CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")
# set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-sources")

# set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Massimiliano Rossi")
# set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT) 
# set(CPACK_COMPONENTS_GROUPING ALL_COMPONENTS_IN_ONE) # Groupp all components
# set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
# set(CPACK_DEB_COMPONENT_INSTALL YES)
# include(CPack)