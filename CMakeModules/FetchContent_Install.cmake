# Install target using cmake standard pipeline
# Usage: FetchContent_Install(<prefix> SOURCE <directory> DESTINATION <directory> BINARY <directory> LOG <directory>)
function (FetchContent_Install target)
    set(options "")
    set(oneValueArgs SOURCE DESTINATION BINARY LOG FILE)
    set(multiValueArgs "")
    cmake_parse_arguments(MY_INSTALL "${options}" "${oneValueArgs}"
                            "${multiValueArgs}" ${ARGN} )

    # Install in the current binary dir
    execute_process(COMMAND cmake -DCMAKE_INSTALL_PREFIX:PATH=${MY_INSTALL_DESTINATION} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH} -DCMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH} -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH} ${MY_INSTALL_SOURCE}
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${MY_INSTALL_BINARY}
            OUTPUT_FILE ${MY_INSTALL_LOG}/${target}-cmake.log
            ERROR_FILE ${MY_INSTALL_LOG}/${target}-cmake.log)
    file(APPEND ${MY_INSTALL_FILE} "MESSAGE(STATUS \"Installing ${target}.\")\n")
    file(APPEND ${MY_INSTALL_FILE} "execute_process(COMMAND cmake -DCMAKE_INSTALL_PREFIX:PATH=\${CMAKE_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH} -DCMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH} -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH} ${MY_INSTALL_SOURCE} 
            WORKING_DIRECTORY ${MY_INSTALL_BINARY}
            OUTPUT_FILE ${MY_INSTALL_LOG}/${target}-cmake.log
            ERROR_FILE ${MY_INSTALL_LOG}/${target}-cmake.log)\n")
    if(result)
        message(FATAL_ERROR "CMake step for ${target} failed: ${result}")
    endif()

    execute_process(COMMAND make
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${MY_INSTALL_BINARY} 
            OUTPUT_FILE ${MY_INSTALL_LOG}/${target}-make.log
            ERROR_FILE ${MY_INSTALL_LOG}/${target}-make.log)

    file(APPEND ${MY_INSTALL_FILE} "execute_process(COMMAND make 
        WORKING_DIRECTORY ${MY_INSTALL_BINARY}
        OUTPUT_FILE ${MY_INSTALL_LOG}/${target}-make.log
        ERROR_FILE ${MY_INSTALL_LOG}/${target}-make.log)\n")
    if(result)
        message(FATAL_ERROR "CMake step for ${target} failed: ${result}")
    endif()

    execute_process(COMMAND make install
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${MY_INSTALL_BINARY} 
            OUTPUT_FILE ${MY_INSTALL_LOG}/${target}-install.log
            ERROR_FILE ${MY_INSTALL_LOG}/${target}-install.log)
    
    file(APPEND ${MY_INSTALL_FILE} "execute_process(COMMAND make install
        WORKING_DIRECTORY ${MY_INSTALL_BINARY}
        OUTPUT_FILE ${MY_INSTALL_LOG}/${target}-install.log
        ERROR_FILE ${MY_INSTALL_LOG}/${target}-install.log)\n")
    if(result)
        message(FATAL_ERROR "CMake step for ${target} failed: ${result}")
    endif()
endfunction()