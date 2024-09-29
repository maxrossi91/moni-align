if (NOT TARGET common)
    add_subdirectory(common)
    target_include_directories(common INTERFACE common)
endif()