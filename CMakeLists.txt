cmake_minimum_required(VERSION 3.10)
project(CaffeineClustering)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(PkgConfig REQUIRED)
find_package(Qt5 REQUIRED COMPONENTS Core Gui)

# Try to find mathfu
find_path(MATHFU_INCLUDE_DIR
    NAMES mathfu/vector.h
    PATHS 
        /usr/include
        /usr/local/include
        ${CMAKE_PREFIX_PATH}/include
        ${MATHFU_ROOT}/include
    DOC "Path to mathfu headers"
)

if(NOT MATHFU_INCLUDE_DIR)
    message(STATUS "mathfu not found, downloading from GitHub...")
    include(FetchContent)
    FetchContent_Declare(
        mathfu
        GIT_REPOSITORY https://github.com/google/mathfu.git
        GIT_TAG master
    )
    FetchContent_GetProperties(mathfu)
    if(NOT mathfu_POPULATED)
        FetchContent_Populate(mathfu)
        # Don't add mathfu subdirectory to avoid CMake compatibility issues
        set(MATHFU_INCLUDE_DIR ${mathfu_SOURCE_DIR}/include)
    endif()
endif()

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${MATHFU_INCLUDE_DIR}
)

# Add preprocessor definitions to disable SIMD in mathfu
add_definitions(-DMATHFU_COMPILE_WITHOUT_SIMD_SUPPORT)

# Source files
set(SOURCES
    src/mathfu_utilities.cpp
    src/math3d.cpp
    src/index3d_hashing.cpp
    src/vector_comparator.cpp
    src/floating_point.cpp
)

# Header files (for IDE organization)
set(HEADERS
    include/SciData/datagrid3d.h
    include/SciData/vectorfield3d.h
    include/SciData/grid_clusterer.h
    include/SciData/vectorclusterer.h
    include/SciData/vector_comparator.h
    include/Utilities/mathfu_utilities.h
    include/Utilities/math3d.h
    include/Utilities/log.h
    include/Utilities/array3d.h
    include/Utilities/iterable_priority_queue.h
    include/Utilities/index3d_hashing.h
)

# Create static library
add_library(clustering_lib STATIC ${SOURCES} ${HEADERS})
target_link_libraries(clustering_lib Qt5::Core Qt5::Gui)

# Example executable
add_executable(clustering_example examples/example_usage.cpp)
target_link_libraries(clustering_example clustering_lib)

# Clustering test executable
add_executable(clustering_test examples/clustering_test.cpp)
target_link_libraries(clustering_test clustering_lib)

# Enable debug symbols
set(CMAKE_BUILD_TYPE Debug)

# Compiler flags
target_compile_options(clustering_lib PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
)

# Installation
install(TARGETS clustering_lib 
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

# Export package information
set(CLUSTERING_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include)
set(CLUSTERING_LIBRARIES clustering_lib)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/ClusteringConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/ClusteringConfig.cmake"
    @ONLY
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ClusteringConfig.cmake"
    DESTINATION lib/cmake/Clustering
)
