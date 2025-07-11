cmake_minimum_required(VERSION 3.10)
project(CaffeineClustering)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(PkgConfig REQUIRED)

# Find Boost (required)
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost not found! Install with: sudo apt install libboost-all-dev")
endif()

# No Python bindings in this branch

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${Boost_INCLUDE_DIRS}
)

# Source files
set(SOURCES
    src/boost_math_utilities.cpp
    src/math3d.cpp
    src/index3d_hashing.cpp
    src/vector_comparator.cpp
    src/vectorclusterer.cpp
    src/floating_point.cpp
)

# Header files (for IDE organization)
set(HEADERS
    include/SciData/datagrid3d.h
    include/SciData/vectorfield3d.h
    include/SciData/grid_clusterer.h
    include/SciData/vectorclusterer.h
    include/SciData/vector_comparator.h
    include/Utilities/boost_math_utilities.h
    include/Utilities/math3d.h
    include/Utilities/log.h
    include/Utilities/array3d.h
    include/Utilities/iterable_priority_queue.h
    include/Utilities/index3d_hashing.h
)

# Create static library
add_library(clustering_lib STATIC ${SOURCES} ${HEADERS})
# No external library linking needed for clustering_lib

# No Python bindings in this branch

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
