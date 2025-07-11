# Caffeine Clustering Library

A standalone C++ library extracted from the Caffeine project originally developed by Andrea Salvadori, containing core DataGrid3D and clustering functionality. The library uses Boost.uBLAS as its mathematical backend.

## Features

- **DataGrid3D**: A 3D grid data structure for storing and manipulating volumetric data
- **VectorField3D**: Specialized 3D grid for vector field data  
- **GridClusterer**: Clustering algorithms for 3D grid data
- **VectorClusterer**: Vector-specific clustering implementations
- **VectorComparator**: Vector comparison utilities for clustering
- **Math Utilities**: Mathematical operations using Boost.uBLAS library


## Project Structure

```
caffeine-clustering/
├── include/
│   ├── SciData/
│   │   ├── datagrid3d.h          # Core 3D grid data structure
│   │   ├── vectorfield3d.h       # 3D vector field specialization
│   │   ├── grid_clusterer.h      # Clustering algorithms
│   │   ├── vectorclusterer.h     # Vector-specific clustering
│   │   └── vector_comparator.h   # Vector comparison utilities
│   └── Utilities/
│       ├── boost_math_utilities.h # Math utilities using Boost.uBLAS
│       ├── math3d.h              # 3D math operations
│       ├── math_backend.h        # Math backend definitions
│       ├── log.h                 # Logging utilities
│       ├── array3d.h             # 3D array template
│       ├── iterable_priority_queue.h # Priority queue implementation
│       ├── index3d_hashing.h     # 3D coordinate hashing
│       └── floating_point.h      # Floating point utilities
├── src/
│   ├── boost_math_utilities.cpp  # Math utility implementations
│   ├── math3d.cpp               # 3D math implementations
│   ├── index3d_hashing.cpp      # Hashing implementations
│   ├── vectorclusterer.cpp      # Vector clustering implementations
│   ├── vector_comparator.cpp    # Vector comparison implementations
│   └── floating_point.cpp      # Floating point utilities
├── examples/
│   ├── example_usage.cpp        # Basic usage examples
│   └── clustering_test.cpp      # Clustering test examples
├── CMakeLists.txt              # Build configuration
├── ClusteringConfig.cmake.in   # CMake package config
├── LICENSE.md                  # License information
└── README.md                   # This file
```

## Dependencies

- **CMake 3.10+**: Build system
- **C++17**: Language standard
- **Boost**: Boost.uBLAS library for mathematical operations
- **pkg-config**: For finding dependencies

## Building

1. **Install dependencies**:


2. **Clone and build**:
   ```bash
   git clone <repository>
   cd caffeine-clustering
   mkdir build && cd build
   cmake ..
   make
   ```

## Usage

The library provides a comprehensive C++ API for 3D data clustering and vector field analysis.

### Basic DataGrid3D Usage

```cpp
#include "SciData/vectorfield3d.h"
#include "SciData/vector_comparator.h"

using namespace CFF::SciData;
using namespace CFF::Utilities;

// Create a 3D vector field
VectorField3D field("my_field", 10, 10, 10);

// Set values using Boost vectors
Vector3d vector = makeVector3d(1.0, 0.0, 0.0);
field.setValue(makeVector3u(5, 5, 5), vector);

// Get values
auto retrieved = field.getValue(makeVector3u(5, 5, 5));

// Create comparator for analysis
VectorComparator comparator(&field);
```

## Build Options

- **Debug build**: `cmake -DCMAKE_BUILD_TYPE=Debug ..`
- **Release build**: `cmake -DCMAKE_BUILD_TYPE=Release ..`

## Installation

```bash
make install
```

This installs headers to `/usr/local/include` and libraries to `/usr/local/lib`.

## Integration

To use this library in your project:

```cmake
find_package(Clustering REQUIRED)
target_link_libraries(your_target clustering_lib)
```

## Status

This is a pure C++ library extracted from the Caffeine project. All Qt dependencies have been removed and replaced with Boost.uBLAS for mathematical operations.

## License

Extracted from the Caffeine project. Same of the original project.


