# Caffeine Clustering Library

This is a standalone C++ library extracted from the Caffeine project, containing the core DataGrid3D and clustering functionality. The library has been migrated from mathfu to use Boost.uBLAS as its mathematical backend.

## Features

- **DataGrid3D**: A 3D grid data structure for storing and manipulating volumetric data
- **VectorField3D**: Specialized 3D grid for vector field data
- **GridClusterer**: Clustering algorithms for 3D grid data
- **VectorComparator**: Vector comparison utilities for clustering
- **Math Utilities**: Mathematical operations using Boost.uBLAS library
- **Python Bindings**: High-performance Python interface with NumPy integration

## Project Structure

```
standalone_clustering/
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
│       ├── log.h                 # Logging utilities
│       ├── array3d.h             # 3D array template
│       ├── iterable_priority_queue.h # Priority queue implementation
│       ├── index3d_hashing.h     # 3D coordinate hashing
│       └── floating_point.h      # Floating point utilities
├── src/
│   ├── boost_math_utilities.cpp
│   ├── math3d.cpp
│   ├── index3d_hashing.cpp
│   ├── vector_comparator.cpp
│   └── floating_point.cpp
├── examples/
│   ├── example_usage.cpp         # Basic usage example
│   └── clustering_test.cpp       # Vector field test
├── python/
│   └── bindings.cpp              # Python bindings with pybind11
├── CMakeLists.txt               # Build configuration
├── ClusteringConfig.cmake.in    # CMake package config
├── PYTHON_BINDINGS.md           # Python bindings documentation
└── README.md                    # This file
```

## Dependencies

- **CMake 3.10+**: Build system
- **C++17**: Language standard
- **Boost**: Boost.uBLAS library for mathematical operations
- **Qt5**: Core and Gui components (for vector conversion utilities)
- **pkg-config**: For finding dependencies

### Optional Dependencies (for Python bindings)

- **pybind11**: Python binding library
- **NumPy**: For Python array integration
- **Python 3.6+**: Python interpreter

## Building

1. **Install dependencies**:
   ```bash
   sudo apt update
   sudo apt install build-essential cmake pkg-config qtbase5-dev libboost-all-dev
   ```

2. **Clone and build**:
   ```bash
   git clone <repository>
   cd standalone_clustering
   mkdir build && cd build
   cmake ..
   make
   ```

3. **Run examples**:
   ```bash
   ./clustering_example      # Basic usage
   ./clustering_test         # Vector field test
   ```

### Python Bindings

The library includes optional Python bindings that provide high-performance access to the C++ functionality:

1. **Quick setup**:
   ```bash
   python3 setup_python_bindings.py --test
   ```

2. **Manual build**:
   ```bash
   mkdir build && cd build
   cmake -DBUILD_PYTHON_BINDINGS=ON ..
   make
   ```

3. **Usage**:
   ```python
   import sys
   sys.path.insert(0, 'build')
   import pycaffeine_clustering as cff
   
   # Create vectors
   v1 = cff.utils.make_vector3d(1.0, 2.0, 3.0)
   v2 = cff.utils.make_vector3d(4.0, 5.0, 6.0)
   
   # Mathematical operations
   norm = cff.utils.norm(v1)
   dot_product = cff.utils.dot(v1, v2)
   
   # Vector field operations
   field = cff.VectorField3D("field", 10, 10, 10)
   # ... work with NumPy arrays
   ```

See `PYTHON_BINDINGS.md` for detailed documentation and examples.

## Usage

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

### Advanced Features

The library provides sophisticated clustering algorithms and vector analysis tools. See the `examples/` directory for detailed usage patterns.

## Build Options

- **Debug build**: `cmake -DCMAKE_BUILD_TYPE=Debug ..`
- **Release build**: `cmake -DCMAKE_BUILD_TYPE=Release ..`
- **Python bindings**: `cmake -DBUILD_PYTHON_BINDINGS=ON ..`

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

✅ **Completed**:
- Core DataGrid3D functionality
- Vector field operations
- Mathematical utilities with Boost.uBLAS
- Basic clustering infrastructure
- CMake build system
- Example programs
- Complete migration from mathfu to Boost
- Python bindings with NumPy integration
- Comprehensive documentation

📝 **Potential Extensions**:
- Full clustering algorithm implementations
- Additional mathematical operations
- Performance optimizations
- Advanced Python features (matplotlib integration, etc.)

## License

Extracted from the Caffeine project. See original project for licensing terms.

## Contributing

This is a standalone extraction. For the full Caffeine project, see the main repository.
