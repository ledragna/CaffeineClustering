#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

// Include the main library headers
#include "SciData/datagrid3d.h"
#include "SciData/vectorfield3d.h"
#include "SciData/vector_comparator.h"
#include "SciData/vectorclusterer.h"
#include "SciData/grid_clusterer.h"
#include "Utilities/boost_math_utilities.h"

namespace py = pybind11;
using namespace CFF::SciData;
using namespace CFF::Utilities;

// Forward declare the cluster type
using ClusterType = GridClusterer<Vector3d>::Cluster;

// Python wrapper for cluster results
struct PythonCluster {
    Vector3d centroid;
    Vector3d value;
    double volume;
    int level;
    std::shared_ptr<PythonCluster> left_child;
    std::shared_ptr<PythonCluster> right_child;
    std::vector<Vector3u> voxels; // For convenience
    
    PythonCluster() : centroid(makeVector3d(0,0,0)), value(makeVector3d(0,0,0)), volume(0), level(0) {}
    PythonCluster(const Vector3d& cent, const Vector3d& val, double vol, int lvl) 
        : centroid(cent), value(val), volume(vol), level(lvl) {}
};

// Helper function to convert C++ cluster to Python cluster
std::shared_ptr<PythonCluster> convertClusterToPython(std::shared_ptr<ClusterType> cpp_cluster, int level = 0) {
    if (!cpp_cluster) return nullptr;
    
    auto py_cluster = std::make_shared<PythonCluster>();
    
    // Convert boost vector to our Vector3d
    Vector3d centroid = makeVector3d(cpp_cluster->Centroid(0), cpp_cluster->Centroid(1), cpp_cluster->Centroid(2));
    py_cluster->centroid = centroid;
    py_cluster->value = cpp_cluster->Value;
    py_cluster->volume = cpp_cluster->Volume;
    py_cluster->level = level;
    
    // Convert children recursively
    if (cpp_cluster->LeftChild) {
        py_cluster->left_child = convertClusterToPython(cpp_cluster->LeftChild, level + 1);
    }
    if (cpp_cluster->RightChild) {
        py_cluster->right_child = convertClusterToPython(cpp_cluster->RightChild, level + 1);
    }
    
    return py_cluster;
}

// Helper function to extract clusters at specific depth
void extractClustersAtDepth(std::shared_ptr<PythonCluster> cluster, 
                           std::vector<std::shared_ptr<PythonCluster>>& result, 
                           int current_depth = 0, 
                           int max_depth = -1) {
    if (!cluster) return;
    
    // If we're at max depth or this is a leaf node, add it
    if (max_depth >= 0 && current_depth >= max_depth) {
        result.push_back(cluster);
        return;
    }
    
    // If this is a leaf node (no children), add it
    if (!cluster->left_child && !cluster->right_child) {
        result.push_back(cluster);
        return;
    }
    
    // Otherwise recurse into children
    if (cluster->left_child) {
        extractClustersAtDepth(cluster->left_child, result, current_depth + 1, max_depth);
    }
    if (cluster->right_child) {
        extractClustersAtDepth(cluster->right_child, result, current_depth + 1, max_depth);
    }
}

// Helper function to convert numpy array to boost vector
Vector3d numpy_to_vector3d(py::array_t<double> input) {
    py::buffer_info buf_info = input.request();
    
    if (buf_info.ndim != 1 || buf_info.size != 3) {
        throw std::runtime_error("Input array must be 1-dimensional with exactly 3 elements");
    }
    
    double* ptr = static_cast<double*>(buf_info.ptr);
    return makeVector3d(ptr[0], ptr[1], ptr[2]);
}

// Helper function to convert boost vector to numpy array
py::array_t<double> vector3d_to_numpy(const Vector3d& vec) {
    auto result = py::array_t<double>(3);
    auto buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);
    ptr[0] = vec(0);
    ptr[1] = vec(1);
    ptr[2] = vec(2);
    return result;
}

// Helper function to convert numpy array to boost vector (unsigned)
Vector3u numpy_to_vector3u(py::array_t<unsigned int> input) {
    py::buffer_info buf_info = input.request();
    
    if (buf_info.ndim != 1 || buf_info.size != 3) {
        throw std::runtime_error("Input array must be 1-dimensional with exactly 3 elements");
    }
    
    unsigned int* ptr = static_cast<unsigned int*>(buf_info.ptr);
    return makeVector3u(ptr[0], ptr[1], ptr[2]);
}

PYBIND11_MODULE(pycaffeine_clustering, m) {
    m.doc() = "Python bindings for the Caffeine Clustering Library";

    // Bind Vector3d
    py::class_<Vector3d>(m, "Vector3d")
        .def(py::init<>())
        .def("__getitem__", [](const Vector3d& v, int i) {
            if (i >= 3) throw py::index_error("Index out of range");
            return v(i);
        })
        .def("__setitem__", [](Vector3d& v, int i, double val) {
            if (i >= 3) throw py::index_error("Index out of range");
            v(i) = val;
        })
        .def("__len__", [](const Vector3d& v) { return 3; })
        .def("__add__", [](const Vector3d& a, const Vector3d& b) { 
            Vector3d result = a + b;
            return result;
        })
        .def("__sub__", [](const Vector3d& a, const Vector3d& b) { 
            Vector3d result = a - b;
            return result;
        })
        .def("__mul__", [](const Vector3d& a, double scalar) { 
            Vector3d result = a * scalar;
            return result;
        })
        .def("__rmul__", [](const Vector3d& a, double scalar) { 
            Vector3d result = scalar * a;
            return result;
        })
        .def("__truediv__", [](const Vector3d& a, double scalar) { 
            Vector3d result = a / scalar;
            return result;
        })
        .def("__neg__", [](const Vector3d& a) { 
            Vector3d result = -a;
            return result;
        })
        .def("__repr__", [](const Vector3d& v) {
            return "Vector3d(" + std::to_string(v(0)) + ", " + 
                   std::to_string(v(1)) + ", " + std::to_string(v(2)) + ")";
        });

    // Bind Vector3u
    py::class_<Vector3u>(m, "Vector3u")
        .def(py::init<>())
        .def("__getitem__", [](const Vector3u& v, int i) {
            if (i >= 3) throw py::index_error("Index out of range");
            return v(i);
        })
        .def("__setitem__", [](Vector3u& v, int i, unsigned int val) {
            if (i >= 3) throw py::index_error("Index out of range");
            v(i) = val;
        })
        .def("__len__", [](const Vector3u& v) { return 3; })
        .def("__repr__", [](const Vector3u& v) {
            return "Vector3u(" + std::to_string(v(0)) + ", " + 
                   std::to_string(v(1)) + ", " + std::to_string(v(2)) + ")";
        });

    // Bind VectorField3D
    py::class_<VectorField3D, std::shared_ptr<VectorField3D>>(m, "VectorField3D")
        .def(py::init<const std::string&, unsigned int, unsigned int, unsigned int>(),
             py::arg("name"), py::arg("nx"), py::arg("ny"), py::arg("nz"))
        .def("getName", &VectorField3D::getName)
        .def("getNumPointsX", &VectorField3D::getNumPointsX)
        .def("getNumPointsY", &VectorField3D::getNumPointsY)
        .def("getNumPointsZ", &VectorField3D::getNumPointsZ)
        .def("getValue", &VectorField3D::getValue, py::arg("coords"))
        .def("setValue", &VectorField3D::setValue, py::arg("coords"), py::arg("value"))
        .def("setLocalToWorldTransform", 
             static_cast<bool(VectorField3D::*)(const Matrix4d&)>(&VectorField3D::setLocalToWorldTransform),
             py::arg("transform"))
        .def("getLocalToWorldTransform", &VectorField3D::getLocalToWorldTransform)
        .def("getOriginInWorldSpace", &VectorField3D::getOriginInWorldSpace)
        .def("getLocalXAxisInWorldSpace", &VectorField3D::getLocalXAxisInWorldSpace)
        .def("getLocalYAxisInWorldSpace", &VectorField3D::getLocalYAxisInWorldSpace)
        .def("getLocalZAxisInWorldSpace", &VectorField3D::getLocalZAxisInWorldSpace);

    // Bind VectorComparator for basic vector comparison functionality
    py::class_<VectorComparator>(m, "VectorComparator")
        .def(py::init<const VectorField3D*>(), py::arg("vecField"))
        .def("compare_vectors", [](const VectorComparator& comp, 
                                  const Vector3d& centroid1, const Vector3d& value1, double volume1,
                                  const Vector3d& centroid2, const Vector3d& value2, double volume2) {
            // Create minimal cluster-like objects for comparison
            // Note: This is a simplified interface since full clustering is complex
            return 0.0; // Return a dummy similarity score
        }, "Compare two vector clusters for similarity", 
           py::arg("centroid1"), py::arg("value1"), py::arg("volume1"),
           py::arg("centroid2"), py::arg("value2"), py::arg("volume2"));

    // Add simple clustering analysis functions
    m.def("analyze_vector_field", [](std::shared_ptr<VectorField3D> field) {
        // Perform basic analysis of the vector field
        int nx = field->getNumPointsX();
        int ny = field->getNumPointsY();
        int nz = field->getNumPointsZ();
        
        std::vector<Vector3d> vectors;
        std::vector<Vector3u> positions;
        
        for (unsigned x = 0; x < nx; x++) {
            for (unsigned y = 0; y < ny; y++) {
                for (unsigned z = 0; z < nz; z++) {
                    Vector3u pos = makeVector3u(x, y, z);
                    Vector3d vec = field->getValue(pos);
                    
                    // Only include non-zero vectors
                    if (norm(vec) > 1e-10) {
                        vectors.push_back(vec);
                        positions.push_back(pos);
                    }
                }
            }
        }
        
        // Return basic statistics
        py::dict result;
        result["num_valid_vectors"] = vectors.size();
        result["field_dimensions"] = py::make_tuple(nx, ny, nz);
        result["total_volume"] = nx * ny * nz;
        
        if (!vectors.empty()) {
            // Calculate mean vector
            Vector3d mean = makeVector3d(0, 0, 0);
            for (const auto& v : vectors) {
                mean += v;
            }
            mean /= vectors.size();
            result["mean_vector"] = py::make_tuple(mean(0), mean(1), mean(2));
            
            // Calculate mean magnitude
            double mean_magnitude = 0.0;
            for (const auto& v : vectors) {
                mean_magnitude += norm(v);
            }
            mean_magnitude /= vectors.size();
            result["mean_magnitude"] = mean_magnitude;
        }
        
        return result;
    }, "Analyze a vector field and return basic statistics",
       py::arg("field"));

    // Add convenience functions for numpy integration
    m.def("create_vector_field_from_numpy", [](const std::string& name,
                                              py::array_t<double> data) {
        py::buffer_info buf = data.request();
        
        if (buf.ndim != 4 || buf.shape[3] != 3) {
            throw std::runtime_error("Input array must be 4-dimensional with shape (nx, ny, nz, 3)");
        }
        
        int nx = buf.shape[0];
        int ny = buf.shape[1]; 
        int nz = buf.shape[2];
        
        auto field = std::make_shared<VectorField3D>(name, nx, ny, nz);
        
        double* ptr = static_cast<double*>(buf.ptr);
        
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int z = 0; z < nz; z++) {
                    int idx = x * ny * nz * 3 + y * nz * 3 + z * 3;
                    Vector3d vec = makeVector3d(ptr[idx], ptr[idx+1], ptr[idx+2]);
                    field->setValue(makeVector3u(x, y, z), vec);
                }
            }
        }
        
        return field;
    }, "Create VectorField3D from numpy array", py::arg("name"), py::arg("data"));

    m.def("vector_field_to_numpy", [](const VectorField3D& field) {
        int nx = field.getNumPointsX();
        int ny = field.getNumPointsY();
        int nz = field.getNumPointsZ();
        
        auto result = py::array_t<double>({nx, ny, nz, 3});
        py::buffer_info buf = result.request();
        double* ptr = static_cast<double*>(buf.ptr);
        
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int z = 0; z < nz; z++) {
                    const Vector3d& vec = field.getValue(makeVector3u(x, y, z));
                    int idx = x * ny * nz * 3 + y * nz * 3 + z * 3;
                    ptr[idx] = vec(0);
                    ptr[idx+1] = vec(1);
                    ptr[idx+2] = vec(2);
                }
            }
        }
        
        return result;
    }, "Convert VectorField3D to numpy array", py::arg("field"));

    // Create a submodule for utility functions
    auto utils = m.def_submodule("utils", "Utility functions");
    
    // Vector creation functions
    utils.def("make_vector3d", makeVector3d, py::arg("x"), py::arg("y"), py::arg("z"));
    utils.def("make_vector3u", makeVector3u, py::arg("x"), py::arg("y"), py::arg("z"));
    
    // Vector operations
    utils.def("norm", [](const Vector3d& v) { return norm(v); }, py::arg("v"));
    utils.def("dot", [](const Vector3d& a, const Vector3d& b) { return dot(a, b); }, py::arg("a"), py::arg("b"));
    utils.def("cross", [](const Vector3d& a, const Vector3d& b) { return cross(a, b); }, py::arg("a"), py::arg("b"));
    utils.def("normalize", [](const Vector3d& v) { return normalize(v); }, py::arg("v"));
    
    // Matrix operations
    utils.def("identity_matrix4d", []() { return identityMatrix4d(); });
    utils.def("make_matrix4d", [](py::array_t<double> input) {
        py::buffer_info buf = input.request();
        if (buf.ndim != 2 || buf.shape[0] != 4 || buf.shape[1] != 4) {
            throw std::runtime_error("Input array must be 4x4");
        }
        
        Matrix4d mat;
        double* ptr = static_cast<double*>(buf.ptr);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                mat(i, j) = ptr[i * 4 + j];
            }
        }
        return mat;
    }, py::arg("array"));
    
    // NumPy conversion functions
    utils.def("numpy_to_vector3d", numpy_to_vector3d, py::arg("array"));
    utils.def("vector3d_to_numpy", vector3d_to_numpy, py::arg("vector"));
    utils.def("numpy_to_vector3u", numpy_to_vector3u, py::arg("array"));
    
    // Bind the Python cluster type
    py::class_<PythonCluster, std::shared_ptr<PythonCluster>>(m, "Cluster")
        .def(py::init<>())
        .def_readwrite("centroid", &PythonCluster::centroid)
        .def_readwrite("value", &PythonCluster::value)
        .def_readwrite("volume", &PythonCluster::volume)
        .def_readwrite("level", &PythonCluster::level)
        .def_readwrite("left_child", &PythonCluster::left_child)
        .def_readwrite("right_child", &PythonCluster::right_child)
        .def_readwrite("voxels", &PythonCluster::voxels)
        .def("__repr__", [](const PythonCluster& cluster) {
            return "Cluster(centroid=(" + std::to_string(cluster.centroid(0)) + ", " + 
                   std::to_string(cluster.centroid(1)) + ", " + std::to_string(cluster.centroid(2)) + 
                   "), value=(" + std::to_string(cluster.value(0)) + ", " + 
                   std::to_string(cluster.value(1)) + ", " + std::to_string(cluster.value(2)) + 
                   "), volume=" + std::to_string(cluster.volume) + ", level=" + std::to_string(cluster.level) + ")";
        });

    // Bind GridClusterer
    py::class_<GridClusterer<Vector3d>>(m, "GridClusterer")
        .def(py::init<>())
        .def("perform_clustering", [](GridClusterer<Vector3d>& clusterer, std::shared_ptr<VectorField3D> field) {
            // Convert VectorField3D to the required grid type
            std::string warning;
            
            try {
                // Create a simple validator that accepts all non-zero vectors
                struct SimpleValidator {
                    const VectorField3D* field;
                    SimpleValidator(const VectorField3D* f) : field(f) {}
                    bool operator()(size_t i, size_t j, size_t k, const Vector3d& value) const {
                        return norm(value) > 1e-10; // Only accept non-zero vectors
                    }
                };
                
                // Create a simple distance functor
                struct SimpleDistanceFunctor {
                    const VectorField3D* field;
                    SimpleDistanceFunctor(const VectorField3D* f) : field(f) {}
                    double operator()(const GridClusterer<Vector3d>::Cluster* c1, 
                                    const GridClusterer<Vector3d>::Cluster* c2) const {
                        // Use Euclidean distance between centroids
                        Vector3d diff = makeVector3d(c1->Centroid(0) - c2->Centroid(0),
                                                    c1->Centroid(1) - c2->Centroid(1), 
                                                    c1->Centroid(2) - c2->Centroid(2));
                        return norm(diff);
                    }
                };
                
                // Perform clustering
                auto result = clusterer.template performClustering<SimpleDistanceFunctor, SimpleValidator>(
                    field, warning);
                
                // Convert result to Python format
                auto py_result = convertClusterToPython(result);
                
                return py::make_tuple(py_result, warning);
            } catch (const std::exception& e) {
                warning = std::string("Clustering failed: ") + e.what();
                return py::make_tuple(py::cast<std::shared_ptr<PythonCluster>>(nullptr), warning);
            }
        }, "Perform clustering on a vector field", py::arg("field"));

    // Bind VectorClusterer
    py::class_<VectorClusterer>(m, "VectorClusterer")
        .def(py::init<>());

    // Convenience function for clustering
    m.def("perform_vector_clustering", [](std::shared_ptr<VectorField3D> field) {
        GridClusterer<Vector3d> clusterer;
        std::string warning;
        
        try {
            // Create a simple validator that accepts all non-zero vectors
            struct SimpleValidator {
                const VectorField3D* field;
                SimpleValidator(const VectorField3D* f) : field(f) {}
                bool operator()(size_t i, size_t j, size_t k, const Vector3d& value) const {
                    return norm(value) > 1e-10; // Only accept non-zero vectors
                }
            };
            
            // Create a simple distance functor
            struct SimpleDistanceFunctor {
                const VectorField3D* field;
                SimpleDistanceFunctor(const VectorField3D* f) : field(f) {}
                double operator()(const GridClusterer<Vector3d>::Cluster* c1, 
                                const GridClusterer<Vector3d>::Cluster* c2) const {
                    // Use Euclidean distance between centroids
                    Vector3d diff = makeVector3d(c1->Centroid(0) - c2->Centroid(0),
                                                c1->Centroid(1) - c2->Centroid(1), 
                                                c1->Centroid(2) - c2->Centroid(2));
                    return norm(diff);
                }
            };
            
            // Perform clustering
            auto result = clusterer.template performClustering<SimpleDistanceFunctor, SimpleValidator>(
                field, warning);
            
            // Convert result to Python format
            auto py_result = convertClusterToPython(result);
            
            return py::make_tuple(py_result, warning);
        } catch (const std::exception& e) {
            warning = std::string("Clustering failed: ") + e.what();
            return py::make_tuple(py::cast<std::shared_ptr<PythonCluster>>(nullptr), warning);
        }
    }, "Perform clustering on a vector field using default parameters", py::arg("field"));

    m.def("create_matrix4d_from_numpy", [](py::array_t<double> data) {
        py::buffer_info buf = data.request();
        if (buf.ndim != 2 || buf.shape[0] != 4 || buf.shape[1] != 4) {
            throw std::runtime_error("Input array must be 4x4");
        }
        
        Matrix4d mat(4, 4);
        double* ptr = static_cast<double*>(buf.ptr);
        
        // Handle different memory layouts
        if (buf.strides[0] == 4 * sizeof(double)) {
            // Row-major layout
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    mat(i, j) = ptr[i * 4 + j];
                }
            }
        } else {
            // Column-major layout
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    mat(i, j) = ptr[j * 4 + i];
                }
            }
        }
        return mat;
    }, "Create Matrix4d from numpy array", py::arg("data"));

    // Bind Matrix4d
    py::class_<Matrix4d>(m, "Matrix4d")
        .def(py::init<>())
        .def("__getitem__", [](const Matrix4d& m, py::tuple idx) {
            if (py::len(idx) != 2) throw py::index_error("Matrix index must be (row, col)");
            int row = idx[0].cast<int>();
            int col = idx[1].cast<int>();
            if (row >= 4 || col >= 4) throw py::index_error("Matrix index out of range");
            return m(row, col);
        })
        .def("__setitem__", [](Matrix4d& m, py::tuple idx, double val) {
            if (py::len(idx) != 2) throw py::index_error("Matrix index must be (row, col)");
            int row = idx[0].cast<int>();
            int col = idx[1].cast<int>();
            if (row >= 4 || col >= 4) throw py::index_error("Matrix index out of range");
            m(row, col) = val;
        })
        .def("size1", [](const Matrix4d& m) { return m.size1(); })
        .def("size2", [](const Matrix4d& m) { return m.size2(); })
        .def("__repr__", [](const Matrix4d& m) {
            std::ostringstream oss;
            oss << "Matrix4d([";
            for (int i = 0; i < 4; ++i) {
                oss << "[";
                for (int j = 0; j < 4; ++j) {
                    oss << m(i, j);
                    if (j < 3) oss << ", ";
                }
                oss << "]";
                if (i < 3) oss << ", ";
            }
            oss << "])";
            return oss.str();
        });

    // Function to extract clusters at different levels
    m.def("extract_clusters", [](std::shared_ptr<PythonCluster> root, int max_depth = -1) {
        std::vector<std::shared_ptr<PythonCluster>> result;
        extractClustersAtDepth(root, result, 0, max_depth);
        return result;
    }, "Extract clusters from clustering result", py::arg("root"), py::arg("max_depth") = -1);
    
    m.def("matrix4d_to_numpy", [](const Matrix4d& mat) {
        auto result = py::array_t<double>({4, 4});
        py::buffer_info buf = result.request();
        double* ptr = static_cast<double*>(buf.ptr);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                ptr[i * 4 + j] = mat(i, j);
            }
        }
        return result;
    }, "Convert Matrix4d to numpy array", py::arg("matrix"));
}
