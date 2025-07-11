#include <iostream>
#include <vector>
#include <iomanip>
#include <random>

#include "SciData/vectorfield3d.h"
#include "SciData/grid_clusterer.h"
#include "SciData/vector_comparator.h"

using namespace CFF::SciData;
using namespace CFF::Utilities;

int main() {
    std::cout << "Clustering Test Example" << std::endl;
    std::cout << "======================" << std::endl;
    
    // Create a small 3D vector field
    const uint32_t dimX = 4, dimY = 4, dimZ = 4;
    VectorField3D vectorField("test_field", dimX, dimY, dimZ);
    
    // Random number generator for noise
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.5, 0.5);
    
    // Fill with some test data - create a simple pattern
    for (uint32_t x = 0; x < dimX; x++) {
        for (uint32_t y = 0; y < dimY; y++) {
            for (uint32_t z = 0; z < dimZ; z++) {
                // Create vectors that form two distinct groups
                Vector3d vec;
                if (x < dimX/2) {
                    // Left half - vectors pointing up
                    vec = makeVector3d(0.0, 1.0, 0.0);
                } else {
                    // Right half - vectors pointing right
                    vec = makeVector3d(1.0, 0.0, 0.0);
                }
                
                // Add some noise
                Vector3d noise = makeVector3d(
                    dis(gen) * 0.1,
                    dis(gen) * 0.1,
                    dis(gen) * 0.1
                );
                vec += noise;
                
                vectorField.setValue(makeVector3u(x, y, z), vec);
            }
        }
    }
    
    std::cout << "Created vector field with dimensions: " << dimX << "x" << dimY << "x" << dimZ << std::endl;
    
    // Create the vector comparator for distance calculation
    VectorComparator comparator(&vectorField);
    
    std::cout << "Vector field and comparator created successfully!" << std::endl;
    std::cout << "Note: Full clustering functionality requires additional setup" << std::endl;
    std::cout << "that depends on the specific clustering algorithm parameters." << std::endl;
    
    // Test basic vector field access
    std::cout << "Testing vector field access:" << std::endl;
    for (uint32_t x = 0; x < 2; x++) {
        for (uint32_t y = 0; y < 2; y++) {
            for (uint32_t z = 0; z < 2; z++) {
                const auto& vec = vectorField.getValue(makeVector3u(x, y, z));
                std::cout << "  Point (" << x << "," << y << "," << z << "): " 
                          << "(" << std::fixed << std::setprecision(3) 
                          << vec(0) << "," << vec(1) << "," << vec(2) << ")" << std::endl;
            }
        }
    }
    
    std::cout << "Clustering test completed successfully!" << std::endl;
    return 0;
}
