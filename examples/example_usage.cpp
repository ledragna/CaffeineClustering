#include <iostream>
#include <memory>
#include <SciData/datagrid3d.h>
#include <SciData/vectorfield3d.h>
#include <SciData/grid_clusterer.h>
#include <SciData/vector_comparator.h>

using namespace CFF::SciData;

int main()
{
    std::cout << "Standalone Clustering Example" << std::endl;
    std::cout << "=============================" << std::endl;
    
    try {
        // Create a simple 3D vector field
        const unsigned nx = 5, ny = 5, nz = 5;
        auto vectorField = std::make_shared<VectorField3D>("Test Vector Field", nx, ny, nz);
        
        // Fill with some sample data
        for(unsigned x = 0; x < nx; x++) {
            for(unsigned y = 0; y < ny; y++) {
                for(unsigned z = 0; z < nz; z++) {
                    CFF::Utilities::Vector3u coords = CFF::Utilities::makeVector3u(x, y, z);
                    
                    // Create a simple spiral pattern
                    double fx = (double(x) - nx/2.0) / nx;
                    double fy = (double(y) - ny/2.0) / ny;
                    double fz = (double(z) - nz/2.0) / nz;
                    
                    CFF::Utilities::Vector3d vector = CFF::Utilities::makeVector3d(-fy, fx, fz * 0.5);
                    vectorField->setValue(coords, vector);
                }
            }
        }
        
        std::cout << "Created vector field with dimensions: " 
                  << nx << "x" << ny << "x" << nz << std::endl;
        
        // Create vector comparator
        VectorComparator comparator(vectorField.get());
        std::cout << "Vector comparator created successfully" << std::endl;
        
        // Test GridClusterer with VectorComparator would go here
        // Note: The GridClusterer is a template class and requires
        // specific instantiation for different data types
        
        std::cout << "Example completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
