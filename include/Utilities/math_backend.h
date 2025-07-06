#ifndef MATH_BACKEND_H
#define MATH_BACKEND_H

// This header provides a unified interface using Boost.uBLAS as the math backend
// (mathfu has been removed from the project)

#include <Utilities/boost_math_utilities.h>

namespace SNS { namespace Math {
    // Type aliases for unified interface
    using Vector3d = SNS::Utilities::Vector3d;
    using Vector4d = SNS::Utilities::Vector4d;
    using Matrix3d = SNS::Utilities::Matrix3d;
    using Matrix4d = SNS::Utilities::Matrix4d;
    using Vector3u = boost::numeric::ublas::vector<unsigned int>;
    
    // Helper functions
    inline Vector3d makeVector3d(double x, double y, double z) {
        return SNS::Utilities::makeVector3d(x, y, z);
    }
    
    inline Vector4d makeVector4d(double x, double y, double z, double w) {
        return SNS::Utilities::makeVector4d(x, y, z, w);
    }
    
    inline Vector3u makeVector3u(unsigned int x, unsigned int y, unsigned int z) {
        Vector3u v(3);
        v(0) = x; v(1) = y; v(2) = z;
        return v;
    }
    
    inline Matrix4d identityMatrix4d() {
        return SNS::Utilities::identityMatrix4d();
    }
    
    // Mathematical operations
    inline double dot(const Vector3d& a, const Vector3d& b) {
        return SNS::Utilities::dot(a, b);
    }
    
    inline Vector3d cross(const Vector3d& a, const Vector3d& b) {
        return SNS::Utilities::cross(a, b);
    }
    
    inline double norm(const Vector3d& v) {
        return SNS::Utilities::norm(v);
    }
    
    inline Vector3d normalize(const Vector3d& v) {
        return SNS::Utilities::normalize(v);
    }
    
    // Access elements
    inline double& x(Vector3d& v) { return v(0); }
    inline double& y(Vector3d& v) { return v(1); }
    inline double& z(Vector3d& v) { return v(2); }
    inline const double& x(const Vector3d& v) { return v(0); }
    inline const double& y(const Vector3d& v) { return v(1); }
    inline const double& z(const Vector3d& v) { return v(2); }
    
    inline unsigned int& x(Vector3u& v) { return v(0); }
    inline unsigned int& y(Vector3u& v) { return v(1); }
    inline unsigned int& z(Vector3u& v) { return v(2); }
    inline const unsigned int& x(const Vector3u& v) { return v(0); }
    inline const unsigned int& y(const Vector3u& v) { return v(1); }
    inline const unsigned int& z(const Vector3u& v) { return v(2); }
}}

#endif // MATH_BACKEND_H
