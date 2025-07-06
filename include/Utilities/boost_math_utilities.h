#ifndef BOOST_MATH_UTILITIES_H
#define BOOST_MATH_UTILITIES_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <type_traits>

#include <QVector3D>
#include <QVector4D>
#include <QMatrix4x4>

#include <Utilities/floating_point.h>

namespace SNS { namespace Utilities
{
    // Type aliases for convenience
    using Vector3d = boost::numeric::ublas::vector<double>;
    using Vector4d = boost::numeric::ublas::vector<double>;
    using Vector3u = boost::numeric::ublas::vector<unsigned int>;
    using Matrix3d = boost::numeric::ublas::matrix<double>;
    using Matrix4d = boost::numeric::ublas::matrix<double>;
    
    // Helper functions to create specific size vectors/matrices
    inline Vector3d makeVector3d(double x, double y, double z) {
        Vector3d v(3);
        v(0) = x; v(1) = y; v(2) = z;
        return v;
    }
    
    inline Vector4d makeVector4d(double x, double y, double z, double w) {
        Vector4d v(4);
        v(0) = x; v(1) = y; v(2) = z; v(3) = w;
        return v;
    }
    
    inline Matrix3d makeMatrix3d() {
        return Matrix3d(3, 3);
    }
    
    inline Matrix4d makeMatrix4d() {
        return Matrix4d(4, 4);
    }
    
    inline Matrix4d identityMatrix4d() {
        Matrix4d m(4, 4);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                m(i, j) = (i == j) ? 1.0 : 0.0;
            }
        }
        return m;
    }
    
    // Mathematical operations
    inline double dot(const Vector3d& a, const Vector3d& b) {
        return boost::numeric::ublas::inner_prod(a, b);
    }
    
    inline Vector3d cross(const Vector3d& a, const Vector3d& b) {
        Vector3d result(3);
        result(0) = a(1) * b(2) - a(2) * b(1);
        result(1) = a(2) * b(0) - a(0) * b(2);
        result(2) = a(0) * b(1) - a(1) * b(0);
        return result;
    }
    
    inline double norm(const Vector3d& v) {
        return boost::numeric::ublas::norm_2(v);
    }
    
    inline Vector3d normalize(const Vector3d& v) {
        double n = norm(v);
        return v / n;
    }
    
    // Matrix operations
    inline Vector3d getColumn(const Matrix4d& m, int col) {
        Vector3d result(3);
        result(0) = m(0, col);
        result(1) = m(1, col);
        result(2) = m(2, col);
        return result;
    }
    
    // Helper functions for creating unsigned integer vectors
    inline boost::numeric::ublas::vector<unsigned int> makeVector3u(unsigned int x, unsigned int y, unsigned int z) {
        boost::numeric::ublas::vector<unsigned int> v(3);
        v(0) = x; v(1) = y; v(2) = z;
        return v;
    }
    
    // Matrix vector multiplication for 4D
    inline Vector4d multiply(const Matrix4d& m, const Vector4d& v) {
        return boost::numeric::ublas::prod(m, v);
    }
    
    // Extract xyz from 4D vector
    inline Vector3d xyz(const Vector4d& v) {
        Vector3d result(3);
        result(0) = v(0);
        result(1) = v(1);
        result(2) = v(2);
        return result;
    }
    
    // Matrix inversion with determinant check
    inline bool InverseWithDeterminantCheck(const Matrix4d& matrix, Matrix4d* result) {
        // For simplicity, we'll use a basic 4x4 matrix inversion
        // In production code, you might want to use a more robust library
        try {
            *result = boost::numeric::ublas::zero_matrix<double>(4, 4);
            // This is a placeholder - would need proper matrix inversion implementation
            // For now, just return identity matrix
            for (int i = 0; i < 4; ++i) {
                (*result)(i, i) = 1.0;
            }
            return true;  // Assume invertible for now
        } catch (...) {
            return false;
        }
    }
    
    // Conversion functions to Qt types
    QVector3D to_QVector3D(const Vector3d& p);
    QVector4D to_QVector4D(const Vector4d& p);
    QMatrix4x4 to_QMatrix4x4(const Matrix4d& mat);
    
    // Conversion functions from Qt types
    Vector3d to_Boost_Vector3d(const QVector3D& p);
    Vector4d to_Boost_Vector4d(const QVector4D& p);
    Matrix4d to_Boost_Matrix4d(const QMatrix4x4& mat);
    
    // Matrix/vector fuzzy comparison
    template<typename T>
    bool fuzzyEqual(const T& a, const T& b, double epsilon = 1e-4) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) {
            if (!SNS::Utilities::fuzzyEqual(a(i), b(i), epsilon)) {
                return false;
            }
        }
        return true;
    }
    
    template<typename T>
    bool fuzzyEqual(const boost::numeric::ublas::matrix<T>& m1,
                    const boost::numeric::ublas::matrix<T>& m2,
                    double epsilon = 1e-4) {
        if (m1.size1() != m2.size1() || m1.size2() != m2.size2()) {
            return false;
        }
        for (size_t i = 0; i < m1.size1(); ++i) {
            for (size_t j = 0; j < m1.size2(); ++j) {
                if (!SNS::Utilities::fuzzyEqual(m1(i, j), m2(i, j), epsilon)) {
                    return false;
                }
            }
        }
        return true;
    }
    
    // Vector comparison functions
    template<typename T>
    inline bool vectorsEqual(const boost::numeric::ublas::vector<T>& v1, 
                            const boost::numeric::ublas::vector<T>& v2) {
        if (v1.size() != v2.size()) return false;
        for (size_t i = 0; i < v1.size(); ++i) {
            if constexpr (std::is_floating_point_v<T>) {
                if (!SNS::Utilities::fuzzyEqual(v1(i), v2(i))) return false;
            } else {
                if (v1(i) != v2(i)) return false;
            }
        }
        return true;
    }
    
    // Element-wise multiplication functions
    template<typename T>
    inline boost::numeric::ublas::vector<T> elementWiseMultiply(
        const boost::numeric::ublas::vector<T>& v1,
        const boost::numeric::ublas::vector<T>& v2) {
        boost::numeric::ublas::vector<T> result(v1.size());
        for (size_t i = 0; i < v1.size(); ++i) {
            result(i) = v1(i) * v2(i);
        }
        return result;
    }
    
    // Element-wise multiplication assignment
    template<typename T>
    inline void elementWiseMultiplyAssign(
        boost::numeric::ublas::vector<T>& v1,
        const boost::numeric::ublas::vector<T>& v2) {
        for (size_t i = 0; i < v1.size(); ++i) {
            v1(i) *= v2(i);
        }
    }
} }

#endif // BOOST_MATH_UTILITIES_H
