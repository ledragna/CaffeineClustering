#include "Utilities/boost_math_utilities.h"

namespace SNS { namespace Utilities
{
    QVector3D to_QVector3D(const Vector3d& p) {
        return QVector3D(p(0), p(1), p(2));
    }
    
    QVector4D to_QVector4D(const Vector4d& p) {
        return QVector4D(p(0), p(1), p(2), p(3));
    }
    
    QMatrix4x4 to_QMatrix4x4(const Matrix4d& mat) {
        QMatrix4x4 result;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                result(i, j) = mat(i, j);
            }
        }
        return result;
    }
    
    Vector3d to_Boost_Vector3d(const QVector3D& p) {
        return makeVector3d(p.x(), p.y(), p.z());
    }
    
    Vector4d to_Boost_Vector4d(const QVector4D& p) {
        return makeVector4d(p.x(), p.y(), p.z(), p.w());
    }
    
    Matrix4d to_Boost_Matrix4d(const QMatrix4x4& mat) {
        Matrix4d result(4, 4);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                result(i, j) = mat(i, j);
            }
        }
        return result;
    }
} }
