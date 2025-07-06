#include "Utilities/mathfu_utilities.h"


QVector3D SNS::Utilities::to_QVector3D(const mathfu::Vector<double,3>& p)
{
	return QVector3D(p.x, p.y, p.z);
}

QVector4D SNS::Utilities::to_QVector4D(const mathfu::Vector<double,4>& p)
{
	return QVector4D(p.x, p.y, p.z, p.w);
}

QMatrix4x4 SNS::Utilities::to_QMatrix4x4(const mathfu::Matrix<double,4>& mat)
{
	return QMatrix4x4(	mat(0,0), mat(0,1), mat(0,2), mat(0,3),
						mat(1,0), mat(1,1), mat(1,2), mat(1,3),
						mat(2,0), mat(2,1), mat(2,2), mat(2,3),
						mat(3,0), mat(3,1), mat(3,2), mat(3,3));
}

mathfu::Vector<double,3> SNS::Utilities::to_Mathfu_dvec3(const QVector3D& p)
{
	return mathfu::Vector<double,3>(p.x(), p.y(), p.z());
}

mathfu::Vector<double,4> SNS::Utilities::to_Mathfu_dvec4(const QVector4D& p)
{
	return mathfu::Vector<double,4>(p.x(), p.y(), p.z(), p.w());
}
