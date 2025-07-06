#include "Utilities/floating_point.h"

#include <cmath>

bool SNS::Utilities::fuzzyEqual(double x, double y, double epsilon)
{
	return	((fabs(x - y)/fabs(x == 0.0 ? 1.0 : x)) <= epsilon) ||
			((fabs(x - y)/fabs(y == 0.0 ? 1.0 : y)) <= epsilon);
}

bool SNS::Utilities::fuzzyEqual(float x, float y, float epsilon)
{
	return	((fabs(x - y)/fabs(x == 0.0f ? 1.0f : x)) <= epsilon) ||
			((fabs(x - y)/fabs(y == 0.0f ? 1.0f : y)) <= epsilon);
}

bool SNS::Utilities::fuzzyEqual(const double& x, const double& y)
{
	return fuzzyEqual(x,y,1e-4);
}

bool SNS::Utilities::fuzzyEqual(const float& x, const float& y)
{
	return fuzzyEqual(x,y,1e-4f);
}


bool SNS::Utilities::fuzzyLessThan(double x, double y, double epsilon)
{
	return (x < y) && (!fuzzyEqual(x, y, epsilon));
}

bool SNS::Utilities::fuzzyLessThan(float x, float y, float epsilon)
{
	return (x < y) && (!fuzzyEqual(x, y, epsilon));
}

bool SNS::Utilities::fuzzyLessThan(const double& x, const double& y)
{
	return fuzzyLessThan(x,y,1e-4);
}

bool SNS::Utilities::fuzzyLessThan(const float& x, const float& y)
{
	return fuzzyLessThan(x,y,1e-4f);
}


bool SNS::Utilities::fuzzyEqualOrLessThan(double x, double y, double epsilon)
{
	return (x < y) || fuzzyEqual(x, y, epsilon);
}

bool SNS::Utilities::fuzzyEqualOrLessThan(float x, float y, float epsilon)
{
	return (x < y) || fuzzyEqual(x, y, epsilon);
}

bool SNS::Utilities::fuzzyEqualOrLessThan(const double& x, const double& y)
{
	return fuzzyEqualOrLessThan(x,y,1e-4);
}

bool SNS::Utilities::fuzzyEqualOrLessThan(const float& x, const float& y)
{
	return fuzzyEqualOrLessThan(x,y,1e-4f);
}


bool SNS::Utilities::fuzzyGreaterThan(double x, double y, double epsilon)
{
	return (x > y) && (!fuzzyEqual(x, y, epsilon));
}

bool SNS::Utilities::fuzzyGreaterThan(float x, float y, float epsilon)
{
	return (x > y) && (!fuzzyEqual(x, y, epsilon));
}

bool SNS::Utilities::fuzzyGreaterThan(const double& x, const double& y)
{
	return fuzzyGreaterThan(x,y,1e-4);
}

bool SNS::Utilities::fuzzyGreaterThan(const float& x, const float& y)
{
	return fuzzyGreaterThan(x,y,1e-4f);
}


bool SNS::Utilities::fuzzyEqualOrGreaterThan(double x, double y, double epsilon)
{
	return (x > y) || fuzzyEqual(x, y, epsilon);
}

bool SNS::Utilities::fuzzyEqualOrGreaterThan(float x, float y, float epsilon)
{
	return (x > y) || fuzzyEqual(x, y, epsilon);
}

bool SNS::Utilities::fuzzyEqualOrGreaterThan(const double& x, const double& y)
{
	return fuzzyEqualOrGreaterThan(x,y,1e-4);
}

bool SNS::Utilities::fuzzyEqualOrGreaterThan(const float& x, const float& y)
{
	return fuzzyEqualOrGreaterThan(x,y,1e-4f);
}

bool SNS::Utilities::fuzzyEqual(const QVector3D& v1,
									  const QVector3D& v2,
									  float epsilon)
{
	return	fuzzyEqual(v1.x(), v2.x(), epsilon) &&
			fuzzyEqual(v1.y(), v2.y(), epsilon) &&
			fuzzyEqual(v1.z(), v2.z(), epsilon);
}

bool SNS::Utilities::fuzzyEqual(const QVector3D& v1, const QVector3D& v2)
{
	return	fuzzyEqual(v1.x(), v2.x(), 1e-4f) &&
			fuzzyEqual(v1.y(), v2.y(), 1e-4f) &&
			fuzzyEqual(v1.z(), v2.z(), 1e-4f);
}

bool SNS::Utilities::fuzzyEqual(const QVector4D& v1,
								const QVector4D& v2,
								float epsilon)
{
	return	fuzzyEqual(v1.x(), v2.x(), epsilon) &&
			fuzzyEqual(v1.y(), v2.y(), epsilon) &&
			fuzzyEqual(v1.z(), v2.z(), epsilon) &&
			fuzzyEqual(v1.w(), v2.w(), epsilon);
}

bool SNS::Utilities::fuzzyEqual(const QVector4D& v1, const QVector4D& v2)
{
	return	fuzzyEqual(v1.x(), v2.x(), 1e-4f) &&
			fuzzyEqual(v1.y(), v2.y(), 1e-4f) &&
			fuzzyEqual(v1.z(), v2.z(), 1e-4f) &&
			fuzzyEqual(v1.w(), v2.w(), 1e-4f);
}


bool SNS::Utilities::fuzzyEqual(const QMatrix4x4& m1,
								const QMatrix4x4& m2,
								float epsilon)
{
	bool equals = true;
	for(int i = 0; (i < 4) && equals; ++i)
		for(int j = 0; (j < 4) && equals; ++j)
			equals = equals &&
					 SNS::Utilities::fuzzyEqual(m1(i,j), m2(i,j), epsilon);

	return equals;
}

bool SNS::Utilities::fuzzyEqual(const QMatrix4x4& m1, const QMatrix4x4& m2)
{
	return SNS::Utilities::fuzzyEqual(m1, m2, 1e-4f);
}
