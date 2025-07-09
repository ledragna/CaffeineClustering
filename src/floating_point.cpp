#include "Utilities/floating_point.h"

#include <cmath>

bool CFF::Utilities::fuzzyEqual(double x, double y, double epsilon)
{
	return	((fabs(x - y)/fabs(x == 0.0 ? 1.0 : x)) <= epsilon) ||
			((fabs(x - y)/fabs(y == 0.0 ? 1.0 : y)) <= epsilon);
}

bool CFF::Utilities::fuzzyEqual(float x, float y, float epsilon)
{
	return	((fabs(x - y)/fabs(x == 0.0f ? 1.0f : x)) <= epsilon) ||
			((fabs(x - y)/fabs(y == 0.0f ? 1.0f : y)) <= epsilon);
}

bool CFF::Utilities::fuzzyEqual(const double& x, const double& y)
{
	return fuzzyEqual(x,y,1e-4);
}

bool CFF::Utilities::fuzzyEqual(const float& x, const float& y)
{
	return fuzzyEqual(x,y,1e-4f);
}


bool CFF::Utilities::fuzzyLessThan(double x, double y, double epsilon)
{
	return (x < y) && (!fuzzyEqual(x, y, epsilon));
}

bool CFF::Utilities::fuzzyLessThan(float x, float y, float epsilon)
{
	return (x < y) && (!fuzzyEqual(x, y, epsilon));
}

bool CFF::Utilities::fuzzyLessThan(const double& x, const double& y)
{
	return fuzzyLessThan(x,y,1e-4);
}

bool CFF::Utilities::fuzzyLessThan(const float& x, const float& y)
{
	return fuzzyLessThan(x,y,1e-4f);
}


bool CFF::Utilities::fuzzyEqualOrLessThan(double x, double y, double epsilon)
{
	return (x < y) || fuzzyEqual(x, y, epsilon);
}

bool CFF::Utilities::fuzzyEqualOrLessThan(float x, float y, float epsilon)
{
	return (x < y) || fuzzyEqual(x, y, epsilon);
}

bool CFF::Utilities::fuzzyEqualOrLessThan(const double& x, const double& y)
{
	return fuzzyEqualOrLessThan(x,y,1e-4);
}

bool CFF::Utilities::fuzzyEqualOrLessThan(const float& x, const float& y)
{
	return fuzzyEqualOrLessThan(x,y,1e-4f);
}


bool CFF::Utilities::fuzzyGreaterThan(double x, double y, double epsilon)
{
	return (x > y) && (!fuzzyEqual(x, y, epsilon));
}

bool CFF::Utilities::fuzzyGreaterThan(float x, float y, float epsilon)
{
	return (x > y) && (!fuzzyEqual(x, y, epsilon));
}

bool CFF::Utilities::fuzzyGreaterThan(const double& x, const double& y)
{
	return fuzzyGreaterThan(x,y,1e-4);
}

bool CFF::Utilities::fuzzyGreaterThan(const float& x, const float& y)
{
	return fuzzyGreaterThan(x,y,1e-4f);
}


bool CFF::Utilities::fuzzyEqualOrGreaterThan(double x, double y, double epsilon)
{
	return (x > y) || fuzzyEqual(x, y, epsilon);
}

bool CFF::Utilities::fuzzyEqualOrGreaterThan(float x, float y, float epsilon)
{
	return (x > y) || fuzzyEqual(x, y, epsilon);
}

bool CFF::Utilities::fuzzyEqualOrGreaterThan(const double& x, const double& y)
{
	return fuzzyEqualOrGreaterThan(x,y,1e-4);
}

bool CFF::Utilities::fuzzyEqualOrGreaterThan(const float& x, const float& y)
{
	return fuzzyEqualOrGreaterThan(x,y,1e-4f);
}

// Note: Qt-specific vector and matrix comparison functions have been removed.
// The library now uses only standard C++ and Boost types.
