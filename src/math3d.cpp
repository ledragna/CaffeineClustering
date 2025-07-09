#include "Utilities/math3d.h"

#include <cmath>

#include <Utilities/boost_math_utilities.h>

// Using declarations for convenience
using CFF::Utilities::Vector3d;
using CFF::Utilities::makeVector3d;
using CFF::Utilities::normalize;
using CFF::Utilities::cross;
using CFF::Utilities::dot;
using CFF::Utilities::fuzzyEqual;


bool CFF::Utilities::rayAABBIntersection( const Vector3d& rayOrigin,
												const Vector3d& rayDir,
												const Vector3d& aabbMinXYZ,
												const Vector3d& aabbMaxXYZ,
												double& outTmin, double& outTmax)
{
	const double EPSILON = 1e-4;

	outTmin = std::numeric_limits<double>::lowest();
	outTmax = std::numeric_limits<double>::max();

	// For each pair of planes orthogonal to X, Y and Z...
	for(short i = 0; i < 3; i++)
	{
		if(std::abs(rayDir[i]) < EPSILON)
		{
			// The ray is orthogonal to the i-th axis.
			// The ray may intersect the box only if the origin lies
			// within the considered pair of planes.
			if((rayOrigin[i] < aabbMinXYZ[i]) || (rayOrigin[i] > aabbMaxXYZ[i]))
				return false;
		}
		else
		{
			// Computes the intersection "times" between the ray and the
			// pair of plains orthogonal to the i-th axis.
			double inv_dir = 1.0 / rayDir[i];
			double t_i_min =  (aabbMinXYZ[i] - rayOrigin[i]) * inv_dir;
			double t_i_max =  (aabbMaxXYZ[i] - rayOrigin[i]) * inv_dir;
			if(t_i_max < t_i_min) std::swap(t_i_min, t_i_max);

			// If needed, updates outTmin and outTmax
			outTmin = std::max(outTmin, t_i_min);
			outTmax = std::min(outTmax, t_i_max);

			// If the intersection between the 3 intervals
			// [t_i_min, t_i_max] is empty, the ray miss the box.
			if(outTmin > outTmax) return false;
		}
	}

	// The ray hits the box. The two intersection "times"
	// are stored in outTmin and outTmax.
	return true;
}


bool CFF::Utilities::rayPlaneIntersection(const Vector3d& rayOrigin,
												const Vector3d& rayDir,
												const Vector3d& planePoint,
												const Vector3d& planeNormal,
												double& out_t)
{
	typedef Vector3d dvec3;
	dvec3 N = normalize(planeNormal);
	double NdotD = dot(N, rayDir);
	// If the ray is (almost) parallel to the plane, there is no intersection
	const double EPSILON = 1e-5;
	if(std::fabs(NdotD) < EPSILON) return false;
	out_t = dot(planePoint - rayOrigin, N) / NdotD;
	return true;
}


void CFF::Utilities::findOrthonormalBasis(Vector3d& inoutPseudoZ,
										  Vector3d& outPseudoY,
										  Vector3d& outPseudoX)
{
	typedef Vector3d dvec3;
	if(fuzzyEqual(inoutPseudoZ, makeVector3d(0,0,0), 1e-12))
		throwAndPrintError<std::invalid_argument>(
			"Utilities::findOrthogonalBase : the input is a null vector!"
		);

	dvec3 v1 = normalize(inoutPseudoZ);

	// Finds, among X/Y/Z, the axis which "most orthogonal" to inoutVector1
	double v1XAbs = fabs(v1(0));
	double v1YAbs = fabs(v1(1));
	double v1ZAbs = fabs(v1(2));
	dvec3 almostOrthogonalToV1;
	if( (v1XAbs <= v1YAbs) && (v1XAbs <= v1ZAbs) )
		almostOrthogonalToV1 = makeVector3d(1,0,0);
	else if( (v1YAbs <= v1XAbs) && (v1YAbs <= v1ZAbs) )
		almostOrthogonalToV1 = makeVector3d(0,1,0);
	else
		almostOrthogonalToV1 = makeVector3d(0,0,1);
	dvec3 v2 = normalize(cross(v1, almostOrthogonalToV1));
	dvec3 v3 = normalize(cross(v2, v1));

	inoutPseudoZ = v1;
	outPseudoY = v2;
	outPseudoX = v3;
}



