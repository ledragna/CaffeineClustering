#include "Utilities/math3d.h"

#include <QVector4D>
#include <QMatrix4x4>
#include <QRectF>

#include <cmath>

#include <Utilities/boost_math_utilities.h>

// Using declarations for convenience
using SNS::Utilities::Vector3d;
using SNS::Utilities::makeVector3d;
using SNS::Utilities::normalize;
using SNS::Utilities::cross;
using SNS::Utilities::dot;
using SNS::Utilities::fuzzyEqual;

QColor SNS::Utilities::mixQColors(double t,
										const QColor& startColor,
										const QColor& endColor)
{
	double r = (t * endColor.redF())   + ((1.0 - t) * startColor.redF());
	double g = (t * endColor.greenF()) + ((1.0 - t) * startColor.greenF());
	double b = (t * endColor.blueF())  + ((1.0 - t) * startColor.blueF());

	QColor finalColor;
	finalColor.setRedF(   qBound(0.0, r, 1.0) );
	finalColor.setGreenF( qBound(0.0, g, 1.0) );
	finalColor.setBlueF(  qBound(0.0, b, 1.0) );
	return finalColor;
}



QVector3D SNS::Utilities::BSplinePoint(const QVector3D controlPoints[], float t)
{
	// Cubic B-Spline equation:
	// P(t) = [ t^3  t^2  t  1 ] * M * C
	// where M (defined below) is a 4x4 matrix that characterizes this type of curve
	// and C is a 4x4 matrix whose rows are the 4 control points.

	QMatrix4x4 M;
	M.setColumn(0, QVector4D(-1,  3, -3,  1));
	M.setColumn(1, QVector4D( 3, -6,  0,  4));
	M.setColumn(2, QVector4D(-3,  3,  3,  1));
	M.setColumn(3, QVector4D( 1,  0,  0,  0));
	M *= (1.0f / 6.0f);

	QMatrix4x4 C;
	C.setRow(0, QVector4D(controlPoints[0], 1));
	C.setRow(1, QVector4D(controlPoints[1], 1));
	C.setRow(2, QVector4D(controlPoints[2], 1));
	C.setRow(3, QVector4D(controlPoints[3], 1));

	QMatrix4x4 BS = M * C;

	QVector4D t_vect(t*t*t, t*t, t, 1);
	QVector4D curvePoint = t_vect * BS;

	return curvePoint.toVector3D();
}

// Return nsegments + 1 points
// nsegments must be > 0
void SNS::Utilities::BSplinePath(const QVector3D controlPoints[],
									   unsigned int nsegments,
									   QVector<QVector3D> &outResult)
{
	if(!nsegments) throw std::invalid_argument("BSplinePath: The number "
											   "of segments must be > 0.");

	outResult.resize(nsegments + 1);

	// Cubic B-Spline equation:
	// P(t) = [ t^3  t^2  t  1 ] * M * C
	// where M (defined below) is a 4x4 matrix that characterizes this type of curve
	// and C is a 4x4 matrix whose rows are the 4 control points.

	QMatrix4x4 M;
	M.setColumn(0, QVector4D(-1,  3, -3,  1));
	M.setColumn(1, QVector4D( 3, -6,  0,  4));
	M.setColumn(2, QVector4D(-3,  3,  3,  1));
	M.setColumn(3, QVector4D( 1,  0,  0,  0));
	M *= (1.0f / 6.0f);

	QMatrix4x4 C;
	C.setRow(0, QVector4D(controlPoints[0], 1));
	C.setRow(1, QVector4D(controlPoints[1], 1));
	C.setRow(2, QVector4D(controlPoints[2], 1));
	C.setRow(3, QVector4D(controlPoints[3], 1));

	QMatrix4x4 BS = M * C;

	float t_step = 1.0 / nsegments;
	for(unsigned int i = 0; i < nsegments + 1; i++)
	{
		double t = i * t_step;
		QVector4D t_vect(t*t*t, t*t, t, 1);
		QVector4D curvePoint = t_vect * BS;
		outResult[i] = curvePoint.toVector3D();
	}
}


bool SNS::Utilities::rayAABBIntersection( const Vector3d& rayOrigin,
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


bool SNS::Utilities::rayPlaneIntersection(const Vector3d& rayOrigin,
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

bool SNS::Utilities::rayPlaneIntersection(const QVector3D& rayOrigin,
										  const QVector3D& rayDir,
										  const QVector3D& planePoint,
										  const QVector3D& planeNormal,
										  float& out_t)
{
	QVector3D N = planeNormal.normalized();
	float NdotD = QVector3D::dotProduct(N, rayDir);
	// If the ray is (almost) parallel to the plane, there is no intersection
	const double EPSILON = 1e-5;
	if(std::fabs(NdotD) < EPSILON) return false;

	out_t = QVector3D::dotProduct(planePoint - rayOrigin, N) / NdotD;
	return true;
}


void SNS::Utilities::findOrthonormalBasis(QVector3D& inoutPseudoZ,
										  QVector3D& outPseudoY,
										  QVector3D& outPseudoX)
{
	if(fuzzyEqual(inoutPseudoZ, QVector3D(0,0,0), 1e-5f))
		throwAndPrintError<std::invalid_argument>(
			"Utilities::findOrthogonalBase : the input is a null vector!"
		);

	QVector3D v1 = inoutPseudoZ.normalized();

	// Finds, among X/Y/Z, the axis which "most orthogonal" to inoutVector1
	float v1XAbs = fabs(v1.x());
	float v1YAbs = fabs(v1.y());
	float v1ZAbs = fabs(v1.z());
	QVector3D almostOrthogonalToV1;
	if( (v1XAbs <= v1YAbs) && (v1XAbs <= v1ZAbs) )
		almostOrthogonalToV1 = QVector3D(1,0,0);
	else if( (v1YAbs <= v1XAbs) && (v1YAbs <= v1ZAbs) )
		almostOrthogonalToV1 = QVector3D(0,1,0);
	else
		almostOrthogonalToV1 = QVector3D(0,0,1);

	QVector3D v2 = QVector3D::crossProduct(v1, almostOrthogonalToV1).normalized();
	QVector3D v3 = QVector3D::crossProduct(v2, v1).normalized();

	inoutPseudoZ = v1;
	outPseudoY = v2;
	outPseudoX = v3;
}

void SNS::Utilities::findOrthonormalBasis(Vector3d& inoutPseudoZ,
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



QVector3D SNS::Utilities::unproject(const QVector3D& windowCoords,
									const QRectF& viewPort,
									const QMatrix4x4& MP)
{
	QVector3D winPnt(
		windowCoords.x(),
		windowCoords.y(),
		std::min(1.0f, std::max(0.0f, windowCoords.z()))
	);

	QVector4D ndc; // Normalized Device Coordinates
	ndc.setX((2.0f * (winPnt.x() - viewPort.x()) / viewPort.width()) -1.0f);
	ndc.setY((2.0f * (winPnt.y() - viewPort.y()) / viewPort.height()) -1.0f);
	ndc.setZ((2.0f * winPnt.z()) -1.0f);
	ndc.setW(1);

	QVector4D worldSpacePnt = MP.inverted() * ndc;
	return worldSpacePnt.toVector3DAffine();
}



SNS::Utilities::AxesSimilarity SNS::Utilities::computeGlobalAxesSimilarity(
														const QVector3D& v1,
														const QVector3D& v2,
														const QVector3D& v3)
{
	const QVector3D normV1 = v1.normalized();
	const QVector3D normV2 = v2.normalized();
	const QVector3D normV3 = v3.normalized();

	AxesSimilarity results;
	if( (fabs(normV1.x()) >= fabs(normV2.x())) &&
		(fabs(normV1.x()) >= fabs(normV3.x())) )
	{
		results.isV1SimilarToX = std::signbit(normV1.x()) ? -1 : 1;
		results.isV2SimilarToX = 0;
		results.isV3SimilarToX = 0;

		if(fabs(normV2.y()) >= fabs(normV3.y()))
		{
			results.isV1SimilarToY = 0;
			results.isV2SimilarToY = std::signbit(normV2.y()) ? -1 : 1;
			results.isV3SimilarToY = 0;

			results.isV1SimilarToZ = 0;
			results.isV2SimilarToZ = 0;
			results.isV3SimilarToZ = std::signbit(normV3.z()) ? -1 : 1;
		}
		else
		{
			results.isV1SimilarToZ = 0;
			results.isV2SimilarToZ = std::signbit(normV2.z()) ? -1 : 1;
			results.isV3SimilarToZ = 0;

			results.isV1SimilarToY = 0;
			results.isV2SimilarToY = 0;
			results.isV3SimilarToY = std::signbit(normV3.y()) ? -1 : 1;
		}
	}
	else if( (fabs(normV1.y()) >= fabs(normV2.y())) &&
			 (fabs(normV1.y()) >= fabs(normV3.y())) )
	{
		results.isV1SimilarToY = std::signbit(normV1.y()) ? -1 : 1;
		results.isV2SimilarToY = 0;
		results.isV3SimilarToY = 0;

		if(fabs(normV2.x()) >= fabs(normV3.x()))
		{
			results.isV1SimilarToX = 0;
			results.isV2SimilarToX = std::signbit(normV2.x()) ? -1 : 1;
			results.isV3SimilarToX = 0;

			results.isV1SimilarToZ = 0;
			results.isV2SimilarToZ = 0;
			results.isV3SimilarToZ = std::signbit(normV3.z()) ? -1 : 1;
		}
		else
		{
			results.isV1SimilarToZ = 0;
			results.isV2SimilarToZ = std::signbit(normV2.z()) ? -1 : 1;
			results.isV3SimilarToZ = 0;

			results.isV1SimilarToX = 0;
			results.isV2SimilarToX = 0;
			results.isV3SimilarToX = std::signbit(normV3.x()) ? -1 : 1;
		}
	}
	else
	{
		results.isV1SimilarToZ = std::signbit(normV1.z()) ? -1 : 1;
		results.isV2SimilarToZ = 0;
		results.isV3SimilarToZ = 0;

		if(fabs(normV2.x()) >= fabs(normV3.x()))
		{
			results.isV1SimilarToX = 0;
			results.isV2SimilarToX = std::signbit(normV2.x()) ? -1 : 1;
			results.isV3SimilarToX = 0;

			results.isV1SimilarToY = 0;
			results.isV2SimilarToY = 0;
			results.isV3SimilarToY = std::signbit(normV3.y()) ? -1 : 1;
		}
		else
		{
			results.isV1SimilarToY = 0;
			results.isV2SimilarToY = std::signbit(normV2.y()) ? -1 : 1;
			results.isV3SimilarToY = 0;

			results.isV1SimilarToX = 0;
			results.isV2SimilarToX = 0;
			results.isV3SimilarToX = std::signbit(normV3.x()) ? -1 : 1;
		}
	}

	return results;
}



bool SNS::Utilities::closestPointsBetweenLines(const QVector3D& ray1Origin,
											   const QVector3D& ray1Dir,
											   const QVector3D& ray2Origin,
											   const QVector3D& ray2Dir,
											   float* out_t1, float* out_t2)
{
	QVector3D w = ray1Origin - ray2Origin;
	float a = QVector3D::dotProduct(ray1Dir,ray1Dir);
	float b = QVector3D::dotProduct(ray1Dir,ray2Dir);
	float c = QVector3D::dotProduct(ray2Dir,ray2Dir);
	float d = QVector3D::dotProduct(ray1Dir,w);
	float e = QVector3D::dotProduct(ray2Dir,w);

	float denom = (a*c)-(b*b);
	if(fuzzyEqual(0.0f, denom, 1e-5f))
	{
		if(out_t1) *out_t1 = 0.0f;
		if(out_t2) *out_t2 = 0.0f;
		return false;
	}

	float t1 = ((b*e)-(c*d)) / denom;
	if(out_t1) *out_t1 = t1;

	float t2 = ((a*e)-(b*d)) / denom;
	if(out_t2) *out_t2 = t2;

	return true;
}
