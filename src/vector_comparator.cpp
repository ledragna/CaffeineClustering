#include "SciData/vector_comparator.h"

using namespace CFF::SciData;

typedef boost::numeric::ublas::vector<double> dvec3;
typedef boost::numeric::ublas::matrix<double> dmat3;
typedef CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster Cluster;

VectorComparator::VectorComparator(const VectorField3D* vecField)
{
	_maxDist = 0.0;
	_maxMagnitude = 0.0;

	uint numPntX = vecField->getNumPointsX();
	uint numPntY = vecField->getNumPointsY();
	uint numPntZ = vecField->getNumPointsZ();

	for(uint x = 0; x < numPntX; ++x)
		for(uint y = 0; y < numPntY; ++y)
			for(uint z = 0; z < numPntZ; ++z)
			{				dvec3 v = vecField->getValue(CFF::Utilities::makeVector3u(x,y,z));
				double mag = CFF::Utilities::norm(v);
				_maxMagnitude = std::max(_maxMagnitude, mag);
			}

	dvec3 endpoint0 = vecField->fromLocalToWorldFrame(CFF::Utilities::makeVector3u(0,0,0));
	dvec3 endpoint1 = vecField->fromLocalToWorldFrame(CFF::Utilities::makeVector3u(numPntX-1,numPntY-1,numPntZ-1));
	_maxDist = CFF::Utilities::norm(endpoint1 - endpoint0);
}

double VectorComparator::operator()(const Cluster* cluster1, const Cluster* cluster2) const
{
	return  magnitudeComparator(cluster1, cluster2) +
			directionComparator(cluster1, cluster2) +
			distanceComparator(cluster1, cluster2);
}

bool VectorComparator::zeroValue(const Cluster* cluster1) const
{
	const double ZERO_VECTOR_PER_COMPONENT_THRESHOLD = 1e-9;
	return	(abs(cluster1->Value(0)) >= ZERO_VECTOR_PER_COMPONENT_THRESHOLD) ||
			(abs(cluster1->Value(1)) >= ZERO_VECTOR_PER_COMPONENT_THRESHOLD) ||
			(abs(cluster1->Value(2)) >= ZERO_VECTOR_PER_COMPONENT_THRESHOLD);
}

double VectorComparator::magnitudeComparator(const Cluster* cluster1, const Cluster* cluster2) const
{
	double magV1 = CFF::Utilities::norm(cluster1->Value);
	double magV2 = CFF::Utilities::norm(cluster2->Value);
	double res = std::abs(magV2 - magV1) / _maxMagnitude;
	return res;
}

double VectorComparator::distanceComparator(const Cluster* cluster1, const Cluster* cluster2) const
{
	return CFF::Utilities::norm(cluster1->Centroid - cluster2->Centroid) / _maxDist;
}

double VectorComparator::directionComparator(const Cluster* cluster1, const Cluster* cluster2) const
{
	double cosAlpha = CFF::Utilities::dot(cluster1->Value, cluster2->Value) /
									   (CFF::Utilities::norm(cluster1->Value) * CFF::Utilities::norm(cluster2->Value));
	return 0.5 * (1.0 - cosAlpha);
}

double VectorComparator::euclideanDirMagComparison(const Cluster* cluster1, const Cluster* cluster2) const
{
	dvec3 v1 = cluster1->Value;
	dvec3 v2 = cluster2->Value;
	return CFF::Utilities::norm(v2-v1) / std::min(CFF::Utilities::norm(v1), CFF::Utilities::norm(v2));
}

double VectorComparator::teleaDirMagComparison(const Cluster* cluster1, const Cluster* cluster2) const
{
	return std::max(teleaDirMagAsymmetric(cluster1,cluster2),
					teleaDirMagAsymmetric(cluster2,cluster1));
}

double VectorComparator::teleaDirMagAsymmetric(const Cluster* cluster1, const Cluster* cluster2) const
{
	dvec3 V1 = cluster1->Value;
	dvec3 V2 = cluster2->Value;

	/* Preparatory phase: express V2 in another reference frame (Xv, Yv, Zv)
			 * in which V1 lies on the X axis */

	dvec3 Xv = V1;
	// Computes 2 vectors orthogonal to Xv
	dvec3 Yv, Zv;
	CFF::Utilities::findOrthonormalBasis(Xv,Yv,Zv); //Note: it also cares to normalize Xv!
	
	// Create a 3x3 matrix for change of basis
	dmat3 T(3, 3);
	T(0,0) = Xv(0); T(0,1) = Xv(1); T(0,2) = Xv(2);
	T(1,0) = Yv(0); T(1,1) = Yv(1); T(1,2) = Yv(2);
	T(2,0) = Zv(0); T(2,1) = Zv(1); T(2,2) = Zv(2);
	
	// For simplicity, we'll use transpose as a rough approximation of inverse for orthogonal matrix
	dmat3 T_inv(3, 3);
	T_inv(0,0) = T(0,0); T_inv(0,1) = T(1,0); T_inv(0,2) = T(2,0);
	T_inv(1,0) = T(0,1); T_inv(1,1) = T(1,1); T_inv(1,2) = T(2,1);
	T_inv(2,0) = T(0,2); T_inv(2,1) = T(1,2); T_inv(2,2) = T(2,2);
	
	// Express V2 in the local reference frame of V1.
	// Let's call the resulting vector 'W'
	dvec3 W = boost::numeric::ublas::prod(T_inv, V2);

	/* Error computation */

	double v1_len = CFF::Utilities::norm(V1);

	// Params of the elliptic similarity function and their squared value
	double a = 0.5  * v1_len; // Radius of the ellipse along the V1 direction (Xv)
	double b = 0.25 * v1_len; // Radius of the ellipse along the a direction (Yv) orthogonal to Xv
	double c = 0.25 * v1_len; // Radius of the ellipse along the a direction (Zv) orthogonal to Xv and Yv
	double d = 0.25 * v1_len; // Offset detween the apex of Xv and the center of the ellipse

	double a2 = a * a;
	double b2 = b * b;
	double c2 = c * c;
	double d2 = d * d;

	double k1 = W(0) - v1_len;
	double k2 = a2 - d2;

	double p1 = (a2 * k1 * k1) + (k2 * ((a2*W(1)*W(1)/b2) + (a2*W(2)*W(2)/c2)));
	double p2 = d * k1 / k2;

	double err = (sqrt(p1) / k2) - p2;
	return err;
}

double VectorComparator::volumeComparison(const Cluster* cluster1, const Cluster* cluster2) const
{
	double vol1 = (double) cluster1->Volume;
	double vol2 = (double) cluster2->Volume;
	return std::abs(vol2-vol1) / std::min(vol2,vol1);
}
