#include "SciData/vector_comparator.h"

using namespace SNS::SciData;

typedef mathfu::Vector<double,3> dvec3;
typedef mathfu::Matrix<double,3,3> dmat3;
typedef SNS::SciData::GridClusterer<mathfu::Vector<double,3>>::Cluster Cluster;

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
			{
				dvec3 v = vecField->getValue(mathfu::Vector<uint,3>(x,y,z));
				double mag = v.Length();
				_maxMagnitude = std::max(_maxMagnitude, mag);
			}

	dvec3 endpoint0 = vecField->fromLocalToWorldFrame(mathfu::Vector<uint,3>(0,0,0));
	dvec3 endpoint1 = vecField->fromLocalToWorldFrame(mathfu::Vector<uint,3>(numPntX-1,numPntY-1,numPntZ-1));
	_maxDist = (endpoint1 - endpoint0).Length();
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
	return	(abs(cluster1->Value.x) >= ZERO_VECTOR_PER_COMPONENT_THRESHOLD) ||
			(abs(cluster1->Value.y) >= ZERO_VECTOR_PER_COMPONENT_THRESHOLD) ||
			(abs(cluster1->Value.z) >= ZERO_VECTOR_PER_COMPONENT_THRESHOLD);
}

double VectorComparator::magnitudeComparator(const Cluster* cluster1, const Cluster* cluster2) const
{
	double magV1 = cluster1->Value.Length();
	double magV2 = cluster2->Value.Length();
	double res = std::abs(magV2 - magV1) / _maxMagnitude;
	return res;
}

double VectorComparator::distanceComparator(const Cluster* cluster1, const Cluster* cluster2) const
{
	return (cluster1->Centroid - cluster2->Centroid).Length() / _maxDist;
}

double VectorComparator::directionComparator(const Cluster* cluster1, const Cluster* cluster2) const
{
	double cosAlpha = dvec3::DotProduct(cluster1->Value, cluster2->Value) /
									   (cluster1->Value.Length() * cluster2->Value.Length());
	return 0.5 * (1.0 - cosAlpha);
}

double VectorComparator::euclideanDirMagComparison(const Cluster* cluster1, const Cluster* cluster2) const
{
	dvec3 v1 = cluster1->Value;
	dvec3 v2 = cluster2->Value;
	return (v2-v1).Length() / std::min(v1.Length(), v2.Length());
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
	SNS::Utilities::findOrthonormalBasis(Xv,Yv,Zv); //Note: it also cares to normalize Xv!
	// Computes a proper change of base matrix from local to world
	dmat3 T(Xv.x, Xv.y, Xv.z,
			Yv.x, Yv.y, Zv.z,
			Zv.x, Zv.y, Zv.z);
	// Invert the matrix to obtain a change of base matrix from world to local
	T = T.Inverse();
	// Express V2 in the local reference frame of V1.
	// Let's call the resulting vector 'W'
	dvec3 W = T * V2;

	/* Error computation */

	double v1_len = V1.Length();

	// Params of the elliptic similarity function and their squared value
	double a = 0.5  * v1_len; // Radius of the ellipse along the V1 direction (Xv)
	double b = 0.25 * v1_len; // Radius of the ellipse along the a direction (Yv) orthogonal to Xv
	double c = 0.25 * v1_len; // Radius of the ellipse along the a direction (Zv) orthogonal to Xv and Yv
	double d = 0.25 * v1_len; // Offset detween the apex of Xv and the center of the ellipse

	double a2 = a * a;
	double b2 = b * b;
	double c2 = c * c;
	double d2 = d * d;

	double k1 = W.x - v1_len;
	double k2 = a2 - d2;

	double p1 = (a2 * k1 * k1) + (k2 * ((a2*W.y*W.y/b2) + (a2*W.z*W.z/c2)));
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
