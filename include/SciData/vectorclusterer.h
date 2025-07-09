#ifndef VECTORCLUSTERER_H
#define VECTORCLUSTERER_H

#include <memory>
#include <unordered_set>

#include <SciData/vectorfield3d.h>
#include <Utilities/boost_math_utilities.h>

/**
 * @author Andrea Salvadori and Marco Fusè
 */
namespace CFF { namespace SciData
{

// TODO: doc
class VectorClusterer
{
public:
	VectorClusterer();

private:

	struct GridCoordsHasher
	{		size_t operator()(CFF::Utilities::Vector3u p) const
		{
			// Computes an esh table for the voxel's coordinates
			// using the "prime number multiplication" algorithm.
			// See: https://commons.apache.org/proper/commons-lang/apidocs/org/apache/commons/lang3/builder/HashCodeBuilder.html
			// See: http://myeyesareblind.com/2017/02/06/Combine-hash-values/
			const size_t k1 = 13; // must be a prime number
			const size_t k2 = 71; // must be a prime number
			
			size_t hx = std::hash<uint>()(p(0));
			size_t hy = std::hash<uint>()(p(1));
			size_t hz = std::hash<uint>()(p(2));

			size_t hash = k1;
			hash = (hash * k2) + hx;
			hash = (hash * k2) + hy;
			hash = (hash * k2) + hz;

			return hash;
		}
	};
	struct GridCoordsComparator
	{
		bool operator()(const CFF::Utilities::Vector3u& p1,
						const CFF::Utilities::Vector3u& p2) const
		{
			return (p1(0) == p2(0)) && (p1(1) == p2(1)) && (p1(2) == p2(2));
		}
	};
	struct Cluster
	{
		CFF::Utilities::Vector3d Centroid;
		CFF::Utilities::Vector3d Vector;
		double Volume;
		std::unordered_set< CFF::Utilities::Vector3u,
							GridCoordsHasher,
							GridCoordsComparator > Voxels;
		std::unordered_set< CFF::Utilities::Vector3u,
							GridCoordsHasher,
							GridCoordsComparator > Neighbours;
		size_t Level;
		std::shared_ptr<Cluster> RightChild;
		std::shared_ptr<Cluster> LeftChild;
		bool Clustered;

	};

	typedef std::shared_ptr<Cluster> ClusterSP;

private:

	ClusterSP mergeClusters(ClusterSP c1, ClusterSP c2, size_t level);
	double computeDirMagError(const CFF::Utilities::Vector3d& V1,
							  const CFF::Utilities::Vector3d& V2);

	// Returns 0 when P1 == P2 and "err_max" when len(P2-P1) = max_distance.
	// Otherwise returns the linear iterpolation between these two values.
	double computePosError_Version1(const CFF::Utilities::Vector3d& P1,
									const CFF::Utilities::Vector3d& P2,
									double max_distance, double err_max);
	// Returns 1 when P1 == P2 and "err_max" when len(P2-P1) = max_distance.
	// Otherwise returns the linear iterpolation between these two values.
	double computePosError_Version2(const CFF::Utilities::Vector3d& P1,
									const CFF::Utilities::Vector3d& P2,
									double max_distance, double err_max);

	double computeError_Version1(const CFF::Utilities::Vector3d& V1,
								 const CFF::Utilities::Vector3d& P1,
								 const CFF::Utilities::Vector3d& V2,
								 const CFF::Utilities::Vector3d& P2,
								 double max_distance, double pos_err_max,
								 double weight);

	double computeError_Version2(const CFF::Utilities::Vector3d& V1,
								 const CFF::Utilities::Vector3d& P1,
								 const CFF::Utilities::Vector3d& V2,
								 const CFF::Utilities::Vector3d& P2,
								 double max_distance, double pos_err_max);

	ClusterSP performClustering(ConstVectorField3DSP grid);

};

} }

#endif // VECTORCLUSTERER_H
