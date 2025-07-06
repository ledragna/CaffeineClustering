#ifndef VECTORCLUSTERER_H
#define VECTORCLUSTERER_H

#include <memory>
#include <unordered_set>

#include <SciData/vectorfield3d.h>
#include <Utilities/mathfu_utilities.h>

/**
 * @author Andrea Salvadori and Marco Fusè
 */
namespace SNS { namespace SciData
{

// TODO: doc
class VectorClusterer
{
public:
	VectorClusterer();

private:

	struct GridCoordsHasher
	{
		size_t operator()(mathfu::Vector<uint,3> p) const
		{
			// Computes an esh table for the voxel's coordinates
			// using the "prime number multiplication" algorithm.
			// See: https://commons.apache.org/proper/commons-lang/apidocs/org/apache/commons/lang3/builder/HashCodeBuilder.html
			// See: http://myeyesareblind.com/2017/02/06/Combine-hash-values/

			const size_t k1 = 13; // must be a prime number
			const size_t k2 = 71; // must be a prime number

			size_t hx = std::hash<uint>()(p.x);
			size_t hy = std::hash<uint>()(p.y);
			size_t hz = std::hash<uint>()(p.z);

			size_t hash = k1;
			hash = (hash * k2) + hx;
			hash = (hash * k2) + hy;
			hash = (hash * k2) + hz;

			return hash;
		}
	};

	struct GridCoordsComparator
	{
		bool operator()(const mathfu::Vector<uint,3>& p1,
						const mathfu::Vector<uint,3>& p2) const
		{
			return (p1.x == p2.x) && (p1.y == p2.y) && (p1.z == p2.z);
		}
	};

	struct Cluster
	{
		mathfu::Vector<double,3> Centroid;
		mathfu::Vector<double,3> Vector;
		double Volume;

		std::unordered_set< mathfu::Vector<uint,3>,
							GridCoordsHasher,
							GridCoordsComparator > Voxels;
		std::unordered_set< mathfu::Vector<uint,3>,
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

	double computeDirMagError(const mathfu::Vector<double,3>& V1,
							  const mathfu::Vector<double,3>& V2);

	// Returns 0 when P1 == P2 and "err_max" when len(P2-P1) = max_distance.
	// Otherwise returns the linear iterpolation between these two values.
	double computePosError_Version1(const mathfu::Vector<double,3>& P1,
									const mathfu::Vector<double,3>& P2,
									double max_distance, double err_max);

	// Returns 1 when P1 == P2 and "err_max" when len(P2-P1) = max_distance.
	// Otherwise returns the linear iterpolation between these two values.
	double computePosError_Version2(const mathfu::Vector<double,3>& P1,
									const mathfu::Vector<double,3>& P2,
									double max_distance, double err_max);

	double computeError_Version1(const mathfu::Vector<double,3>& V1,
								 const mathfu::Vector<double,3>& P1,
								 const mathfu::Vector<double,3>& V2,
								 const mathfu::Vector<double,3>& P2,
								 double max_distance, double pos_err_max,
								 double weight);

	double computeError_Version2(const mathfu::Vector<double,3>& V1,
								 const mathfu::Vector<double,3>& P1,
								 const mathfu::Vector<double,3>& V2,
								 const mathfu::Vector<double,3>& P2,
								 double max_distance, double pos_err_max);

	ClusterSP performClustering(ConstVectorField3DSP grid);

};

} }

#endif // VECTORCLUSTERER_H
