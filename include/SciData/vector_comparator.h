#ifndef VECTORCOMPARATOR_H
#define VECTORCOMPARATOR_H

#include <SciData/grid_clusterer.h>
#include <Utilities/boost_math_utilities.h>

/**
 * @author Andrea Salvadori and Marco Fusè
 */
namespace CFF { namespace SciData
{

	/**
	 * @brief	Functor designed to be used with CFF::SciData::GridClusterer,
	 *			and in particular as "DistFunctor" template of the
	 *			performClustering() method.
	 */
	class VectorComparator
	{
	private:
		double _maxDist;
		double _maxMagnitude;

	public:

		/**
		 * @brief Constructor. Requires the analyzed vector field as parameter.
		 */
		VectorComparator(const VectorField3D* vecField);

		/**
		 * @brief	Returns a grade of dissimilarity between the two clusters
		 *			of vectors passed as parameter. The more the returned value
		 *			is high, the more the two cluster are dissimilar.
		 */
		double operator()(const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster1,
						  const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster2) const;

	private:

		double magnitudeComparator(const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster1,
								   const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster2) const;
		double distanceComparator(const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster1,
								  const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster2) const;
		double directionComparator(const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster1,
								   const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster2) const;
		double euclideanDirMagComparison(const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster1,
										 const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster2) const;
		double teleaDirMagComparison(const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster1,
									 const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster2) const;
		double teleaDirMagAsymmetric(const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster1,
									 const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster2) const;
		double volumeComparison(const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster1,
								const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster2) const;
		bool zeroValue(const CFF::SciData::GridClusterer<CFF::Utilities::Vector3d>::Cluster* cluster1) const;

	};

} }

#endif // VECTORCOMPARATOR_H
