#ifndef GRIDCLUSTERER_H
#define GRIDCLUSTERER_H

#include <atomic>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <vector>
#include <functional>
#include <map>

#include <SciData/vectorfield3d.h>
#include <Utilities/array3d.h>
#include <Utilities/iterable_priority_queue.h>
#include <Utilities/index3d_hashing.h>
#include <Utilities/boost_math_utilities.h>
#include <Utilities/log.h>

/**
 * @author Andrea Salvadori and Marco Fusè
 */
namespace CFF { namespace SciData
{
	//TODO: doc
	// T must support operator+(T,T) , operator*(double,T) and default constructor T()
	template<typename T>
	class GridClusterer
	{

	private:
		/*** PRIVATE DATA STRUCTURES ***/

		// Used to ask the termination of the performClustering()
		// procedure from a different thread.
		// See the interruptClustering() method.
		std::atomic_flag _can_continue = ATOMIC_FLAG_INIT;


	public:
		/*** PUBLIC DATA STRUCTURES ***/

		struct Cluster
		{
			friend class GridClusterer;
		public:
			boost::numeric::ublas::vector<double> Centroid;
			T Value;

			/*** FOR INTERNAL AND DEBUG USE ONLY ***/
			bool Clustered;
			size_t ErrorLevel;
			size_t Volume; // MUST BE > 0
			std::shared_ptr<Cluster> RightChild;
			std::shared_ptr<Cluster> LeftChild;
			std::unordered_set<Cluster*> Neighbours;
			/**************************************/			/* Default constructor */
			Cluster()
				: Centroid(3), Value(T()),
				  Clustered(false), ErrorLevel(0), Volume(1)
			{
				Centroid(0) = 0.0;
				Centroid(1) = 0.0;
				Centroid(2) = 0.0;
			}

			//~Cluster() { qDebug() << "~Cluster"; }
		};


	private:
		/*** PRIVATE DATA STRUCTURES ***/

		struct MergingCandidate
		{
			Cluster* first;
			Cluster* second;
			double distance;

			MergingCandidate(Cluster* a, Cluster* b, double dist)
				: first(a), second(b), distance(dist) {}

			MergingCandidate(const MergingCandidate& other) = default;
			MergingCandidate& operator=(const MergingCandidate& other) = default;

			MergingCandidate(MergingCandidate&& other)
				: first(other.first), second(other.second), distance(other.distance)
			{
				other.first = NULL;
				other.second = NULL;
				other.distance = 0.0;
			}

			MergingCandidate& operator=(MergingCandidate&& other)
			{
				first = other.first;
				second = other.second;
				distance = other.distance;

				other.first = NULL;
				other.second = NULL;
				other.distance = 0.0;

				return *this;
			}

			bool operator<(const MergingCandidate& other) const
			{ return this->distance < other.distance; }

			bool operator>(const MergingCandidate& other) const
			{ return this->distance > other.distance; }
		};


	public:
		/*** INSTANCE METHODS ***/

		/**
		 * @brief	Allows to interrupt the execution of the performClustering()
		 *			procedure from another thread.
		 */
		void interruptClustering()
		{ _can_continue.clear(); }

		// DistFunctor must implement the following constructor and operator():
		// ClassName(const DataGrid3D<T>* grid) // constructor
		// double operator()(const Cluster*, const Cluster*) const
		//
		// Validator must implement the following constructor and operator():
		// ClassName(const DataGrid3D<T>* grid) // constructor
		// bool operator()(size_t i, size_t j, size_t k, const T& value) const
		//
		// If the the shared pointer passed as parameter is the only reference
		// to the grid, that is if the grid is going to be destroyed, interrupts
		// the clustering procedure and returns a warning message explaining the
		// reasons (completing the procedure is useless...).
		//
		// Special case: The dataset does not contains any valid elements
		// -> returns NULL and a warning message.
		template<typename DistFunctor, typename Validator>
		std::shared_ptr<Cluster> performClustering( std::shared_ptr<const DataGrid3D<T>> grid,
													std::string& warning,
													bool freeUnusedNeighbours = true)
		{
			typedef boost::numeric::ublas::vector<double> dvec3;
			typedef boost::numeric::ublas::vector<unsigned> uvec3;
			typedef std::shared_ptr<GridClusterer<T>::Cluster> ClusterSP;

			if(!grid) Utilities::throwAndPrintError<std::invalid_argument>("GridClusterer::performClustering - Null grid!");
			const DataGrid3D<T>* gridPtr = grid.get();

			// Initialize the termination flag
			_can_continue.test_and_set();
			// Clears the warning message
			warning.clear();

			// Dictionary associating each cluster's raw pointer to its shared_ptr.
			// Needed to keep cluster objects in life until the end of the procedure.
			std::unordered_map<Cluster*, ClusterSP> clustersSet;

			// 3D array used to store a cluster for each valid element of the grid.
			// If the element of the grid is not valid (i.e. it must to be excluded
			// from the clustering), the corresponding element of the array is set to NULL.
			Utilities::Array3D<Cluster*> leafNodes( grid->getNumPointsX(),
													grid->getNumPointsY(),
													grid->getNumPointsZ());

			// Function object used to test the elements of the grid for validity
			Validator validator(gridPtr);

			//qDebug() << "Creating a cluster for each voxel of the grid";

			// Creates a cluster for each voxel of the grid
			for(uint i=0; i < grid->getNumPointsX(); ++i)
				for(uint j=0; j < grid->getNumPointsY(); ++j)
					for(uint k=0; k < grid->getNumPointsZ(); ++k)
					{
						// Checks if the procedure must be interrupted
						if(!_can_continue.test_and_set())
						{
							warning = "Clustering procedure interrupted correctly :-)";
							return NULL;
						}

						uvec3 coords(3);
						coords(0) = i;
						coords(1) = j;
						coords(2) = k;
						T value = grid->getValue(coords);
						if(validator(i,j,k,value))
						{							ClusterSP C = std::make_shared<Cluster>();
							C->Centroid = dvec3(3);
							C->Centroid(0) = i;
							C->Centroid(1) = j;
							C->Centroid(2) = k;
							C->Value = value;
							C->Volume = 1;
							C->ErrorLevel = 0;
							C->Clustered = false;
							// Note: neighbours will be initialized later

							leafNodes(i,j,k) = C.get();
							clustersSet.insert({C.get(), C});
						}
						else leafNodes(i,j,k) = NULL;
					}

			/* Special case: The dataset does not contains any valid elements -> returns NULL */
			if(clustersSet.empty())
			{
				warning = "The clustering can't take place since the dataset "
						  "does not contains any valid item according to the "
						  "specified validator.";
				return NULL;
			}
			/***********************************************************/
			/* Special case: The dataset contains only 1 valid element */
			if(clustersSet.size() == 1) return clustersSet.begin()->second;
			/***********************************************************/

			// Ordered queue. Each item is a pair of Clusters to be merged
			// plus their distance. The "top" method returns the candidate
			// with the LOWER distance.
			Utilities::MyPriorityQueue<	MergingCandidate,
										std::vector<MergingCandidate>,
										std::greater<MergingCandidate>> mergingCandidatesQueue;

			// Keeps tracks of the lists nodes that, due to the presence
			// of invalid elements in the dataset, are isolated (i.e. have
			// no neighbours).
			// The booleans are used in the following neighburns search process
			// to verify if it is performed on each isolatedleaf
			std::vector<std::pair<bool,Cluster*>> isolatedLeaves;

			//qDebug() << "Initializing the neighbours of each leaf node";

			// Initializes the neighbours of each leaf node
			DistFunctor distF(gridPtr);
			for(auto iter : clustersSet)
			{
				// Checks if the procedure must be interrupted
				if(!_can_continue.test_and_set())
				{
					warning = "Clustering procedure interrupted correctly :-)";
					return NULL;
				}				Cluster* C = iter.first;
				// Coordinates of the cluster
				uint i = static_cast<uint>(std::round(C->Centroid(0)));
				uint j = static_cast<uint>(std::round(C->Centroid(1)));
				uint k = static_cast<uint>(std::round(C->Centroid(2)));

				if( (i > 0) && (leafNodes(i-1,j,k) != NULL) )
					C->Neighbours.insert(leafNodes(i-1,j,k));

				if( (j > 0) && (leafNodes(i,j-1,k) != NULL) )
					C->Neighbours.insert(leafNodes(i,j-1,k));

				if( (k > 0) && (leafNodes(i,j,k-1) != NULL) )
					C->Neighbours.insert(leafNodes(i,j,k-1));

				// Note: to avoid duplicates, only candidates with
				// indices >= (i,j,k) are inserted to the queue.

				if( (i < grid->getNumPointsX()-1) && (leafNodes(i+1,j,k) != NULL) )
				{
					Cluster* D = leafNodes(i+1,j,k);
					C->Neighbours.insert(D);
					MergingCandidate candidate(C, D, distF(C,D));
					mergingCandidatesQueue.push(candidate);
				}

				if( (j < grid->getNumPointsY()-1) && (leafNodes(i,j+1,k) != NULL) )
				{
					Cluster* F = leafNodes(i,j+1,k);
					C->Neighbours.insert(F);
					MergingCandidate candidate(C, F, distF(C,F));
					mergingCandidatesQueue.push(candidate);
				}

				if( (k < grid->getNumPointsZ()-1) && (leafNodes(i,j,k+1) != NULL) )
				{
					Cluster* G = leafNodes(i,j,k+1);
					C->Neighbours.insert(G);
					MergingCandidate candidate(C, G, distF(C,G));
					mergingCandidatesQueue.push(candidate);
				}

				/*
				 * Special case: due to invalid elements in the dataset
				 * (according to the Validator), the leaf node is isolated
				 * (i.e. it has no neighbours).
				 * That will cause this object to be discarted.
				 * To avoid that, we must consider it as neighbour of any other
				 * leaf node (see later).
				 */
				if(C->Neighbours.empty())
				{
//					qDebug() << C->Centroid.x << C->Centroid.y << C->Centroid.z << "is isolated.";
					isolatedLeaves.push_back({false,C});
				}
			}


			/*
			 * Special case: due to invalid elements in the dataset
			 * (according to the Validator), some leaf node gets isolated
			 * (i.e. it has no neighbours).
			 * That will cause this object to be discarted.
			 * To avoid that, a neighbourns search is performed.
			 * The process searches for not NULL neighbourns in cubic surfaces
			 * of increasing size
			 */
			for(auto isolatedIter : isolatedLeaves)
			{				Cluster* isolatedCluster = isolatedIter.second;
				uint i = static_cast<uint>(std::round(isolatedCluster->Centroid(0)));
				uint j = static_cast<uint>(std::round(isolatedCluster->Centroid(1)));
				uint k = static_cast<uint>(std::round(isolatedCluster->Centroid(2)));
				//Variable to go in concentric cubes
				int radius = 0;
				int step = 1;
				while(!isolatedIter.first)
				{
					radius++;
					for(int v=-radius; v <= radius; v++)
					{
						if((i+v < 0) || (i+v > grid->getNumPointsX()-1) )
							continue;
						for(int w=-radius; w <= radius; w++)
						{
							if((j+w < 0) || (j+w > grid->getNumPointsY()-1) )
								continue;
							// to avoid the complete sampling if abs(v) and
							// abs(w) are different to radius the z leaps
							// between the edges
							if((abs(v) != radius) && (abs(w) != radius))
								step = radius * 2;
							else
								step = 1;
							for(int z=-radius; z <= radius; z+=step)
							{
								if((k+z < 0) || (k+z > grid->getNumPointsZ()-1) )
									continue;
								if (leafNodes(i+v,j+w,k+z) != NULL)
								{
									Cluster* G = leafNodes(i+v,j+w,k+z);
									isolatedCluster->Neighbours.insert(G);
									G->Neighbours.insert(isolatedCluster);
									MergingCandidate candidate(isolatedCluster, G, distF(isolatedCluster,G));
									mergingCandidatesQueue.push(candidate);
									isolatedIter.first = true;
								}
							}
						}
					}
				}

			}

			// Frees some memory
			leafNodes.clear();
			isolatedLeaves.clear();

			//qDebug() << "Iteratively merging the closest clusters";

			// Iteratively merges the closest clusters until
			// mergingCandidatesQueue is not empty.
			ClusterSP root(NULL);
			size_t errLevel = 0;
			while(!mergingCandidatesQueue.empty())
			{
				//qDebug() << "mergingCandidatesQueue size:" << mergingCandidatesQueue.size();

				// Checks if the procedure must be interrupted
				if(!_can_continue.test_and_set())
				{
					warning = "Clustering procedure interrupted correctly :-)";
					return NULL;
				}

				MergingCandidate candidate = mergingCandidatesQueue.top(); // Candidate with min dist
				mergingCandidatesQueue.pop();

				if( candidate.first->Clustered || candidate.second->Clustered ) continue;

				errLevel++;
				ClusterSP candidateA = clustersSet.at(candidate.first);
				ClusterSP candidateB = clustersSet.at(candidate.second);
				ClusterSP newCluster = mergeClusters(candidateA, candidateB,
													 errLevel, freeUnusedNeighbours);
				clustersSet.insert({newCluster.get(), newCluster});
				root = newCluster;

				/*
				 * Special case:
				 * due to invalid elements in the dataset (according to the
				 * Validator), the new cluster is not adjacent to any other
				 * cluster (i.e. it has no neighbours).
				 * That will cause this object to be discarted.
				 * To avoid that, we must consider it as neighbour of any other
				 * cluster (still unclustered).
				 */
				if(newCluster->Neighbours.empty())
				{
					//Number of neighbourns to add
					uint MAP_SIZE = 10u;
					//the map contains the first MAP_SIZE closest neighbourns
					//multimap is used to avoid the discard of cluster with the same euclidean distance
					//the supportSet contains the actual pointer to the cluster, once the search is finished
					//supportSet is moved to the Neighbourns in newCluster
					std::multimap<double, Cluster*> neighSearch;
					std::unordered_set<Cluster*> supportSet;

//					qDebug() << newCluster->Centroid.x << newCluster->Centroid.y << newCluster->Centroid.z << "is isolated.";
					for(const MergingCandidate& aMergingCandidate : mergingCandidatesQueue)
					{
						if(!aMergingCandidate.first->Clustered)
						{
							double tmp = boost::numeric::ublas::norm_2(aMergingCandidate.first->Centroid - newCluster->Centroid);
							if (neighSearch.size() <= MAP_SIZE)
							{
								supportSet.insert(aMergingCandidate.first);
								neighSearch.insert({tmp,aMergingCandidate.first});
							}
							else if ((tmp < std::prev(neighSearch.end())->first)
								&& (supportSet.find(aMergingCandidate.first) != supportSet.end()))
							{
								supportSet.erase(std::prev(neighSearch.end())->second);
								neighSearch.erase(std::prev(neighSearch.end()));
								neighSearch.insert({tmp,aMergingCandidate.first});
							}

//							qDebug() << newCluster->Centroid.x << newCluster->Centroid.y << newCluster->Centroid.z << "and"
//									 << aMergingCandidate.first->Centroid.x << aMergingCandidate.first->Centroid.y << aMergingCandidate.first->Centroid.z << "are now neighbours.";
						}

						if(!aMergingCandidate.second->Clustered)
						{
//							qDebug() << newCluster->Centroid.x << newCluster->Centroid.y << newCluster->Centroid.z << "and"
//									 << aMergingCandidate.second->Centroid.x << aMergingCandidate.second->Centroid.y << aMergingCandidate.second->Centroid.z << "are now neighbours.";
							double tmp = boost::numeric::ublas::norm_2(aMergingCandidate.second->Centroid - newCluster->Centroid);
							if (neighSearch.size() <= MAP_SIZE)
							{
								supportSet.insert(aMergingCandidate.second);
								neighSearch.insert({tmp,aMergingCandidate.second});
							}
							else if ((tmp < std::prev(neighSearch.end())->first)
								&& (supportSet.find(aMergingCandidate.second) != supportSet.end()))
							{
								supportSet.erase(std::prev(neighSearch.end())->second);
								neighSearch.erase(std::prev(neighSearch.end()));
								neighSearch.insert({tmp,aMergingCandidate.second});
							}

						}
					}
					neighSearch.clear();
					newCluster->Neighbours = std::move(supportSet);
//					for(Cluster* isolatedCluster = supportSet.begin(); isolatedCluster != supportSet.end(); ++isolatedCluster )
//					{
//						newCluster->Neighbours.insert();
//						//aMergingCandidate.first->Neighbours.insert(newCluster.get());
//					}
				}

				for(Cluster* neighbour : newCluster->Neighbours)
				{
					MergingCandidate candidate(newCluster.get(), neighbour,
											   distF(newCluster.get(),neighbour));
					mergingCandidatesQueue.push(candidate);
				}
			}

			// Consistency check
			if(root)
			{
				size_t numNodesInTree = (2 * root->ErrorLevel) + 1;
				if(numNodesInTree != clustersSet.size())
				{
					warning = "BUG: GridClusterer::performClustering : the "
							  "number of clusters in the tree is different "
							  "from the total number of generated numbers. "
							  "Some clusters will be erroneously discarded! "
							  "Please report this notification!";
				}
			}

			return root;
		}



		// FOR INTERNAL AND DEBUG USE ONLY
		std::shared_ptr<Cluster> mergeClusters(std::shared_ptr<Cluster> c1,
											   std::shared_ptr<Cluster> c2,
											   size_t level,
											   bool freeUnusedNeighbours = true) const
		{
			std::shared_ptr<Cluster> newCluster = std::make_shared<Cluster>();
			newCluster->Volume = c1->Volume + c2->Volume;
			double w1 = double(c1->Volume) / double(newCluster->Volume);
			double w2 = double(c2->Volume) / double(newCluster->Volume);
			newCluster->Centroid = (w1 * c1->Centroid) + (w2 * c2->Centroid);
			newCluster->Value =    (w1 * c1->Value)    + (w2 * c2->Value);

			if(c1->ErrorLevel <= c2->ErrorLevel)
			{
				newCluster->LeftChild = c1;
				newCluster->RightChild = c2;
			}
			else
			{
				newCluster->LeftChild = c2;
				newCluster->RightChild = c1;
			}

			newCluster->ErrorLevel = level;
			newCluster->Clustered = false;

			c1->Clustered = true;
			c2->Clustered = true;

			// The set of neighbours of the new cluster is given by the union of the
			// neighbours of the two merged clusters (c1 and c2) minus themselves.
			newCluster->Neighbours.insert(c1->Neighbours.begin(), c1->Neighbours.end());
			newCluster->Neighbours.insert(c2->Neighbours.begin(), c2->Neighbours.end());
			newCluster->Neighbours.erase(c1.get());
			newCluster->Neighbours.erase(c2.get());

			// To save some memory
			if(freeUnusedNeighbours)
			{
				c1->Neighbours.clear();
				c2->Neighbours.clear();
			}

			// For each adjacent cluster, update its set of neighbours by removing the two
			// merged clusters (c1 and c2) from the set and adding the new cluster to it.
			for(Cluster* n : newCluster->Neighbours)
			{
				n->Neighbours.erase(c1.get());
				n->Neighbours.erase(c2.get());
				n->Neighbours.insert(newCluster.get());
			}

			return newCluster;
		}


public:
		/*** STATIC METHODS ***/

		static std::vector<std::shared_ptr<const Cluster>> extractClusters(
									std::shared_ptr<Cluster> clusteringTreeRoot,
									size_t num_clusters)
		{
			if(!clusteringTreeRoot)
				Utilities::throwAndPrintError<std::invalid_argument>("GridClusterer::extractClusters - Null clustering tree!");

			if(num_clusters == 0) return {};

			std::vector<std::shared_ptr<const Cluster>> resultSet;
			size_t numLeaves = clusteringTreeRoot->Volume;
			num_clusters = std::min(numLeaves, num_clusters);
			size_t targetLevel = numLeaves - num_clusters;

			/* Iterative depth first tree traversal */
			std::stack<std::shared_ptr<Cluster>> stack;
			std::shared_ptr<Cluster> current = clusteringTreeRoot;

			while(current)
			{
				if(current->ErrorLevel <= targetLevel)
				{
					resultSet.push_back(current);
					current.reset();
				}
				else
				{
					if(current->RightChild) stack.push(current->RightChild);
					current = current->LeftChild;
				}

				if( (!current) && (!stack.empty()) )
				{
					current = stack.top();
					stack.pop();
				}
			}

//			qDebug() << "Extracted Clusters:" << num_clusters;
//			for(const std::shared_ptr<Cluster>& c : resultSet)
//				qDebug() << "C:" << c->Centroid.x << c->Centroid.y << c->Centroid.z
//						 << "; V:" << c->Value.x << c->Value.y << c->Value.z ;

			return resultSet;
		}
	};

} }

#endif // GRIDCLUSTERER_H
