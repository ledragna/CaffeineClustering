#ifndef DATAGRID3D_H
#define DATAGRID3D_H

#include <Utilities/boost_math_utilities.h>
#include <Utilities/math3d.h>
#include <Utilities/log.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <memory>
#include <array>
#include <vector>
#include <stdexcept>

// Type aliases for easier migration
using Vector3d = SNS::Utilities::Vector3d;
using Vector4d = SNS::Utilities::Vector4d;
using Matrix4d = SNS::Utilities::Matrix4d;
using Vector3u = SNS::Utilities::Vector3u;

/**
 * @author Andrea Salvadori
 */
namespace SNS { namespace SciData
{
	/**
	 *	Instances of this class represent volumetric datasets
	 *	whose voxels are equally distanced in space.
	 *
	 *	A volumetric dataset is a set of pairs <Pi, Vi> called “voxels”
	 *	(short for “volume elements”), where Pi is a point in space and Vi
	 *	is its associated value. In DataGrid3D, the location of the voxels
	 *  are equally spaced locations in space, so to obtain a regular
	 *	three-dimensional grid of values. Note that the number of voxels
	 *  can differ among the three dimensions.
	 *
	 *	The value for an arbitrary point lying within a cell of the grid
	 *  can be approximated by interpolating the values of the eight voxels
	 *	delimiting the cell (see the sampleValue() method of this class).
	 *
	 *  The operations provided by this class that need a point of the grid
	 *	as parameter, require the coordinates of the point to be expressed in
	 *  the local reference frame of the grid. Let <l,n,m> be the number of
	 *	voxels along the local <X,Y,Z> axis, the range of the local coordinates
	 *	is [0,l-1] for the X component, [0,n-1] for the Y component and [0,m-1]
	 *	for the Z component. The overloaded method setLocalToWorldTransform()
	 *	allows to set a change of reference frame matrixfrom the local
	 *	coordinate system to the global (world) one.
	 *	The method getLocalToWorldTransform() returns the local to world matrix,
	 *	while getWorldToLocalTransform() returns the world to local matrix.
	 *	Furthermore, the methods fromLocalToWorldFrame() and
	 *	fromWorldToLocalFrame() are provided to convert a point from local to
	 *	world coordinates and vice versa.
	 *
	 *  Finally, there are some requirements for the type parameter T:
	 *  - Must provide a default constructor.
	 *  - Must provide a copy constructor.
	 *  - Must provide an assignement operator (=).
	 *
	 *  Additional requirements are possible, according to the methods you
	 *	plan to use. Check the documentation of the methods you need for
	 *	further details.
	 */
	template <typename T> class DataGrid3D
	{
	protected:

		/**
		 * @brief A descriptive name for the dataset.
		 */
		std::string _name;

		/**
		 * @brief Number of voxels (vertices) of the grid along the 3 dimensions
		 */
		std::array<uint,3> _nPoints;

		/**
		 * @brief Change of reference frame matrix from local coordinates
		 *		  of the grid (from 0 to _nPoints[i]-1 for each dimension)
		 *		  to world (global) coordinates.
		 */
		Matrix4d _localToWorld;

		/**
		 * @brief Change of reference frame matrix from world (global)
		 *		  coordinates to the local reference frame of the grid.
		 */
		Matrix4d _worldToLocal;

		/**
		 * @brief 3D matrix of values
		 */
		std::vector< std::vector< std::vector<T>>> _values;

	protected:

		/** @brief  Relative position of the 8 corners of the cube */
		const unsigned cubeCornersRelativePos[8][3] = {
			{0u,0u,0u}, {0u,0u,1u}, {0u,1u,0u}, {0u,1u,1u},
			{1u,0u,0u}, {1u,0u,1u}, {1u,1u,0u}, {1u,1u,1u},
		};

		/**
		 * @brief	The Cube struct represents a region of space bounded by 8
		 *			adjacent voxels of the 3D grid. Each instance of this
		 *			structure stores the local coordinates of its 8 voxels
		 *			and their associated value, in the following order:
		 *			{V000, V001, V010, V011, V100, V101, V110, V111}
		 */
		struct Cube
		{
			/** @brief Local coordinates of the 8 vertices of the cube within the grid */
			std::array<std::array<uint,3>,8> voxelsCoords;

			/** @brief Values associated to the 8 voxels */
			std::array<T,8> voxelsValues;
		};

		/**
		 * @brief	Helper method that, given the LOCAL coordinates of a point
		 *			within the grid (localCoords), returns a Cube structure
		 *			(outCube) containing the the LOCAL coordinates of the eight
		 *			vertices of the cube (bounding the specified point) and
		 *			their associated values. Furthermore, this method also
		 *			returns the noormalized coordinates (in the range [0,1])
		 *			of the specified point within the cube (outNormCoords).
		 */
		void getCubeInfo(	const Vector3d& localCoords,
							Cube& outCube, Vector3d& outNormCoords) const
		{
			unsigned cubeOriginX = unsigned(localCoords(0));
            if(cubeOriginX == (_nPoints[0]-1u))
			{
				cubeOriginX -= 1u;
				outNormCoords[0] = 1.0;
			}
			else
				outNormCoords[0] = localCoords(0) - cubeOriginX;

			unsigned cubeOriginY = unsigned(localCoords(1));
            if(cubeOriginY == (_nPoints[1]-1u))
			{
				cubeOriginY -= 1u;
				outNormCoords[1] = 1.0;
			}
			else
				outNormCoords[1] = localCoords(1) - cubeOriginY;

			unsigned cubeOriginZ = unsigned(localCoords(2));
            if(cubeOriginZ == (_nPoints[2]-1u))
			{
				cubeOriginZ -= 1u;
				outNormCoords[2] = 1.0;
			}
			else
				outNormCoords[2] = localCoords(2) - cubeOriginZ;

			for(unsigned i = 0; i < 8; i++)
			{
				uint voxelCoordX = cubeOriginX + cubeCornersRelativePos[i][0];
				uint voxelCoordY = cubeOriginY + cubeCornersRelativePos[i][1];
				uint voxelCoordZ = cubeOriginZ + cubeCornersRelativePos[i][2];

				outCube.voxelsCoords[i][0]= voxelCoordX;
				outCube.voxelsCoords[i][1]= voxelCoordY;
				outCube.voxelsCoords[i][2]= voxelCoordZ;
				outCube.voxelsValues[i] = _values[voxelCoordX][voxelCoordY][voxelCoordZ];
			}
		}

	public:

		/**
		 * @brief Constructor.
		 * @param name		A descriptive name for the dataset.
		 * @param nPointsX	Number of voxels along X. Must be > 1.
		 * @param nPointsY	Number of voxels along Y. Must be > 1.
		 * @param nPointsZ	Number of voxels along Z. Must be > 1.
		 */
		DataGrid3D(	const std::string& name,
					unsigned nPointsX, unsigned nPointsY, unsigned nPointsZ)
		{
			if((nPointsX < 2) || (nPointsY < 2) || (nPointsZ < 2))
				Utilities::throwAndPrintError<std::invalid_argument>(
								"DataGrid3D::DataGrid3D - Invalid grid size!");

			_name = name;
			_nPoints[0] = nPointsX;
			_nPoints[1] = nPointsY;
			_nPoints[2] = nPointsZ;
			_localToWorld = SNS::Utilities::identityMatrix4d();
			_worldToLocal = SNS::Utilities::identityMatrix4d();

			_values.resize(nPointsX);
			for(unsigned x = 0; x < nPointsX; x++)
			{
				_values[x].resize(nPointsY);
				for(unsigned y = 0; y < nPointsY; y++)
				{
					_values[x][y].resize(nPointsZ);
					for(unsigned z = 0; z < nPointsZ; z++)
					{
						_values[x][y][z] = T();
					}
				}
			}
		}

		/** @brief Destructor. */
		virtual ~DataGrid3D() = default;

		/**
		 * @brief Copy constructor.
		 */
		DataGrid3D(const DataGrid3D<T>& other)
			: _name(other._name),
			  _nPoints{other._nPoints},
			  _localToWorld(other._localToWorld),
			  _worldToLocal(other._worldToLocal),
			  _values(other._values)
		{}

		/**
		 * @brief Move constructor.
		 */
		DataGrid3D(DataGrid3D<T>&& other)
			: _name(std::move(other._name)),
			  _nPoints{std::move(other._nPoints)},
			  _localToWorld(std::move(other._localToWorld)),
			  _worldToLocal(std::move(other._worldToLocal)),
			  _values(std::move(other._values))
		{
			// Clears the source
			other._name.clear();
			other._nPoints = {0,0,0};
			other._localToWorld = SNS::Utilities::identityMatrix4d();
			other._worldToLocal = SNS::Utilities::identityMatrix4d();
			other._values.clear();
		}

		/**
		 * @brief Copy assignement operator.
		 */
		virtual DataGrid3D<T>& operator=(const DataGrid3D<T>& other)
		{
			if(this != &other)
			{
				_name = other._name;
				_nPoints = other._nPoints;
				_localToWorld = other._localToWorld;
				_worldToLocal = other._worldToLocal;
				_values = other._values;
			}

			return (*this);
		}

		/**
		 * @brief Move assignement operator.
		 */
		virtual DataGrid3D<T>& operator=(DataGrid3D<T>&& other)
		{
			if(this != &other)
			{
				// Steals the data
				_name = std::move(other._name);
				_nPoints = std::move(other._nPoints);
				_localToWorld = std::move(other._localToWorld);
				_worldToLocal = std::move(other._worldToLocal);
				_values = std::move(other._values);

				// Clears the source
				other._name.clear();
				other._nPoints = {0,0,0};
				other._localToWorld = SNS::Utilities::identityMatrix4d();
				other._worldToLocal = SNS::Utilities::identityMatrix4d();
				other._values.clear();
			}

			return (*this);
		}

		/**
		 * @brief   Factory method: Creates an heap allocated copy of
		 *          this object and returns a smart pointer to it.
		 */
		virtual std::shared_ptr<DataGrid3D<T>> clone() const
		{
			std::shared_ptr<DataGrid3D<T>> clone =
					std::make_shared<DataGrid3D<T>> (*this);
			return clone;
		}

		/**
		 * @brief Returns the name associated to this dataset.
		 */
		virtual std::string getName() const { return _name; }

		/**
		 * @brief Returns the number of voxels along the X axis.
		 */
		virtual unsigned getNumPointsX() const { return _nPoints[0]; }

		/**
		 * @brief Returns the number of voxels along the Y axis.
		 */
		virtual unsigned getNumPointsY() const { return _nPoints[1]; }

		/**
		 * @brief Returns the number of voxels along the Z axis.
		 */
		virtual unsigned getNumPointsZ() const { return _nPoints[2]; }


		/**
		 * @brief	Returns the origin of the local reference frame of the grid
		 *			(i.e. the point having [0,0,0] as local coordinates)
		 *			expressed in the global reference frame.
		 */		virtual Vector3d getOriginInWorldSpace() const
		{
			return SNS::Utilities::xyz(SNS::Utilities::getColumn(_localToWorld,3));
		}

		/**
		 * @brief	Returns the X axis of the local reference frame of the grid
		 *			expressed in the global reference frame.
		 */
		virtual Vector3d getLocalXAxisInWorldSpace() const
		{
			return SNS::Utilities::getColumn(_localToWorld,0);
		}

		/**
		 * @brief	Returns the Y axis of the local reference frame of the grid
		 *			expressed in the global reference frame.
		 */
		virtual Vector3d getLocalYAxisInWorldSpace() const
		{
			return SNS::Utilities::getColumn(_localToWorld,1);
		}

		/**
		 * @brief	Returns the Z axis of the local reference frame of the grid
		 *			expressed in the global reference frame.
		 */
		virtual Vector3d getLocalZAxisInWorldSpace() const
		{
			return SNS::Utilities::getColumn(_localToWorld,2);
		}

		/**
		 * @brief Returns the change of reference frame matrix from local
		 *		  coordinates of the grid to world (global) coordinates.
		 */
		virtual const Matrix4d& getLocalToWorldTransform() const
		{
			return _localToWorld;
		}

		/**
		 * @brief	Returns the change of reference frame matrix from
		 *			world (global) coordinates to local coordinates of the grid.
		 */
		virtual const Matrix4d& getWorldToLocalTransform() const
		{
			return _worldToLocal;
		}

		/**
		 * @brief Sets a new name for this dataset.
		 */
		virtual void setName(const std::string& newName) { _name = newName;}

		/**
		 * @brief	Sets the change of reference frame transformation from local
		 *			coordinates of the grid to world (global) coordinates.
		 *			Returns true if the procedure succeeded, false if it fails
		 *			since the matrix is not invertible.
		 *
		 * @param newLocalToWorld	The new change of reference frame matrix.
		 */		virtual bool setLocalToWorldTransform(const Matrix4d& newLocalToWorld)
		{
			Matrix4d newWorldToLocal;
			bool invertible = SNS::Utilities::InverseWithDeterminantCheck(newLocalToWorld, &newWorldToLocal);
			if(!invertible) return false;

			_localToWorld = newLocalToWorld;
			_worldToLocal = newWorldToLocal;
			return true;
		}

		/**
		 * @brief	Sets the change of reference frame transformation from local
		 *			coordinates of the grid to world (global) coordinates.
		 *			Returns true if the procedure succeeded, false if it fails
		 *			since the resulting matrix is not invertible.
		 *
		 * @param worldOrigin	Point in the global reference frame to set as
		 *						Origin of the local reference frame of the grid
		 *						(i.e. the point having [0,0,0] as local coords).
		 *
		 * @param localXInWorld Vector in the global reference frame to set as
		 *						X axis of the local reference frame of the grid.
		 *
		 * @param localYInWorld Vector in the global reference frame to set as
		 *						Y axis of the local reference frame of the grid.
		 *
		 * @param localZInWorld Vector in the global reference frame to set as
		 *						Z axis of the local reference frame of the grid.
		 *
		 * @note	The three specified vectors must form a basis,
		 *			although it is not required for the basis to be orthogonal.
		 */		virtual bool setLocalToWorldTransform(const Vector3d& worldOrigin,
											  const Vector3d& localXInWorld,
											  const Vector3d& localYInWorld,
											  const Vector3d& localZInWorld)
		{
			Matrix4d newLocalToWorld = SNS::Utilities::makeMatrix4d();
			// Set columns of transformation matrix (Boost uses row-major order)
			newLocalToWorld(0, 0) = localXInWorld(0); newLocalToWorld(0, 1) = localYInWorld(0); newLocalToWorld(0, 2) = localZInWorld(0); newLocalToWorld(0, 3) = worldOrigin(0);
			newLocalToWorld(1, 0) = localXInWorld(1); newLocalToWorld(1, 1) = localYInWorld(1); newLocalToWorld(1, 2) = localZInWorld(1); newLocalToWorld(1, 3) = worldOrigin(1);
			newLocalToWorld(2, 0) = localXInWorld(2); newLocalToWorld(2, 1) = localYInWorld(2); newLocalToWorld(2, 2) = localZInWorld(2); newLocalToWorld(2, 3) = worldOrigin(2);
			newLocalToWorld(3, 0) = 0;                newLocalToWorld(3, 1) = 0;                newLocalToWorld(3, 2) = 0;                newLocalToWorld(3, 3) = 1;

			return this->setLocalToWorldTransform(newLocalToWorld);
		}

		/**
		 * @brief	Sets the change of reference frame transformation from local
		 *			coordinates of the grid to world (global) coordinates.
		 *			Returns true if the procedure succeeded, false if it fails
		 *			since the resulting matrix is not invertible.
		 *
		 * @param worldOrigin_0, worldOrigin_1, worldOrigin_2
		 *			Coordinates expressed in the global reference frame
		 *			of a point to be set as origin of the local reference frame
		 *			of the grid (i.e. the point having [0,0,0] as local coords).
		 *
		 * @param localXInWorld_0, localXInWorld_1, localXInWorld_2
		 *			Coordinates expressed in the global reference frame of the
		 *			vector to be set as X axis of the local reference frame
		 *			of the grid.
		 *
		 * @param localYInWorld_0, localYInWorld_1, localYInWorld_2
		 *			Coordinates expressed in the global reference frame of the
		 *			vector to be set as Y axis of the local reference frame
		 *			of the grid.
		 *
		 * @param localZInWorld_0, localZInWorld_1, localZInWorld_2
		 *			Coordinates expressed in the global reference frame of the
		 *			vector to be set as Z axis of the local reference frame
		 *			of the grid.
		 *
		 * @note	The three specified vectors must form a basis,
		 *			although it is not required for the basis to be orthogonal.
		 */
		virtual bool setLocalToWorldTransform(double worldOrigin_0,
											  double worldOrigin_1,
											  double worldOrigin_2,
											  double localXInWorld_0,
											  double localXInWorld_1,
											  double localXInWorld_2,
											  double localYInWorld_0,
											  double localYInWorld_1,
											  double localYInWorld_2,
											  double localZInWorld_0,
											  double localZInWorld_1,
											  double localZInWorld_2)
		{			return this->setLocalToWorldTransform(
				SNS::Utilities::makeVector3d(worldOrigin_0, worldOrigin_1, worldOrigin_2),
				SNS::Utilities::makeVector3d(localXInWorld_0, localXInWorld_1, localXInWorld_2),
				SNS::Utilities::makeVector3d(localYInWorld_0, localYInWorld_1, localYInWorld_2),
				SNS::Utilities::makeVector3d(localZInWorld_0, localZInWorld_1, localZInWorld_2) );
		}

		/**
		 * @brief	Given the local coordinates passed as parameters,
		 *			this method computes the corresponding point expressed
		 *			in the global (world) reference frame.
		 *
		 * @param localCoords	Coordinates of a point in the LOCAL
		 *						reference frame of the grid.
		 *
		 * @throw std::out_of_range	If the coordinates are out of
		 *							the grid's bounds.
		 */
		virtual Vector3d fromLocalToWorldFrame(
								const Vector3u& localCoords) const
		{
			if( (localCoords(0) > (_nPoints[0]-1u)) ||
				(localCoords(1) > (_nPoints[1]-1u)) ||
				(localCoords(2) > (_nPoints[2]-1u))   )
			{
				Utilities::throwAndPrintError<std::out_of_range>(
					"DataGrid3D::fromLocalToWorldFrame() : "
					"Coordinates are out of bounds!"
				);
			}			Vector3d d_localCoords = SNS::Utilities::makeVector3d(localCoords(0),
																 localCoords(1),
																 localCoords(2));
			// Convert to homogeneous coordinates (4D vector with w=1)
			Vector4d homogeneous = SNS::Utilities::makeVector4d(d_localCoords(0), d_localCoords(1), d_localCoords(2), 1.0);
			Vector4d result4d = boost::numeric::ublas::prod(_localToWorld, homogeneous);
			// Extract 3D result
			return SNS::Utilities::makeVector3d(result4d(0), result4d(1), result4d(2));
		}

		/**
		 * @brief	Given the world coordinates passed as parameters,
		 *			this method computes the corresponding point or vector
		 *			expressed in the local reference frame of the grid.
		 *
		 *
		 * @param worldCoordinates	Coordinates of a point or vector in the
		 *							WORLD (global) reference frame. The last
		 *							coordinate (w) must be 1 in the case of a
		 *							point and 0 in the case of a vector.
		 *
		 * @note	The returned coordinates may refer to a point
		 *			lying outside the grid!
		 */		virtual Vector3d fromWorldToLocalFrame(
						const Vector4d& worldCoordinates) const
		{
			return SNS::Utilities::xyz(SNS::Utilities::multiply(_worldToLocal, worldCoordinates));
		}


        /**
         * @brief	Returns the value stored at the specified
         *			LOCAL coordinates of the grid.
         *
		 * @param localCoords The LOCAL coordinate of the element to retrieve.
		 *
		 * @throw std::out_of_range	If the coordinates are out of
		 *							the grid's bounds.
         */
		virtual T getValue(const Vector3u& localCoords) const
        {
			if( (localCoords(0) > (_nPoints[0]-1u)) ||
				(localCoords(1) > (_nPoints[1]-1u)) ||
				(localCoords(2) > (_nPoints[2]-1u)) )
            {
				Utilities::throwAndPrintError<std::out_of_range>(
					"DataGrid3D::getValue() : Coordinates are out of bounds!"
				);
            }

			return _values[localCoords(0)][localCoords(1)][localCoords(2)];
        }

        /**
         * @brief	Stores the value passed as parameter in the grid
         *			at the specified LOCAL coordinates.
         *
		 * @param localCoords	The LOCAL coordinate of the voxel to change.
         *
         * @param value		The new value to be stored in the grid.
         *
         * @throw std::out_of_range If local coordinates are out of
         *							the bounds of the grid.
		 *
		 * @note	Use this method only to change the value of a few elements
		 *			of the grid. If you need to set a new value for large part
		 *			of the volume, consider the use of the fill() method.
         */
		virtual void setValue(const Vector3u& localCoords,
							  const T& value)
        {
			if( (localCoords(0) > (_nPoints[0]-1u)) ||
				(localCoords(1) > (_nPoints[1]-1u)) ||
				(localCoords(2) > (_nPoints[2]-1u)) )
            {
                Utilities::throwAndPrintError<std::out_of_range>(
                    "DataGrid3D::setValue() : Coordinates are out of bounds!");
            }

			_values[localCoords(0)][localCoords(1)][localCoords(2)] = value;
        }

		/**
		 * @brief	Fills the grid with the data contained in the
		 *			three-dimensional vector passed as parameter.
		 *			The content of the vector is COPIED into the grid,
		 *			replacing the previous values.
		 *
		 * @throw std::invalid_argument	If the three-dimensional vector passed
		 *								as parameter and this grid contain a
		 *								different number of elements.
		 */
		virtual void fill(const std::vector<std::vector<std::vector<T>>>& data)
		{
			// Checks the size of the outer vector.
			bool sameSize = (data.size() == _nPoints[0]);
			// Checks the size of the inner vectors.
			for(uint i = 0; (i < _nPoints[0]) && sameSize; i++)
			{
				sameSize = sameSize && (data[i].size() == _nPoints[1]);
				for(uint j = 0; (j < _nPoints[1]) && sameSize; j++)
					sameSize = sameSize && (data[i][j].size() == _nPoints[2]);
			}

			if(!sameSize)
			{
				Utilities::throwAndPrintError<std::invalid_argument>(
					"DataGrid3D::fill() : The 3D vector passed as parameter "
					"has different sizes with respect to this grid!");
			}

			_values = data;
		}

		/**
		 * @brief	Fills the grid with the data contained in the
		 *			three-dimensional vector passed as parameter.
		 *			The content of the vector is MOVED into the grid,
		 *			replacing the previous values.
		 *
		 * @note	Since the content of the vector is MOVED (instead of
		 *			copied), after calling this method its state may not be
		 *			valid anymore.
		 *
		 * @throw std::invalid_argument	If the three-dimensional vector passed
		 *								as parameter and this grid contain a
		 *								different number of elements.
		 */
		virtual void fill(std::vector<std::vector<std::vector<T>>>&& data)
		{
			// Checks the size of the outer vector.
			bool sameSize = (data.size() == _nPoints[0]);
			// Checks the size of the inner vectors.
			for(uint i = 0; (i < _nPoints[0]) && sameSize; i++)
			{
				sameSize = sameSize && (data[i].size() == _nPoints[1]);
				for(uint j = 0; (j < _nPoints[1]) && sameSize; j++)
					sameSize = sameSize && (data[i][j].size() == _nPoints[2]);
			}

			if(!sameSize)
			{
				Utilities::throwAndPrintError<std::invalid_argument>(
					"DataGrid3D::fill() : The 3D vector passed as parameter "
					"has different sizes with respect to this grid!");
			}

			// Steals the data
			// BUG: For some unkown reason, the cost of move assignement
			// of std::vector is >= then copy assignement...
			_values = std::move(data);

			// Clears the source
			data.clear();
		}


		/**
		 * @brief	Returns true if this grid and the one passed as
		 *			parameter contains the same number of voxels along
		 *			the three dimensions, false otherwise.
		 */
		virtual bool sameSize(const DataGrid3D<T>& other) const
		{
			if(this->_nPoints[0] != other._nPoints[0]) return false;
			if(this->_nPoints[1] != other._nPoints[1]) return false;
			if(this->_nPoints[2] != other._nPoints[2]) return false;

			return true;
		}

		/**
		 * @brief	Returns true if this grid and the one passed as parameter
		 *			are "conformable", i.e. they have the same number of voxels
		 *			and the same local->world change of reference frame matrix.
		 */
		virtual bool isConformableTo(const DataGrid3D<T>& other) const
		{
			// Same number of voxels?
			if(!this->sameSize(other)) return false;

			// Same tranformation matrix
			for(int i = 0; i < 4; i++)
				for(int j = 0; j < 4; j++)
					if( !Utilities::fuzzyEqual(this->_localToWorld(i,j),
											   other._localToWorld(i,j)) )
						return false;

			return true;
		}

		/**
		 * @brief	Compares this grid with the one passed as parameter,
		 *			and returns true if the two grids have the same size
		 *			and contain the same values, false otherwise.
		 *
		 * @note	This method can be used only if the type parameter T
		 *			implements the equality operator (==).
		 */
		virtual bool sameContent(const DataGrid3D<T>& other) const
		{
			// Same number of voxels?
			if(!this->sameSize(other)) return false;			// Same values?
			for(unsigned x = 0; x < this->_nPoints[0]; x++)
				for(unsigned y = 0; y < this->_nPoints[1]; y++)
					for(unsigned z = 0; z < this->_nPoints[2]; z++)
					{
						if(compareValues(this->_values[x][y][z], other._values[x][y][z]))
							return false;
					}

			return true;
		}

		/**
		 * @brief	Compares this grid with the one passed as parameter,
		 *			and returns true if the two grids have the same size
		 *			and contain the same values, false otherwise.
		 *			The values are compared for equality with the function
		 *			passed as parameter.
		 */
		virtual bool sameContent(
						const DataGrid3D<T>& other,
						bool (*areEquals)(const T& v1, const T& v2) ) const
		{
			// Same number of voxels?
			if(!this->sameSize(other)) return false;

			// Same values?
			for(unsigned x = 0; x < this->_nPoints[0]; x++)
				for(unsigned y = 0; y < this->_nPoints[1]; y++)
					for(unsigned z = 0; z < this->_nPoints[2]; z++)
					{
						if(!areEquals(this->_values[x][y][z],
									  other._values[x][y][z])) return false;
					}

			return true;
		}

		/**
		 * @brief	Equality operator.
		 *			Two instances of DataGrid3D are considered equals if
		 *			they have the same name, the same content and the same
		 *			change of reference frame matrices.
		 *
		 * @note	This method can be used only if the type parameter T
		 *			implements the equality operator (==).
		 */
		virtual bool operator==(const DataGrid3D<T>& other) const
		{
			// If "other" is an instance of a different class return false.
			if( typeid(*this) != typeid(other) ) return false;

			return	(this->_name == other._name) &&
					(this->isConformableTo(other)) &&
					(this->sameContent(other));
		}

		/**
		 * @brief	Inequality operator.
		 *			Two instances of DataGrid3D are considered equals if
		 *			they have the same name, the same content and the same
		 *			change of reference frame matrices.
		 *
		 * @note	This method can be used only if the type parameter T
		 *			implements the equality operator (==).
		 */
		virtual bool operator!=(const DataGrid3D<T>& other) const
		{
			return	!((*this) == other);
		}


		/**
		 * @brief	Samples the value of a point within the grid.
		 *			The returned value is computed by linearly interpolating
		 *			the eight values of	the voxels surrounding the point.
		 *
		 * @param localCoords	LOCAL coordinates of the point for which an
		 *						associated value is requested.
		 *
		 * @note	This method can be used only if the type parameter T
		 *			implements the binary addition operator (T+T) and
		 *			and the multiplication with a scalar value (double*T).
		 *
		 * @throw std::out_of_range	If the coordinates are out
		 *							of the grid's bounds.
		 */
		T sampleValue(const Vector3d& localCoords) const
		{
			if( (localCoords(0) < 0.0) ||
				(localCoords(1) < 0.0) ||
				(localCoords(2) < 0.0) ||
				(localCoords(0) > double(_nPoints[0]-1u)) ||
				(localCoords(1) > double(_nPoints[1]-1u)) ||
				(localCoords(2) > double(_nPoints[2]-1u)) )
			{
				Utilities::throwAndPrintError<std::out_of_range>(
					"DataGrid3D::sampleValue() : Coordinates are out of bounds!"
				);
			}

			Cube cube;
			Vector3d normCoords;
			this->getCubeInfo(localCoords, cube, normCoords);

			return Utilities::trilinearInterpolation(cube.voxelsValues, normCoords);
		}


		/**
		 * @brief	Adds at each voxel of this grid the value of the
		 *			same voxel in the grid passed as parameter.
		 *			Finally returns this grid.
		 *
		 * @note	This method can be used only if the type parameter T
		 *			implements the addition assignment operator (T += T).
		 *
		 * @throw	std::invalid_argument If the two grids have different sizes.
		 */
		virtual DataGrid3D<T>& operator+=(const DataGrid3D<T>& other)
		{
			if(!this->sameSize(other))
			{
				Utilities::throwAndPrintError<std::invalid_argument>(
							"DataGrid3D::operator+= : Trying to "
							"sum two grids with different sizes!");
			}

			for(unsigned i = 0; i < _nPoints[0]; i++)
				for(unsigned j = 0; j < _nPoints[1]; j++)
					for(unsigned k = 0; k < _nPoints[2]; k++)
						this->_values[i][j][k] += other._values[i][j][k];

			return *this;
		}

		/**
		 * @brief	Subtract at each voxel of this grid the value of the
		 *			same voxel in the grid passed as parameter.
		 *			Finally returns this grid.
		 *
		 * @note	This method can be used only if the type parameter T
		 *			implements the subtraction assignment operator (T -= T).
		 *
		 * @throw	std::invalid_argument If the two grids have different sizes.
		 */
		virtual DataGrid3D<T>& operator-=(const DataGrid3D<T>& other)
		{
			if(!this->sameSize(other))
			{
				Utilities::throwAndPrintError<std::invalid_argument>(
							"DataGrid3D::operator-= : Trying to "
							"sum two grids with different sizes!");
			}

			for(unsigned i = 0; i < _nPoints[0]; i++)
				for(unsigned j = 0; j < _nPoints[1]; j++)
					for(unsigned k = 0; k < _nPoints[2]; k++)
						this->_values[i][j][k] -= other._values[i][j][k];

			return *this;
		}

		/**
		 * @brief	Multiplies at each voxel of this grid the value of the
		 *			same voxel in the grid passed as parameter.
		 *			Finally returns this grid.
		 *
		 * @note	This method can be used only if the type parameter T
		 *			implements the multiplication assignement operator (T *= T).
		 *
		 * @throw	std::invalid_argument If the two grids have different sizes.
		 */
		virtual DataGrid3D<T>& operator*=(const DataGrid3D<T>& other)
		{
			if(!this->sameSize(other))
			{
					Utilities::throwAndPrintError<std::invalid_argument>(
								"DataGrid3D::operator*= : Trying to "
								"sum two grids with different sizes!");
			}			for(unsigned i = 0; i < _nPoints[0]; i++)
				for(unsigned j = 0; j < _nPoints[1]; j++)
					for(unsigned k = 0; k < _nPoints[2]; k++)
						multiplyAssignValues(this->_values[i][j][k], other._values[i][j][k]);

			return *this;
		}

		/**
		 * @brief	Multiplies each voxel of this grid with the scalar
		 *			passed as parameter. Finally returns this grid.
		 *
		 * @note	This method can be used only if the type parameter T
		 *			implements the multiplication assignement operator between
		 *			T instances and scalar values (T *= double).
		 */
		virtual DataGrid3D<T>& operator*=(double factor)
		{
			for(unsigned i = 0; i < _nPoints[0]; i++)
				for(unsigned j = 0; j < _nPoints[1]; j++)
					for(unsigned k = 0; k < _nPoints[2]; k++)
						this->_values[i][j][k] *= factor;

			return *this;
		}

		// Template specialization helpers for vector operations
		template<typename U>
		static bool compareValues(const U& v1, const U& v2) {
			return v1 != v2;
		}
		
		template<typename U>
		static void multiplyAssignValues(U& v1, const U& v2) {
			v1 *= v2;
		}
		
		// Specializations for boost vectors
		static bool compareValues(const boost::numeric::ublas::vector<double>& v1, 
								const boost::numeric::ublas::vector<double>& v2) {
			return !SNS::Utilities::vectorsEqual(v1, v2);
		}
		
		static void multiplyAssignValues(boost::numeric::ublas::vector<double>& v1, 
										const boost::numeric::ublas::vector<double>& v2) {
			SNS::Utilities::elementWiseMultiplyAssign(v1, v2);
		}
		
		static bool compareValues(const boost::numeric::ublas::vector<unsigned int>& v1, 
								const boost::numeric::ublas::vector<unsigned int>& v2) {
			return !SNS::Utilities::vectorsEqual(v1, v2);
		}
		
		static void multiplyAssignValues(boost::numeric::ublas::vector<unsigned int>& v1, 
										const boost::numeric::ublas::vector<unsigned int>& v2) {
			SNS::Utilities::elementWiseMultiplyAssign(v1, v2);
		}
	};

} }

#endif // DATAGRID3D_H
