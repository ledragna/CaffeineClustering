#ifndef ARRAY3D_H
#define ARRAY3D_H

#include <stdexcept>
#include <QDebug>

#include <Utilities/log.h>

/**
 * @author Andrea Salvadori
 */
namespace CFF { namespace Utilities
{

/**
 * @brief	Convenient class implementing a 3D array whose dimensions
 *			can be specified at runtime.
 *
 *			This class takes care of allocating the required memory
 *			on the heap and deteting the memory when the objects gets
 *			distroyed. Move constructor and assignement operator are
 *			provided to efficiently return an Array3D from a function
 *			without the need to heap-allocate Array3D instances and
 *			avoiding copies.
 */
template<typename T> class Array3D
{
private:
	T* _data;
	size_t _dimx;
	size_t _dimy;
	size_t _dimz;

public:

	/**
	 * @brief	Creates a 3D array with the specified number
	 *			of elements along each dimension.
	 *
	 * @param dimx Number of elements along the first dimension.
	 * @param dimy Number of elements along the second dimension.
	 * @param dimz Number of elements along the third dimension.
	 */
	Array3D(size_t dimx, size_t dimy, size_t dimz)
		: _dimx(dimx), _dimy(dimy), _dimz (dimz)
	{
		if( (dimx == 0) || (dimy == 0) || (dimz == 0) )
			Utilities::throwAndPrintError<std::invalid_argument>("Array3D: invalid size!");

		size_t dataSize = _dimx*_dimy*_dimz;
		_data = new T[dataSize];
		for(size_t i = 0; i < dataSize; ++i) _data[i] = T();
	}

	/** @brief Destructor */
	~Array3D() { delete[] _data; }

	/** @brief Copy constructor */
	Array3D(const Array3D& other)
		: _dimx(other._dimx), _dimy(other._dimy), _dimz (other._dimz)
	{
		//qDebug() << "Array3D copy constructor";

		size_t dataSize = _dimx*_dimy*_dimz;
		if(dataSize > 0)
		{
			_data = new T[dataSize];
			std::copy_n(other._data, dataSize, _data);
		}
		else _data = NULL;
	}

	/** @brief Move constructor */
	Array3D(Array3D&& other)
	{
		//qDebug() << "Array3D move constructor";

		// Shallow copy
		_dimx = other._dimx;
		_dimy = other._dimy;
		_dimz = other._dimz;
		_data = other._data;

		// Revert the parameter to a minimal consistent state by clearing it
		other._dimx = 0;
		other._dimy = 0;
		other._dimz = 0;
		other._data = NULL;
	}

	/** @brief Copy assignement operator */
	Array3D& operator=(const Array3D& other)
	{
		//qDebug() << "Array3D copy assignement";

		// Deletes the previously allocated memory
		delete[] _data;

		// Deep copy
		_dimx = other._dimx;
		_dimy = other._dimy;
		_dimz = other._dimz;
		size_t dataSize = _dimx*_dimy*_dimz;
		if(dataSize > 0)
		{
			_data = new T[dataSize];
			std::copy_n(other._data, dataSize, _data);
		}
		else _data = NULL;

		// Returns a reference to this object
		return *this;
	}

	/** @brief Move assignement operator */
	Array3D& operator=(Array3D&& other)
	{
		//qDebug() << "Array3D move assignement";

		// Deletes the previously allocated memory
		delete[] _data;

		// Shallow copy
		_dimx = other._dimx;
		_dimy = other._dimy;
		_dimz = other._dimz;
		_data = other._data;

		// Revert the parameter to a minimal consistent state by clearing it
		other._dimx = 0;
		other._dimy = 0;
		other._dimz = 0;
		other._data = NULL;

		// Returns a reference to this object
		return *this;
	}

	/** @brief Subscript operator. */
	T& operator() (size_t i, size_t j, size_t k)
	{
		if((i >= _dimx) || (j >= _dimy) || (k >= _dimz))
			Utilities::throwAndPrintError<std::out_of_range>("Array3D: out of bounds!");

		size_t index = k + (_dimz * j) + (_dimz * _dimy * i);
		return _data[index];
	}

	/** @brief Subscript operator. Const version. */
	const T& operator() (size_t i, size_t j, size_t k) const
	{
		if((i >= _dimx) || (j >= _dimy) || (k >= _dimz))
			Utilities::throwAndPrintError<std::out_of_range>("Array3D: out of bounds!");

		size_t index = k + (_dimz * j) + (_dimz * _dimy * i);
		return _data[index];
	}

	/** @brief Returns the number of elements of the array along the first dimension. */
	size_t sizeX() const { return _dimx; }

	/** @brief Returns the number of elements of the array along the first dimension. */
	size_t sizeY() const { return _dimy; }

	/** @brief Returns the number of elements of the array along the first dimension. */
	size_t sizeZ() const { return _dimz; }

	/** @brief Returns the total number of elements contained in the array. */
	size_t totalSize() const { return _dimx*_dimy*_dimz; }

	/**
	 * @brief	Clears the memory managed by the array.
	 *			Trying to access to a cleared array results in an exception.
	 */
	void clear()
	{
		delete[] _data;
		_dimx = 0;
		_dimy = 0;
		_dimz = 0;
		_data = NULL;
	}

	/**
	 * @brief	Returs true is the array has been cleared
	 *			(and it is not usable anymore), false otherwise.
	 * @see Array3D::clear() method.
	 */
	bool isCleared() const
	{ return (_dimx == 0) || (_dimy == 0) || (_dimz == 0); }
};

} }

#endif // ARRAY3D_H
