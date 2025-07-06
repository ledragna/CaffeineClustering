#ifndef MATHFU_UTILITIES_H
#define MATHFU_UTILITIES_H

#include <mathfu/matrix.h>
#include <mathfu/vector.h>

#include <QVector3D>
#include <QVector4D>
#include <QMatrix4x4>

#include <Utilities/floating_point.h>

namespace SNS { namespace Utilities
{
	/**
	 * @brief	Compares two mathfu::Matrix for equality.
	 *			In particular, m1 and m2 are considered equals if the relative
	 *			differerence between their elements is smaller than a epsilon.
	 *			The default value for epsilon is 1e-4.
	 *
	 * @note    For details see:
	 *			http://www.boost.org/doc/libs/1_59_0/libs/test/doc/html/boost_test/testing_tools/extended_comparison/floating_point/floating_points_comparison_theory.html
	 */
	///@{
	template<class T, int rows, int columns>
	bool fuzzyEqual(const mathfu::Matrix<T, rows, columns>& m1,
					 const mathfu::Matrix<T, rows, columns>& m2,
					 T epsilon)
	{
		bool equal = true;
		for(int i = 0; (i < (rows*columns)) && equal; i++)
			equal = equal && fuzzyEqual(m1[i], m2[i], epsilon);

		return equal;
	}

	template<class T, int rows, int columns>
	bool fuzzyEqual( const mathfu::Matrix<T, rows, columns>& m1,
					 const mathfu::Matrix<T, rows, columns>& m2)
	{
		return fuzzyEqual(m1, m2, (T) 1e-4f);
	}
	///!@}


	/**
	 * @brief	Returns the specified column of the
	 *			mathfu::Matrix passed as parameter.
	 */
	template<class T, int rows, int columns>
	mathfu::Vector<T, rows> getColumn(const mathfu::Matrix<T, rows, columns>& m,
									  int columnNumber)
	{
		mathfu::Vector<T, rows> columnVect;
		for(int r = 0; r < rows; r++)
			columnVect[r] = m(r, columnNumber);
		return columnVect;
	}


	/**
	 * @brief	Compares two mathfu::Vector for equality.
	 *			In particular, v1 and v2 are considered equals if the relative
	 *			differerence between each pair of components is smaller than a
	 *			epsilon. The default value for epsilon is 1e-4.
	 *
	 * @note    For details see:
	 *			http://www.boost.org/doc/libs/1_59_0/libs/test/doc/html/boost_test/testing_tools/extended_comparison/floating_point/floating_points_comparison_theory.html
	 */
	///@{
	template<class T, int n>
	bool fuzzyEqual( const mathfu::Vector<T, n>& v1,
					 const mathfu::Vector<T, n>& v2,
					 T epsilon)
	{
		bool equal = true;
		for(int i = 0; (i < n) && equal; i++)
			equal = equal && fuzzyEqual(v1[i], v2[i], epsilon);

		return equal;
	}

	template<class T, int n>
	bool fuzzyEqual( const mathfu::Vector<T, n>& v1,
					 const mathfu::Vector<T, n>& v2)
	{
		return fuzzyEqual(v1, v2, (T) 1e-4f);
	}
	///!@}


	/**
	 *  Conversion functions from Mathfu vectors/matrices
	 *  to Qt's vectors/matrices.
	 */
	///@{
	QVector3D to_QVector3D(const mathfu::Vector<double,3>& p);
	QVector4D to_QVector4D(const mathfu::Vector<double,4>& p);
	QMatrix4x4 to_QMatrix4x4(const mathfu::Matrix<double,4>& mat);
	///!@}


	/**
	 *  Conversion functions from Qt's vectors/matrices
	 *  to Mathfu vectors/matrices.
	 */
	///@{
	mathfu::Vector<double,3> to_Mathfu_dvec3(const QVector3D& p);
	mathfu::Vector<double,4> to_Mathfu_dvec4(const QVector4D& p);
	///!@}
} }

#endif // MATHFU_UTILITIES_H
