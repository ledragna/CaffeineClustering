#ifndef FLOATING_POINT_H
#define FLOATING_POINT_H

/**
 * @author Andrea Salvadori
 */
namespace CFF { namespace Utilities
{

/**
 * @brief	Compares two floating point numbers for equality.
 *			In particular, x and y are considered equals if the
 *			relative differerence is smaller than epsilon.
 *			The default value for epsilon is 1e-4.
 *
 * @note    For details see:
 *			http://www.boost.org/doc/libs/1_59_0/libs/test/doc/html/boost_test/testing_tools/extended_comparison/floating_point/floating_points_comparison_theory.html
 */
///@{
bool fuzzyEqual(double x, double y, double epsilon);
bool fuzzyEqual(const double& x, const double& y); // const T& has been used to work with template stuff
bool fuzzyEqual(float x, float y, float epsilon);
bool fuzzyEqual(const float& x, const float& y); // const T& has been used to work with template stuff
///!@}


/**
 * @brief	These functions implement the binary comparing
 *			operators between two floating point numbers.
 *			The comparison is performed in an approximate
 *			("fuzzy") manner, by exploiting the fuzzyEqual()
 *			function.
 */
///@{

/**
 * @brief	Approximate < operator. Compares two floating point
 *			numbers according to the following expression:
 *			(x < y) && (!fuzzyEqual(x, y, epsilon))
 *			The default value for epsilon is 1e-4.
 */
///@{
bool fuzzyLessThan(double x, double y, double epsilon);
bool fuzzyLessThan(float x, float y, float epsilon);
bool fuzzyLessThan(const double& x, const double& y); // const T& has been used to work with template stuff
bool fuzzyLessThan(const float& x, const float& y); // const T& has been used to work with template stuff
///!@}

/**
 * @brief	Approximate <= operator. Compares two floating point
 *			numbers according to the following expression:
 *			(x < y) || fuzzyEqual(x, y, epsilon)
 *			The default value for epsilon is 1e-4.
 */
///@{
bool fuzzyEqualOrLessThan(double x, double y, double epsilon);
bool fuzzyEqualOrLessThan(float x, float y, float epsilon);
bool fuzzyEqualOrLessThan(const double& x, const double& y); // const T& has been used to work with template stuff
bool fuzzyEqualOrLessThan(const float& x, const float& y); // const T& has been used to work with template stuff
///!@}

/**
 * @brief	Approximate > operator. Compares two floating point
 *			numbers according to the following expression:
 *			(x > y) && (!fuzzyEqual(x, y, epsilon))
 *			The default value for epsilon is 1e-4.
 */
///@{
bool fuzzyGreaterThan(double x, double y, double epsilon);
bool fuzzyGreaterThan(float x, float y, float epsilon);
bool fuzzyGreaterThan(const double& x, const double& y); // const T& has been used to work with template stuff
bool fuzzyGreaterThan(const float& x, const float& y); // const T& has been used to work with template stuff
///!@}

/**
 * @brief	Approximate <= operator. Compares two floating point
 *			numbers according to the following expression:
 *			(x > y) || fuzzyEqual(x, y, epsilon)
 *			The default value for epsilon is 1e-4.
 */
///@{
bool fuzzyEqualOrGreaterThan(double x, double y, double epsilon);
bool fuzzyEqualOrGreaterThan(float x, float y, float epsilon);
bool fuzzyEqualOrGreaterThan(const double& x, const double& y); // const T& has been used to work with template stuff
bool fuzzyEqualOrGreaterThan(const float& x, const float& y); // const T& has been used to work with template stuff
///!@}

///!@}

} }

#endif // FLOATING_POINT_H
