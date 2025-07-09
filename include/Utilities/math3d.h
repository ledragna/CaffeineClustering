#ifndef MATH3D_H
#define MATH3D_H

#include <vector>
#include <array>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>

#include <Utilities/floating_point.h>
#include <Utilities/log.h>

// Type aliases for easier migration
using Vector3d = boost::numeric::ublas::vector<double>;

/**
 * @author Andrea Salvadori
 */
namespace CFF { namespace Utilities
{



/**
 * @brief	This structure is used to return the results of the
 *			computeGlobalAxesSimilarity() function, which computes
 *			an association between these vectors and the 3 axes of the
 *			global reference frame, by searching for each axis the vector
 *			who have the direction more "similar" to the considered axis.
 *
 *			All the fields of this structure are shorts.
 *			A generic field named is_VN_Almost_A_ will contain:
 *			1	if the vector _VN_ is oriented in a way "similar" to the
 *				global _A_ axis.
 *		   -1	if the vector _VN_ is oriented in a way "similar" to the
 *				opposite of the global _A_ axis.
 *			0	if the vector _VN_ is not consideref to be oriented in a
 *				way "similar" to the global _A_ axis.
 */
struct AxesSimilarity
{
	short isV1SimilarToX;
	short isV1SimilarToY;
	short isV1SimilarToZ;
	short isV2SimilarToX;
	short isV2SimilarToY;
	short isV2SimilarToZ;
	short isV3SimilarToX;
	short isV3SimilarToY;
	short isV3SimilarToZ;
};



/** * @brief	Given 3 linearly indipendent vectors, this function computes
 *			an association between these vectors and the 3 axes of the
 *			global reference frame, by searching for each axis the vector
 *			who have the direction more "similar" to the considered axis.
 *
 * @return An AxesSimilarity structure encoding the computed association.
 */

/**
 * @brief	Given an input vector "inoutPseudoZ" in the 3D space, this function
 *			normalize "inoutPseudoZ" and then computes other two unit vectors
 *			such to form an orthonormal basis. The two computed vectors are
 *			finally stored in "outPseudoY" and "outPseudoX".
 *
 * @throw std::invalid_argument If inoutVector1 is a null vector.
 */

/**
 * @brief	Given an input vector "inoutPseudoZ" in the 3D space, this function
 *			normalize "inoutPseudoZ" and then computes other two unit vectors
 *			such to form an orthonormal basis. The two computed vectors are
 *			finally stored in "outPseudoY" and "outPseudoX".
 *
 * @throw std::invalid_argument If inoutVector1 is a null vector.
 */
void findOrthonormalBasis(Vector3d& inoutPseudoZ,
						  Vector3d& outPseudoY,
						  Vector3d& outPseudoX);

/**
 * @brief	Converts the provided window coordinates in the 3D space.
 *
 *			The MP parameter is a transform matrix containing at least the
 *			projection matrix. The nature of the matrix passed as parameter
 *			will determine the reference frame of the returned point.
 *			In particular, if MP is:
 *			- A projection matrix, then the coordinates of the returned point
 *			  will be in camera space.
 *			- A view-projection matrix, then the coordinates of the returned
 *			  point will be in world space.
 *			- A model-view-projection matrix, then the coordinates of the
 *			  returned point will be in world space.
 *
 * @param windowCoords	Coordinates in window space of the point to convert.
 *						The z component must be in the range [0,1].
 *
 * @param viewPort	Origin and size of the viewport, in window coordinates.
 *
 * @param MP	Can be a projecton matrix, a view-projection matrix or a
 *				model-view-projection matrix.
 * * @return A point in the 3D space corresponding to the provided window coordinates.
 */


/**
 * @brief Computes and returns a linear combination of elements.
 *
 * @param elements Vector containing the elements to be linearly combinated.
 * @param factors  Vector containing the scalar factors of the linear combination.
 *
 * @throw std::invalid_argument	If the two vector passed as parameters
 *								have different size.
 *
 * @note	This method can be used only if the type parameter T
 *			implements the addition assignment operator (T += T) and
 *			the multiplication assignment operator between T instances
 *			and scalar values (T *= double).
 */
template<typename T>
T linearCombination(const std::vector<const T*>& elements,
					const std::vector<double>& factors)
{
	/* Parameters cheking */

	if((elements.size() == 0) || (factors.size() == 0))
	{
		Utilities::throwAndPrintError<std::invalid_argument>(
					"Utilities::linearCombination : Empty vector!");
	}

	if(elements.size() != factors.size())
	{
		Utilities::throwAndPrintError<std::invalid_argument>(
				"Utilities::linearCombination : 'elements' and"
				" 'factors' vectors have different size!");
	}

	/* Computes the linear combination */
	T res( *(elements[0]) );
	if( !Utilities::fuzzyEqual(1.0,factors[0],1e-9) ) res *= factors[0];

    for(size_t i = 1; i < elements.size(); i++)
	{
		if( !Utilities::fuzzyEqual(0.0,factors[i],1e-9) )
		{
			if( !Utilities::fuzzyEqual(1.0,factors[i],1e-9) )
			{
				T tmp( *(elements[i]) );
				tmp *= factors[i];
				res += tmp;
			}
			else
				res += *(elements[i]);
		}
	}

	return res;
}

/**
 *	@brief	Given 8 values associated to the vertices of a cube and
 *			the normalized coordinates of a point within the cube,
 *			this function computes an approximated value for the
 *			specified point by interpolating linearly the values of
 *			the vertices of the cube.
 *
 * @param	cubeValues	Array containing the values associated to the 8
 *						vertices of the cube. The values must be provided
 *						in the following order:
 *						{V000, V001, V010, V011, V100, V101, V110, V111}
 *
 * @param	normCoords	Normalized coordinates of a point within the cube.
 *						Must be in the range [0,1]. If not, they will be
 *						claped to that range.
 *
 * @note	The type parameter T must implement the binary sum operator
 *			(T+T) and the multiplication with a scalar value (double*T).
 */
template<typename T>
T trilinearInterpolation(const std::array<T,8>& cubeValues,
						 const Vector3d& normCoords)
{	// Clamps the provided coordinates in the range [0,1]
	double normCoordsX = std::min(1.0, std::max(0.0, normCoords(0)));
	double normCoordsY = std::min(1.0, std::max(0.0, normCoords(1)));
	double normCoordsZ = std::min(1.0, std::max(0.0, normCoords(2)));

	// Little optimization: precomputes 1-normCoordsX/Y/Z
	double normCoordsNegX = 1.0-normCoordsX;
	double normCoordsNegY = 1.0-normCoordsY;
	double normCoordsNegZ = 1.0-normCoordsZ;

	/*
	 * Computes and returns the tri-linear interpolation of the cube's values
	 *
	 * Vxyz =	V000 (1 - x) (1 - y) (1 - z) +
	 *			V001 (1 - x) (1 - y)    z    +
	 *			V010 (1 - x)    y    (1 - z) +
	 *			V011 (1 - x)    y       z    +
	 *			V100    x    (1 - y) (1 - z) +
	 *			V101    x    (1 - y)    z    +
	 *			V110    x       y    (1 - z) +
	 *			V111    x       y       z
	 */
	return	( (normCoordsNegX * normCoordsNegY * normCoordsNegZ) * cubeValues[0] ) +
			( (normCoordsNegX * normCoordsNegY * normCoordsZ   ) * cubeValues[1] ) +
			( (normCoordsNegX * normCoordsY    * normCoordsNegZ) * cubeValues[2] ) +
			( (normCoordsNegX * normCoordsY    * normCoordsZ   ) * cubeValues[3] ) +
			( (normCoordsX    * normCoordsNegY * normCoordsNegZ) * cubeValues[4] ) +
			( (normCoordsX    * normCoordsNegY * normCoordsZ   ) * cubeValues[5] ) +
			( (normCoordsX    * normCoordsY    * normCoordsNegZ) * cubeValues[6] ) +
			( (normCoordsX    * normCoordsY    * normCoordsZ   ) * cubeValues[7] );
}


/**
 * @brief	Computes the intersection between a ray and a plane.
 *
 *			In particular, let's consider a plane passing through the point
 *			"planePoint" and having the versor "planeNormal" as normal. Let's
 *			also consider the ray with equation Q = rayOrigin + [t * rayDir],
 *			where "rayOrigin" is a point on the ray, "rayDir" is the
 *			direction of the ray, t is a scalar parameter and Q is the point
 *			on the ray corresponding to a specific value of t.
 *			This function returns true if the ray hits the plane, false otherwise.
 *			In the case of hit, the value of "t" corresponding to the
 *			intersection point will be stored in the "out_t" output parameter.
 *			The actual intersection point can the be computed with the equation
 *			of the ray presented above.
 *
 * @param rayOrigin		The coordinates of a point lying on the ray.
 *
 * @param rayDirection	The direction of the ray. May or may not have unit lenght.
 *						The resulting values of "t" will depend on the lenght
 *						of such vector.
 *
 * @param planePoint	The coordinates of a point lying on the plane.
 *
 * @param planeNormal	The normal versor defining the orientation of the plane.
 *						If its length is > 1, it will be normalized.
 *
 * @param out_t			Output parameter. In the case of success, will contain
 *						the value of "t" for the intersection point.
 *						In the case of miss, the value stored in this parameter
 *						is meaningless.
 *
 * @return	Returns true if the ray intersects the plane, false otherwise.
 */


/**
 * @brief	Computes the intersection between a ray and a plane.
 *
 *			In particular, let's consider an plane passing through the point
 *			"planePoint" and having the versor "planeNormal" as normal. Let's
 *			also consider the ray with equation Q = rayOrigin + [t * rayDir],
 *			where "rayOrigin" is a point on the ray, "rayDir" is the
 *			direction of the ray, t is a scalar parameter and Q is the point
 *			on the ray corresponding to a specific value of t.
 *			This function returns true if the ray hits the plane, false otherwise.
 *			In the case of hit, the value of "t" corresponding to the
 *			intersection point will be stored in the "out_t" output parameter.
 *			The actual intersection point can the be computed with the equation
 *			of the ray presented above.
 *
 * @param rayOrigin		The coordinates of a point lying on the ray.
 *
 * @param rayDirection	The direction of the ray. May or may not have unit lenght.
 *						The resulting values of "t" will depend on the lenght
 *						of such vector.
 *
 * @param planePoint	The coordinates of a point lying on the plane.
 *
 * @param planeNormal	The normal versor defining the orientation of the plane.
 *						If its length is > 1, it will be normalized.
 *
 * @param out_t			Output parameter. In the case of success, will contain
 *						the value of "t" for the intersection point.
 *						In the case of miss, the value stored in this parameter
 *						is meaningless.
 *
 * @return	Returns true if the ray intersects the plane, false otherwise.
 */
bool rayPlaneIntersection(	const Vector3d& rayOrigin,
							const Vector3d& rayDir,
							const Vector3d& planePoint,
							const Vector3d& planeNormal,
							double& out_t);


/**
 * @brief	Computes the intersection between a ray and an
 *			Axis-Aligned Bounding Box (AABB).
 *
 *			In particular, let's consider an AABB ranging from aabbMinXYZ to
 *			aabbMaxXYZ, and a ray with equation Q = rayOrigin + [t * rayDir],
 *			where "rayOrigin" is a point on the ray, "rayDir" is the
 *			direction of the ray, t is a scalar parameter and Q is the point
 *			on the ray corresponding to a specific value of t.
 *			This function returns true if the ray hits the AABB, false otherwise.
 *			In the case of hit, the values of "t" corresponding to the two
 *			intersection points will be stored in outTmin and outTmax.
 *			The actual intersection points can be computed with the equation
 *			of the ray presented above.
 *
 *			The implementation is based on the one presented in the book:
 *			"Real-time collision detection" by Ericson, Christer.
 *			CRC Press, 2004. Pag. 179-181
 *
 *
 * @param rayOrigin		The coordinates of a point lying on the ray.
 *
 * @param rayDirection	The direction of the ray. May or may not have unit lenght.
 *						The resulting values of "t" will depend on the lenght
 *						of such vector.
 *
 * @param aabbMinXYZ	The vertex of the AABB having the minimum coordinates.
 *
 * @param aabbMaxXYZ	The vertex of the AABB having the maximum coordinates.
 *
 * @param outTmin		Output parameter. In the case of success, will contain
 *						the value of "t" for the first intersection point.
 *						In the case of miss, the value stored in this parameter
 *						is meaningless.
 *
 * @param outTmax		Output parameter. In the case of success, will contain
 *						the value of "t" for the second intersection point.
 *						In the case of miss, the value stored in this parameter
 *						is meaningless.
 *
 * @return	Returns true if the ray intersects the AABB, false otherwise.
 *
 * @see		"Real-time collision detection" by Ericson, Christer.
 *			CRC Press, 2004. Pag. 179-181
 */
bool rayAABBIntersection(const Vector3d& rayOrigin,
						 const Vector3d& rayDir,
						 const Vector3d& aabbMinXYZ,
						 const Vector3d& aabbMaxXYZ,
						 double& outTmin, double& outTmax);

} }

#endif // MATH3D_H
