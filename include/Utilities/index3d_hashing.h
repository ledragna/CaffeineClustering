#ifndef INDEX3D_HASHING_H
#define INDEX3D_HASHING_H

#include <array>
#include <cstddef>

/**
 * @author Andrea Salvadori and Marco Fusè
 */
namespace CFF { namespace Utilities
{

	/**
	 * @brief	Function object implementing an hash function for a 3D index.
	 *			Can be used to create unsorted sets or unsorted maps of 3D indices.
	 *
	 *			The hash values are computed using the "prime number multiplication" algorithm,
	 *			like the one implemented by the HashCodeBuilder class of the Apache framework.
	 *
	 * @see	https://commons.apache.org/proper/commons-lang/apidocs/org/apache/commons/lang3/builder/HashCodeBuilder.html
	 * @see http://myeyesareblind.com/2017/02/06/Combine-hash-values/
	 */
	struct Index3DHasher
	{
		size_t operator()(std::array<unsigned,3> p) const;
	};

	/**
	 * @brief	Function object implementing a function which tests two 3D indices for equality.
	 *			Can be used to create unsorted sets or unsorted maps of 3D indices.
	 */
	struct Index3DComparator
	{
		bool operator()(const std::array<unsigned,3>& p1,
						const std::array<unsigned,3>& p2) const;
	};

} }

#endif // INDEX3D_HASHING_H
