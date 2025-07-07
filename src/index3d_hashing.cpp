#include "Utilities/index3d_hashing.h"
#include <functional>

using namespace CFF::Utilities;

size_t Index3DHasher::operator()(std::array<unsigned,3> p) const
{
	const size_t k1 = 13; // must be a prime number
	const size_t k2 = 71; // must be a prime number

	size_t hx = std::hash<unsigned>()(p[0]);
	size_t hy = std::hash<unsigned>()(p[1]);
	size_t hz = std::hash<unsigned>()(p[2]);

	size_t hash = k1;
	hash = (hash * k2) + hx;
	hash = (hash * k2) + hy;
	hash = (hash * k2) + hz;

	return hash;
}

bool Index3DComparator::operator()(const std::array<unsigned,3>& p1,
								   const std::array<unsigned,3>& p2) const
{
	return (p1[0] == p2[0]) && (p1[1] == p2[1]) && (p1[2] == p2[2]);
}
