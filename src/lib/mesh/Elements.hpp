#ifndef LIBGCM_ELEMENTS_HPP
#define LIBGCM_ELEMENTS_HPP

#include <lib/linal/linal.hpp>

namespace gcm {
namespace elements {
/**
 * Convex geometric primitives by sets of their vertices
 * @tparam Point representation of point
 * @tparam TN number of points
 */
template<typename Point, int TN>
struct Element {
	static const int N = TN;
	Point p[N];
	bool inner = false;
};

template<typename Point> using Triangle = Element<Point, 3>;
}
}

#endif // LIBGCM_ELEMENTS_HPP
