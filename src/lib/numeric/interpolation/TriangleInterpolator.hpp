#ifndef LIBGCM_TRIANGLEINTERPOLATOR_HPP
#define LIBGCM_TRIANGLEINTERPOLATOR_HPP

#include <vector>
#include <lib/linal/linal.hpp>
#include <lib/util/Types.hpp>


namespace gcm {
template<class TValue>
class TriangleInterpolator {
// TODO - some non-dimensionalization?

public:

	/**
	 * Linear interpolation in plain triangle
	 * @param c_i and v_i - points and values
	 * @param q point to interpolate
	 */
	static TValue interpolate(const Real2& c1, const TValue v1,
	                          const Real2& c2, const TValue v2,
	                          const Real2& c3, const TValue v3,
	                          const Real2& q) {

		real det = linal::determinant(c1(0), c1(1), 1.0,
		                              c2(0), c2(1), 1.0,
		                              c3(0), c3(1), 1.0);

		TValue det1 = linal::determinant(v1, c1(1), 1,
		                                 v2, c2(1), 1,
		                                 v3, c3(1), 1);

		TValue det2 = linal::determinant(c1(0), v1, 1,
		                                 c2(0), v2, 1,
		                                 c3(0), v3, 1);

		TValue det3 = linal::determinant(c1(0), c1(1), v1,
		                                 c2(0), c2(1), v2,
		                                 c3(0), c3(1), v3);

		// v(x, y) = ax + by + c
		auto a = det1 / det;
		auto b = det2 / det;
		auto c = det3 / det;

		return a * q(0) + b* q(1) + c;
	}

};


}

#endif // LIBGCM_TRIANGLEINTERPOLATOR_HPP
