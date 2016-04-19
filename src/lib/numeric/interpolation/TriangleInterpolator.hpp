#ifndef LIBGCM_TRIANGLEINTERPOLATOR_HPP
#define LIBGCM_TRIANGLEINTERPOLATOR_HPP

#include <lib/linal/linal.hpp>


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

		// v(x, y) = ax + by + c

		linal::Matrix33 A = {c1(0), c1(1), 1.0,
		                     c2(0), c2(1), 1.0,
		                     c3(0), c3(1), 1.0};
		
		linal::VECTOR<3, TValue> b = {v1, v2, v3};
		
		auto abc = linal::solveLinearSystem(A, b);
		
		return abc(0) * q(0) + abc(1) * q(1) + abc(2);
	}

};


}

#endif // LIBGCM_TRIANGLEINTERPOLATOR_HPP
