#ifndef LIBGCM_TRIANGLEINTERPOLATOR_HPP
#define LIBGCM_TRIANGLEINTERPOLATOR_HPP

#include <libgcm/linal/linal.hpp>


namespace gcm {
template<class TValue>
class TriangleInterpolator {
public:
	typedef linal::VECTOR<2, TValue>              Gradient;
	typedef linal::SYMMETRIC_MATRIX<2, TValue>    Hessian;
	
	/** Interpolation or extrapolation */
	static bool isInterpolation(const Real3 lambda) {
		return lambda(0) > -EQUALITY_TOLERANCE &&
		       lambda(1) > -EQUALITY_TOLERANCE &&
		       lambda(2) > -EQUALITY_TOLERANCE;
	}
	
	/**
	 * Linear interpolation in plain non-degenerate triangle
	 * @param c_i and v_i - points and values
	 * @param q point to interpolate
	 */
	static TValue interpolate(const Real2& c0, const TValue v0,
	                          const Real2& c1, const TValue v1,
	                          const Real2& c2, const TValue v2,
	                          const Real2& q) {
		Real3 lambda = linal::barycentricCoordinates(c0, c1, c2, q);
		assert_true(isInterpolation(lambda));
		return lambda(0) * v0 + 
		       lambda(1) * v1 + 
		       lambda(2) * v2;
	}
	
	
	/**
	 * Quadratic interpolation in plain non-degenerate triangle
	 * @param c_i points
	 * @param v_i values
	 * @param g_i gradients
	 * @param q point to interpolate
	 */
	static TValue interpolate(
			const Real2& c0, const TValue v0, const Gradient g0,
			const Real2& c1, const TValue v1, const Gradient g1,
			const Real2& c2, const TValue v2, const Gradient g2,
			const Real2& q) {
		Real3 lambda = linal::barycentricCoordinates(c0, c1, c2, q);
		assert_true(isInterpolation(lambda));
		return lambda(0) * (v0 + linal::dotProduct(g0, q - c0) / 2.0) +
		       lambda(1) * (v1 + linal::dotProduct(g1, q - c1) / 2.0) +
		       lambda(2) * (v2 + linal::dotProduct(g2, q - c2) / 2.0);
	}
	
	
	/**
	 * Quadratic interpolation in plain non-degenerate triangle,
	 * limited with min-max limiter to avoid oscillations
	 * @param c_i points
	 * @param v_i values
	 * @param g_i gradients
	 * @param q point to interpolate
	 */
	static TValue minMaxInterpolate(
			const Real2& c0, const TValue v0, const Gradient g0,
			const Real2& c1, const TValue v1, const Gradient g1,
			const Real2& c2, const TValue v2, const Gradient g2,
			const Real2& q) {
		TValue quadratic = interpolate(
				c0, v0, g0, c1, v1, g1, c2, v2, g2, q);
		return linal::limiterMinMax(quadratic, v0, v1, v2);
	}
	
	
	/**
	 * Quadratic interpolation in plain non-degenerate triangle,
	 * limited with hybrid limiter to avoid oscillations.
	 * "Hybrid" means if an ocsillation occurs,
	 * interpolation is switched to first-order.
	 * @param c_i points
	 * @param v_i values
	 * @param g_i gradients
	 * @param q point to interpolate
	 */
	static TValue hybridInterpolate(
			const Real2& c0, const TValue v0, const Gradient g0,
			const Real2& c1, const TValue v1, const Gradient g1,
			const Real2& c2, const TValue v2, const Gradient g2,
			const Real2& q) {
		TValue quadratic = interpolate(
				c0, v0, g0, c1, v1, g1, c2, v2, g2, q);
		TValue minMaxLimited = linal::limiterMinMax(quadratic, v0, v1, v2);
		return (quadratic == minMaxLimited) ?
				quadratic : interpolate(c0, v0, c1, v1, c2, v2, q);
	}
	
	
	/**
	 * Given with 4 point-value pairs, determine among them triangle that
	 * contains the query point inside and perform linear interpolation in it
	 * @param c_i and v_i - points and values
	 * @param q point to interpolate
	 */
	static TValue interpolateInOwner(const Real2& c0, const TValue v0,
	                                 const Real2& c1, const TValue v1,
	                                 const Real2& c2, const TValue v2,
	                                 const Real2& c3, const TValue v3,
	                                 const Real2& q) {
		Real3 lambda;
		
		#define TRY_TRIANGLE(a, b, d) \
			lambda = linal::barycentricCoordinates(c##a, c##b, c##d, q); \
			if (isInterpolation(lambda)) { \
				return lambda(0) * v##a + lambda(1) * v##b + lambda(2) * v##d; \
			}
		
		TRY_TRIANGLE(0, 1, 2)
		TRY_TRIANGLE(0, 1, 3)
		TRY_TRIANGLE(0, 2, 3)
		TRY_TRIANGLE(1, 2, 3)
		
		#undef TRY_TRIANGLE
		
		THROW_INVALID_ARG("Containing triangle is not found");
	}
};


}

#endif // LIBGCM_TRIANGLEINTERPOLATOR_HPP
