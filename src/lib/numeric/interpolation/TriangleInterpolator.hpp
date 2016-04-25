#ifndef LIBGCM_TRIANGLEINTERPOLATOR_HPP
#define LIBGCM_TRIANGLEINTERPOLATOR_HPP

#include <lib/linal/linal.hpp>


namespace gcm {
template<class TValue>
class TriangleInterpolator {
public:
	typedef linal::VECTOR<2, TValue>              Gradient;
	typedef linal::SYMMETRIC_MATRIX<2, TValue>    Hessian;

	/**
	 * Computes barycentric coordinates of point q in triangle {a, b, c}
	 * The order is {lambda_a, lambda_b, lambda_c}
	 */
	static Real3 barycentricCoordinates(const Real2& a,
	                                    const Real2& b,
	                                    const Real2& c,
	                                    const Real2& q) {
	                                    
		/// @see https://en.wikipedia.org/wiki/Barycentric_coordinate_system 
		linal::Matrix22 T = {a(0) - c(0), b(0) - c(0),
		                     a(1) - c(1), b(1) - c(1)};
		Real2 lambda = linal::solveLinearSystem(T, q - c);
		
		return {lambda(0), lambda(1), 1 - lambda(0) - lambda(1)};
	}
	

	/**
	 * Linear interpolation in plain triangle
	 * @param c_i and v_i - points and values
	 * @param q point to interpolate
	 */
	static TValue interpolate(const Real2& c0, const TValue v0,
	                          const Real2& c1, const TValue v1,
	                          const Real2& c2, const TValue v2,
	                          const Real2& q) {

		Real3 lambda = barycentricCoordinates(c0, c1, c2, q);
		return lambda(0) * v0 + 
		       lambda(1) * v1 + 
		       lambda(2) * v2;		
	}
	
	
	/**
	 * Quadratic interpolation in plain triangle
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

		Real3 lambda = barycentricCoordinates(c0, c1, c2, q);
		return lambda(0) * (v0 + linal::dotProduct(g0, q - c0) / 2.0) +
		       lambda(1) * (v1 + linal::dotProduct(g1, q - c1) / 2.0) +
		       lambda(2) * (v2 + linal::dotProduct(g2, q - c2) / 2.0);
	}
	
	
};


}

#endif // LIBGCM_TRIANGLEINTERPOLATOR_HPP
