#ifndef ELASTIC2D_INTERPOLATOR_HPP
#define ELASTIC2D_INTERPOLATOR_HPP

#include <vector>

#include "lib/PDEMatrices.hpp"


class Interpolator {
public:
	/**
	 * Interpolation by %interpolate function with minmax limiter.
	 * @param res value to intepolate
	 * @param src known values at equal distances, will be overwritten(!)
	 * @param q relative distance between first source point and point to interpolate.
	 */
	void minMaxInterpolate(Vector& res, std::vector<Vector>& src, const real& q) const;
	/**
	 * Interpolation by Newton polynomials (forward interpolation).
	 * Interpolate %res by given values %src. The distance between known values is constant,
	 * they are all on the same line.
	 * Point of unknown value is at distance q (relatively to distance between known points)
 	 * from the first point of %src.
 	 * The order of interpolation is equal to (%src.size() - 1).
	 * @param res value to intepolate
	 * @param src known values at equal distances, will be overwritten(!)
	 * @param q relative distance between first source point and point to interpolate.
	 */
	void interpolate(Vector& res, std::vector<Vector>& src, const real& q) const;
};


#endif //ELASTIC2D_INTERPOLATOR_HPP
