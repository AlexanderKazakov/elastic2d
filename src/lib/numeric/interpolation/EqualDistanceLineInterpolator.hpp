#ifndef LIBGCM_EQUALDISTANCELINEINTERPOLATOR_HPP
#define LIBGCM_EQUALDISTANCELINEINTERPOLATOR_HPP

#include <vector>
#include <lib/linal/linal.hpp>
#include <lib/util/Types.hpp>


namespace gcm {
template<class TVector>
class EqualDistanceLineInterpolator {
public:
	/**
	 * Interpolation with minmax limiter.
	 * @param res value to interpolate
	 * @param src known values at equal distances, @warning it will be overwritten(!)
	 * @param q relative distance between first source point and point to interpolate.
	 */
	static void minMaxInterpolate(TVector& res, std::vector<TVector>& src, const real& q) {
		/// where is the point to interpolate
		unsigned long k = (unsigned long) q;
		/// check that we perform interpolation, not extrapolation
		assert_le(k, src.size() - 1);
		assert_ge(q, 0);
		/// yield values for limiter
		TVector maximum, minimum;
		for (int i = 0; i < TVector::M; i++) {
			maximum(i) = fmax(src[k](i), src[k + 1](i));
			minimum(i) = fmin(src[k](i), src[k + 1](i));
		}
		/// interpolate
		interpolate(res, src, q);
		/// minmax limiter
		for (int i = 0; i < TVector::M; i++) {
			if (res(i) > maximum(i)) {
				res(i) = maximum(i);
			} else if (res(i) < minimum(i)) {
				res(i) = minimum(i);
			}
		}
	}

	/**
	 * Interpolation by Newton polynomials (forward interpolation).
	 * Interpolate res by given values src. The distance between known values is constant,
	 * they are all on the same line.
	 * Point of unknown value is at distance q (relatively to distance between known points)
	 * from the first point of src.
	 * The order of interpolation is equal to (src.size() - 1).
	 * @param res value to interpolate
	 * @param src known values at equal distances, @warning it will be overwritten(!)
	 * @param q relative distance between first source point and point to interpolate.
	 */
	static void interpolate(TVector& res, std::vector<TVector>& src, const real& q) {
		/// Newton interpolation
		res = src[0];
		const int p = (int)src.size() - 1; ///< order of interpolation
		for (int i = 1; i <= p; i++) {
			for (int j = 0; j < p - i + 1; j++) {
				// TODO - make all this linal operations faster, replace std::vector
				src[(unsigned long)j] =
				        /* recurrent finite differences */
				        (src[(unsigned long)j + 1] - src[(unsigned long)j]) *
				        /* coordinate and factorial */
				        ((q - i + 1) / i);
			}
			res += src[0];
		}
	}

};


}

#endif // LIBGCM_EQUALDISTANCELINEINTERPOLATOR_HPP
