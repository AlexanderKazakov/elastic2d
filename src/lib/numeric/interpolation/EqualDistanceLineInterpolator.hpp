#ifndef LIBGCM_EQUALDISTANCELINEINTERPOLATOR_HPP
#define LIBGCM_EQUALDISTANCELINEINTERPOLATOR_HPP

#include <lib/linal/linal.hpp>
#include <lib/util/Types.hpp>


namespace gcm {
template<class TVector>
class EqualDistanceLineInterpolator {
public:
	/**
	 * Interpolation with minmax limiter
	 * @param src known values at equal distances, @warning it will be overwritten(!)
	 * @param q relative distance between first source point and point to interpolate
	 * @return interpolated value
	 */
	static TVector minMaxInterpolate(std::vector<TVector>& src, const real& q) {
		/// where is the point to interpolate
		size_t k = (size_t) q;
		/// check that perform interpolation, not extrapolation
		assert_le(k, src.size() - 1);
		assert_ge(q, 0);
		
		/// yield values for limiter
		TVector maximum, minimum;
		for (int i = 0; i < TVector::M; i++) {
			maximum(i) = fmax(src[k](i), src[k + 1](i));
			minimum(i) = fmin(src[k](i), src[k + 1](i));
		}
		
		TVector ans = interpolate(src, q);
		
		/// minmax limiter
		for (int i = 0; i < TVector::M; i++) {
			if (ans(i) > maximum(i)) {
				ans(i) = maximum(i);
			} else if (ans(i) < minimum(i)) {
				ans(i) = minimum(i);
			}
		}
		return ans;
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
	static TVector interpolate(std::vector<TVector>& src, const real& q) {
		/// Newton interpolation
		TVector ans = src[0];
		const int p = (int)src.size() - 1; ///< order of interpolation
		for (int i = 1; i <= p; i++) {
			for (int j = 0; j < p - i + 1; j++) {
				src[(size_t)j] =
				        /* recurrent finite differences */
				        (src[(size_t)j + 1] - src[(size_t)j]) *
				        /* coordinate and factorial */
				        ((q - i + 1) / i);
			}
			ans += src[0];
		}
		return ans;
	}

};


}

#endif // LIBGCM_EQUALDISTANCELINEINTERPOLATOR_HPP
