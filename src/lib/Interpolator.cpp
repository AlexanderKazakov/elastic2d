#include "Interpolator.hpp"

#include <cmath>

void Interpolator::minMaxInterpolate(Vector &res, std::vector<Vector> &src, const real &q) const {
	// where is the point to interpolate
	uint k = (uint) q;
	// check that we perform interpolation, not extrapolation
	if (k > src.size() - 1 || q < 0) throw "The point to interpolate is out of known area";
	// yield values for limiter
	Vector maximum, minimum;
	for (uint i = 0; i < N; i++) {
		maximum(i) = fmax(src[k].get(i), src[k + 1].get(i));
		minimum(i) = fmin(src[k].get(i), src[k + 1].get(i));
	}
	// interpolate
	interpolate(res, src, q);
	// minmax limiter
	for (uint i = 0; i < N; i++) {
		if (res.get(i) > maximum.get(i)) {
			res(i) = maximum.get(i);
		} else if (res.get(i) < minimum.get(i)) {
			res(i) = minimum.get(i);
		}
	}
}


void Interpolator::interpolate(Vector &res, std::vector<Vector> &src, const real &q) const {
	// Newton interpolation
	res = src[0];
	const uint p = src.size() - 1; // order of interpolation
	for (uint i = 1; i <= p; i++) {
		for (uint j = 0; j < p - i + 1; j++) {
			// TODO - make all this linal operations faster
			src[j] = (src[j + 1] - src[j]) /* recurrent finite differences */
			         * ((q - i + 1) / i) /* coordinate and factorial */ ;
		}
		res += src[0];
	}
}
