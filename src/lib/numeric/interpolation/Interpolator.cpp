#include <lib/numeric/interpolation/Interpolator.hpp>
#include <lib/rheology/variables/variables.hpp>

using namespace gcm;

template<class TVector>
void Interpolator<TVector>::minMaxInterpolate(TVector &res, std::vector<TVector> &src, const real &q) const {
	// where is the point to interpolate
	unsigned long k = (unsigned long) q;
	// check that we perform interpolation, not extrapolation
	assert_le(k, src.size() - 1);
	assert_ge(q, 0);
	// yield values for limiter
	TVector maximum, minimum;
	for (int i = 0; i < TVector::M; i++) {
		maximum(i) = fmax(src[k](i), src[k + 1](i));
		minimum(i) = fmin(src[k](i), src[k + 1](i));
	}
	// interpolate
	interpolate(res, src, q);
	// minmax limiter
	for (int i = 0; i < TVector::M; i++) {
		if (res(i) > maximum(i)) {
			res(i) = maximum(i);
		} else if (res(i) < minimum(i)) {
			res(i) = minimum(i);
		}
	}
}

template<class TVector>
void Interpolator<TVector>::interpolate(TVector &res, std::vector<TVector> &src, const real &q) const {
	// Newton interpolation
	res = src[0];
	const int p = (int)src.size() - 1; // order of interpolation
	for (int i = 1; i <= p; i++) {
		for (int j = 0; j < p - i + 1; j++) {
			// TODO - make all this linal operations faster, replace std::vector
			src[(unsigned long)j] = (src[(unsigned long)j + 1] - src[(unsigned long)j]) /* recurrent finite differences */
			         * ((q - i + 1) / i) /* coordinate and factorial */ ;
		}
		res += src[0];
	}
}

template class Interpolator<linal::Vector<2>>;
template class Interpolator<linal::Vector<5>>;
template class Interpolator<linal::Vector<9>>;

template class Interpolator<linal::Vector<2, VelocitySigmaVariables<1>>>;
template class Interpolator<linal::Vector<5, VelocitySigmaVariables<2>>>;
template class Interpolator<linal::Vector<9, VelocitySigmaVariables<3>>>;
