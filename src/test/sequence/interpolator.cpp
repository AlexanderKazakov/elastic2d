#include <gtest/gtest.h>
#include <cmath>

#include "lib/linal/Vector.hpp"
#include "lib/interpolation/Interpolator.hpp"

#define MAX_ACCURACY_ORDER 20

using namespace gcm;
using namespace gcm::linal;

TEST(Interpolator, Const)
{
	const int N = 5;
	Interpolator<Vector<N>> interpolator;
	for (int i = 0; i <= MAX_ACCURACY_ORDER; i++) {
		Vector<N> res;
		std::vector<Vector<N>> src(i + 1);
		for(real q = 0; q <= (real) i; q+= 0.2) {
			for (auto &item : src) {
				for (int j = 0; j < N; j++) {
					item(j) = sinh(j - (real) (N) / 2);
				}
			}

			interpolator.interpolate(res, src, q);
			for (int j = 0; j < N; j++) {
				ASSERT_NEAR(res(j), sinh(j - (real) (N) / 2), EQUALITY_TOLERANCE);
			}
		}
	}
}


TEST(Interpolator, Linear)
{
	const int N = 9;
	Interpolator<Vector<N>> interpolator;
	for (int k = 0; k <= MAX_ACCURACY_ORDER; k++) {
		Vector<N> res;
		std::vector<Vector<N>> src(k + 1);
		for(real q = 0; q <= (real) k; q+= 0.2) {
			for (int i = 0; i < src.size(); i++) {
				for (int j = 0; j < N; j++) {
					src[i](j) = (j - (real) (N) / 2) * i + 2 * (j - (real) (N) / 2);
				}
			}

			interpolator.interpolate(res, src, q);
			for (int j = 0; j < N; j++) {
				ASSERT_NEAR(res(j), (j - (real) (N) / 2) * q + 2 * (j - (real) (N) / 2), EQUALITY_TOLERANCE);
			}
		}
	}
}


TEST(Interpolator, Quadratic)
{
	const int N = 9;
	Interpolator<Vector<N>> interpolator;

	for(real q = 0.0; q <= 2.0; q+= 0.1) {
		Vector<N> res;
		std::vector<Vector<N>> src(3);
		src[0](0) = 0; src[1](0) = 1; src[2](0) = 4;
		src[0](1) = 0; src[1](1) = -1; src[2](1) = -4;
		src[0](2) = -3; src[1](2) = 15; src[2](2) = 89;
		interpolator.interpolate(res, src, q);
		ASSERT_NEAR(res(0), q*q, EQUALITY_TOLERANCE) << q;
		ASSERT_NEAR(res(1), -q*q, EQUALITY_TOLERANCE) << q;
		ASSERT_NEAR(res(2), 7*(2*q)*(2*q) - 5*(2*q) - 3, EQUALITY_TOLERANCE) << q;
	}
}


TEST(Interpolator, MinMax)
{
	const int N = 2;
	Interpolator<Vector<N>> interpolator;
	Vector<N> res;
	std::vector<Vector<N>> src(3);
	src[0](0) = -9.0; src[1](0) = -1.0; src[2](0) = -1.0;
	src[0](1) = 9.0; src[1](1) = 1.0; src[2](1) = 1.0;
	real q = 1.5;
	interpolator.minMaxInterpolate(res, src, q);
	ASSERT_EQ(res(0), -1.0);
	ASSERT_EQ(res(1), 1.0);
}


TEST(Interpolator, MinMaxFifthOrder)
{
	const int N = 2;
	Interpolator<Vector<N>> interpolator;

	// minmax interpolation
	for (real q = 0.0; q <= 5.0; q+= 0.1) {
		Vector<N> res;
		std::vector<Vector<N>> src(6);
		src[0](0) = 0; src[1](0) = 0; src[2](0) = 0;
		src[3](0) = 1; src[4](0) = 1; src[5](0) = 1;

		interpolator.minMaxInterpolate(res, src, q);
		if (q <= 2) {
			ASSERT_NEAR(res(0), 0, EQUALITY_TOLERANCE) << q;
		} else if (q < 3) {
			ASSERT_LT(res(0), 1);
			ASSERT_GT(res(0), 0);
		} else {
			ASSERT_NEAR(res(0), 1, EQUALITY_TOLERANCE) << q;
		}
	}
}


TEST(Interpolator, TenthOrder)
{
	const int N = 2;
	auto func = [](real x) {
		return x * x * x * x * x * x * x * x * x * x + 4 * x * x * x * x * x * x * x * x * x -
		       25 * x * x * x * x * x * x * x + 2 * x * x * x * x * x - 3 * x * x + 5;
	};
	Interpolator<Vector<N>> interpolator;
	Vector<N> res;
	std::vector<Vector<N>> src(11);
	for (real q = 0; q < 11.0; q += 0.1) {
		for (int i = 0; i < 11; i++) {
			real x = i - 5;
			src[i](0) = func(i);
		}
		interpolator.interpolate(res, src, q);
		ASSERT_NEAR(res(0), func(q), fabs(res(0)) * EQUALITY_TOLERANCE * 1e+3);
	}
}


TEST(Interpolator, Exceptions)
{
	const int N = 5;
	Interpolator<Vector<N>> interpolator;
	Vector<N> res;
	std::vector<Vector<N>> src(2);
	ASSERT_ANY_THROW(interpolator.minMaxInterpolate(res, src, -0.3));
	ASSERT_ANY_THROW(interpolator.minMaxInterpolate(res, src, 2.5));
}