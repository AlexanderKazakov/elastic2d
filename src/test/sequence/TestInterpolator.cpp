#include <gtest/gtest.h>
#include <cmath>

#include <lib/linal/linal.hpp>
#include <lib/numeric/interpolation/interpolation.hpp>

using namespace gcm;
using namespace gcm::linal;


TEST(EqualDistanceLineInterpolator, Const) {
	const int N = 5;
	for (int i = 0; i <= 20; i++) {
		std::vector<Vector<N>> src((unsigned int)(i + 1));
		for (real q = 0; q <= (real) i; q += 0.2) {
			for (auto& item : src) {
				for (int j = 0; j < N; j++) {
					item(j) = sinh(j - (real) (N) / 2);
				}
			}

			auto res = EqualDistanceLineInterpolator<Vector<N>>::interpolate(src, q);
			for (int j = 0; j < N; j++) {
				ASSERT_NEAR(res(j), sinh(j - (real) (N) / 2), EQUALITY_TOLERANCE);
			}
		}
	}
}


TEST(EqualDistanceLineInterpolator, Linear) {
	const int N = 9;
	for (int k = 0; k <= 20; k++) {
		std::vector<Vector<N>> src((unsigned int)(k + 1));
		for (real q = 0; q <= (real) k; q += 0.2) {
			for (int i = 0; i < (int) src.size(); i++) {
				for (int j = 0; j < N; j++) {
					src[(unsigned long)i](j) = (j - (real) (N) / 2) * i + 2 * (j - (real) (N) / 2);
				}
			}

			auto res = EqualDistanceLineInterpolator<Vector<N>>::interpolate(src, q);
			for (int j = 0; j < N; j++) {
				ASSERT_NEAR(res(j), (j - (real) (N) / 2) * q + 2 * (j - (real) (N) / 2),
				            EQUALITY_TOLERANCE);
			}
		}
	}
}


TEST(EqualDistanceLineInterpolator, Quadratic) {
	const int N = 9;
	for (real q = 0.0; q <= 2.0; q += 0.1) {
		std::vector<Vector<N>> src(3);
		src[0](0) = 0; src[1](0) = 1; src[2](0) = 4;
		src[0](1) = 0; src[1](1) = -1; src[2](1) = -4;
		src[0](2) = -3; src[1](2) = 15; src[2](2) = 89;
		auto res = EqualDistanceLineInterpolator<Vector<N>>::interpolate(src, q);
		ASSERT_NEAR(res(0), q * q, EQUALITY_TOLERANCE) << q;
		ASSERT_NEAR(res(1), -q * q, EQUALITY_TOLERANCE) << q;
		ASSERT_NEAR(res(2), 7 * (2 * q) * (2 * q) - 5 * (2 * q) - 3,
		            EQUALITY_TOLERANCE) << q;
	}
}


TEST(EqualDistanceLineInterpolator, MinMax) {
	const int N = 2;
	std::vector<Vector<N>> src(3);
	src[0](0) = -9.0; src[1](0) = -1.0; src[2](0) = -1.0;
	src[0](1) = 9.0; src[1](1) = 1.0; src[2](1) = 1.0;
	real q = 1.5;
	auto res = EqualDistanceLineInterpolator<Vector<N>>::minMaxInterpolate(src, q);
	ASSERT_EQ(res(0), -1.0);
	ASSERT_EQ(res(1), 1.0);
}


TEST(EqualDistanceLineInterpolator, MinMaxFifthOrder) {
	const int N = 2;
	// minmax interpolation
	for (real q = 0.0; q <= 5.0; q += 0.1) {
		std::vector<Vector<N>> src(6);
		src[0](0) = 0; src[1](0) = 0; src[2](0) = 0;
		src[3](0) = 1; src[4](0) = 1; src[5](0) = 1;

		auto res = EqualDistanceLineInterpolator<Vector<N>>::minMaxInterpolate(src, q);
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


TEST(EqualDistanceLineInterpolator, TenthOrder) {
	const int N = 2;
	auto func = [] (real x) {
		return x*x*x*x*x*x*x*x*x*x + 4 * x*x*x*x*x*x*x*x*x -
		       25 * x*x*x*x*x*x*x + 2 *x*x*x*x*x - 3 *x*x + 5;
	};

	std::vector<Vector<N>> src(11);
	for (real q = 0; q < 11.0; q += 0.1) {
		for (int i = 0; i < 11; i++) {
			src[(unsigned long)i](0) = func(i);
		}
		auto res = EqualDistanceLineInterpolator<Vector<N>>::interpolate(src, q);
		ASSERT_NEAR(res(0), func(q), fabs(res(0)) * EQUALITY_TOLERANCE * 1e+3);
	}
}


TEST(EqualDistanceLineInterpolator, Exceptions) {
	const int N = 5;
	std::vector<Vector<N>> src(2);
	ASSERT_THROW(EqualDistanceLineInterpolator<Vector<N>>::
			minMaxInterpolate(src, -0.3), Exception);
	ASSERT_THROW(EqualDistanceLineInterpolator<Vector<N>>::
			minMaxInterpolate(src, 2.5), Exception);
}


TEST(TriangleInterpolator, linear) {
	auto f = [](Real2 x) { return 5*x(0) + 8*x(1) - 2; };
	
	Utils::seedRand();
	for (int i = 0; i < 1000; i++) {
		Real2 a = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real2 b = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real2 c = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real2 q = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		
		ASSERT_NEAR(f(q), 
		            TriangleInterpolator<real>::interpolate(a, f(a),
		                                                    b, f(b),
		                                                    c, f(c),
		                                                    q),
		            EQUALITY_TOLERANCE * fabs(f(q)));
	}
}


TEST(TriangleInterpolator, quadratic) {
	auto f = [](Real2 x) {
		return 8*x(0)*x(0) + 10*x(0)*x(1) - 15*x(1)*x(1) + 5*x(0) + 8*x(1) - 2;
	};
	auto g = [](Real2 x) { // = Df / D(x, y)
		return linal::VECTOR<2, real>({ 16*x(0) + 10*x(1) + 5,
		                               -30*x(1) + 10*x(0) + 8});
	};
		
	Utils::seedRand();
	for (int i = 0; i < 1000; i++) {
		Real2 a = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real2 b = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real2 c = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real2 q = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		
		ASSERT_NEAR(f(q), 
		            TriangleInterpolator<real>::interpolate(a, f(a), g(a),
		                                                    b, f(b), g(b),
		                                                    c, f(c), g(c),
		                                                    q),
		            EQUALITY_TOLERANCE * fabs(f(q)));
	}
}


TEST(TriangleInterpolator, interpolateInOwner) {
	ASSERT_EQ(1, TriangleInterpolator<real>::interpolateInOwner(
			{0, 0}, 1,
			{0, 1}, 1,
			{1, 0}, 1,
			{1, 1}, 1e+100,
		    {0.2, 0.2}));
}


TEST(TetrahedronInterpolator, linear) {
	auto f = [](Real3 x) { return 5*x(0) + 8*x(1) - 4*x(2) - 2; };
	
	Utils::seedRand();
	for (int i = 0; i < 1000; i++) {
		Real3 a = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real3 b = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real3 c = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real3 d = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real3 q = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		
		ASSERT_NEAR(f(q), 
		            TetrahedronInterpolator<real>::interpolate(a, f(a),
		                                                       b, f(b),
		                                                       c, f(c),
		                                                       d, f(d),
		                                                       q),
		            EQUALITY_TOLERANCE * fabs(f(q)));
	}
}


TEST(TetrahedronInterpolator, quadratic) {
	auto f = [](Real3 x) {
		return  8*x(0)*x(0) + 10*x(0)*x(1) - 9*x(0)*x(2) -
		       15*x(1)*x(1) + 6*x(1)*x(2) - 7*x(2)*x(2) + 
		        5*x(0) + 8*x(1) - 7*x(2) - 2;
	};
	auto g = [](Real3 x) { // = Df / D(x, y)
		return linal::VECTOR<3, real>({ 16*x(0) + 10*x(1) - 9*x(2) + 5,
		                               -30*x(1) + 10*x(0) + 6*x(2) + 8,
		                               -14*x(2) + 6*x(1)  - 9*x(0) - 7});
	};
		
	Utils::seedRand();
	for (int i = 0; i < 1000; i++) {
		Real3 a = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real3 b = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real3 c = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real3 d = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		Real3 q = {Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6), Utils::randomReal(-1e+6, 1e+6)};
		
		ASSERT_NEAR(f(q), 
		            TetrahedronInterpolator<real>::interpolate(a, f(a), g(a),
		                                                       b, f(b), g(b),
		                                                       c, f(c), g(c),
		                                                       d, f(d), g(d),
		                                                       q),
		            EQUALITY_TOLERANCE * fabs(f(q)));
	}
}


TEST(TetrahedronInterpolator, interpolateInOwner) {
	ASSERT_EQ(1, TetrahedronInterpolator<real>::interpolateInOwner(
			{0, 0, 0}, 1,
			{0, 1, 0}, 1,
			{1, 0, 0}, 1,
			{0, 0, 1}, 1,
			{0, 1, 1}, 1e+100,
			{1, 0, 1}, 1e+100,
		    {0.1, 0.1, 0.1}));
}





