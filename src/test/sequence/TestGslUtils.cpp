#include <gtest/gtest.h>

#include <algorithm>
#include <lib/util/Utils.hpp>
#include <lib/util/GslUtils.hpp>
#include <lib/linal/linal.hpp>

using namespace gcm;


TEST(GslUtils, solveThirdOrderPolynomial) {
	ASSERT_EQ(Real3({0, 0, 0}), GslUtils::solveThirdOrderPolynomial({0, 0, 0}));
	ASSERT_EQ(Real3({1, 1, 1}), GslUtils::solveThirdOrderPolynomial({-3, 3, -1}));
	ASSERT_EQ(Real3({-5, 1, 1}), GslUtils::solveThirdOrderPolynomial({3, -9, 5}));
	
	const real bound = 1e+12;
	const int ITER_NUM = 100000;
	Utils::seedRand();
	for (int i = 0; i < ITER_NUM; i++) {
		real a = Utils::randomReal(-bound, bound);
		real b = Utils::randomReal(-bound, bound);
		real c = Utils::randomReal(-bound, bound);
		real p1 = a + b + c, p2 =  a * b + a * c + b * c, p3 =  a * b * c;
		
		auto x = GslUtils::solveThirdOrderPolynomial({p1, p2, p3});
		
		std::vector<real> expected = {-a, -b, -c};
		std::sort(expected.begin(), expected.end());
		std::vector<real> actual = {x(0), x(1), x(2)};
		std::sort(actual.begin(), actual.end());
		
		for (size_t k = 0; k < 3; k++) {
			ASSERT_NEAR(expected[k], actual[k], fabs(expected[k]) * GslUtils::eps)
					<< expected[0] << "\t" << expected[1] << "\t" << expected[2] << "\n"
					<< actual[0] << "\t" << actual[1] << "\t" << actual[2] << "\n";
		}
	}
	
	Utils::seedRand();
	for (int i = 0; i < ITER_NUM; i++) {
		real a = Utils::randomReal(-bound, bound);
		real b = Utils::randomReal(-bound, bound);
		real p1 = a + 2 * b, p2 = 2 * a * b + b*b, p3 = a * b*b;
		
		auto x = GslUtils::solveThirdOrderPolynomial({p1, p2, p3});
		
		ASSERT_EQ(x(1), x(2))
					<< x(0) << "\t" << x(1) << "\t" << x(2) << "\n";
		ASSERT_NEAR(-a, x(0), fabs(a) * GslUtils::eps)
					<< -a << "\t" << -b << "\t" << -b << "\n"
					<< x(0) << "\t" << x(1) << "\t" << x(2) << "\n";
		ASSERT_NEAR(-b, x(1), fabs(a) * GslUtils::eps)
					<< -a << "\t" << -b << "\t" << -b << "\n"
					<< x(0) << "\t" << x(1) << "\t" << x(2) << "\n";
	}
}


TEST(GslUtils, invertMatrix) {
	const real bound = 1e+12;
	const int ITER_NUM = 1000;
	Utils::seedRand();
	for (int i = 0; i < ITER_NUM; i++) {
		auto m = linal::random<linal::Matrix<9, 9>>(-bound, bound);
		auto m1 = GslUtils::invert(m);
		auto u = m * m1;
		for (int p = 0; p < 9; p++) {
			for (int q = 0; q < 9; q++) {
				ASSERT_NEAR((p == q), u(p, q), GslUtils::eps);
			}
		}
	}
}



