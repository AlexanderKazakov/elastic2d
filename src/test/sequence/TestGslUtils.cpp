#include <gtest/gtest.h>

#include <algorithm>
#include <libgcm/util/Utils.hpp>
#include <libgcm/util/math/GslUtils.hpp>
#include <libgcm/linal/linal.hpp>

using namespace gcm;


TEST(GslUtils, solveThirdOrderPolynomial) {
	ASSERT_EQ(Real3({0, 0, 0}), gsl_utils::solveThirdOrderPolynomial({0, 0, 0}));
	ASSERT_EQ(Real3({1, 1, 1}), gsl_utils::solveThirdOrderPolynomial({-3, 3, -1}));
	ASSERT_EQ(Real3({-5, 1, 1}), gsl_utils::solveThirdOrderPolynomial({3, -9, 5}));
	
	const real bound = 1e+12;
	const int ITER_NUM = 100000;
	Utils::seedRand();
	for (int i = 0; i < ITER_NUM; i++) {
		real a = Utils::randomReal(-bound, bound);
		real b = Utils::randomReal(-bound, bound);
		real c = Utils::randomReal(-bound, bound);
		real p1 = a + b + c, p2 =  a * b + a * c + b * c, p3 =  a * b * c;
		
		auto x = gsl_utils::solveThirdOrderPolynomial({p1, p2, p3});
		
		std::vector<real> expected = {-a, -b, -c};
		std::sort(expected.begin(), expected.end());
		std::vector<real> actual = {x(0), x(1), x(2)};
		std::sort(actual.begin(), actual.end());
		
		for (size_t k = 0; k < 3; k++) {
			ASSERT_NEAR(expected[k], actual[k], fabs(expected[k]) * gsl_utils::eps)
					<< expected[0] << "\t" << expected[1] << "\t" << expected[2] << "\n"
					<< actual[0] << "\t" << actual[1] << "\t" << actual[2] << "\n";
		}
	}
	
	Utils::seedRand();
	for (int i = 0; i < ITER_NUM; i++) {
		real a = Utils::randomReal(-bound, bound);
		real b = Utils::randomReal(-bound, bound);
		real p1 = a + 2 * b, p2 = 2 * a * b + b*b, p3 = a * b*b;
		
		auto x = gsl_utils::solveThirdOrderPolynomial({p1, p2, p3});
		
		ASSERT_EQ(x(1), x(2))
					<< x(0) << "\t" << x(1) << "\t" << x(2) << "\n";
		ASSERT_NEAR(-a, x(0), fabs(a) * gsl_utils::eps)
					<< -a << "\t" << -b << "\t" << -b << "\n"
					<< x(0) << "\t" << x(1) << "\t" << x(2) << "\n";
		ASSERT_NEAR(-b, x(1), fabs(a) * gsl_utils::eps)
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
		auto m1 = gsl_utils::invert(m);
		auto u = m * m1;
		for (int p = 0; p < 9; p++) {
			for (int q = 0; q < 9; q++) {
				ASSERT_NEAR((p == q), u(p, q), gsl_utils::eps);
			}
		}
	}
}


TEST(GslUtils, determinant) {
	ASSERT_EQ( 0, gsl_utils::determinant( linal::Matrix<5, 5>::Zeros()));
	ASSERT_EQ( 0, gsl_utils::determinant( linal::Matrix<6, 6>::Ones()));
	ASSERT_EQ( 1, gsl_utils::determinant( linal::Matrix<9, 9>::Identity()));
	ASSERT_EQ(-1, gsl_utils::determinant(-linal::Matrix<7, 7>::Identity()));
	
	const int N = 9;
	linal::Matrix<N, N> V; // Vandermonde matrix
	const linal::Vector<N> a = linal::random<linal::Vector<N>>();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			V(i, j) = std::pow(a(i), real(j));
		}
	}
	real det = 1;
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			det *= a(j) - a(i);
		}
	}
	ASSERT_NEAR(det, gsl_utils::determinant(V), std::fabs(det) * gsl_utils::eps);
}


TEST(GslUtils, solveLinearSystem) {
	ASSERT_EQ(Real4::Ones(), gsl_utils::solveLinearSystem(
			5 * linal::Matrix<4, 4>::Identity(), 5 * Real4::Ones()));
	ASSERT_THROW(gsl_utils::solveLinearSystem(
			linal::Matrix<5, 5>::Ones(), linal::Vector<5>::Ones()), Exception);
	
	const real bound = 1e+12;
	const int ITER_NUM = 1000;
	Utils::seedRand();
	for (int i = 0; i < ITER_NUM; i++) {
		auto A = linal::random<linal::Matrix<6, 6>>(-bound, bound);
		auto b = linal::random<linal::Vector<6>>(-bound, bound);
		auto x = gsl_utils::solveLinearSystem(A, b);
		ASSERT_TRUE(linal::approximatelyEqual(A * x, b));
	}
}



