#include <gtest/gtest.h>

#include "lib/gcm_matrices/IdealElastic3DGcmMatrices.hpp"
#include "lib/gcm_matrices/IdealElastic2DGcmMatrices.hpp"
#include "lib/gcm_matrices/IdealElastic1DGcmMatrices.hpp"

using namespace gcm;
using namespace gcm::linal;

template<class TGcmMatrices>
class TestGcmMatrices : public testing::Test {

	const int NUMBER_ITERATIONS = 1000;

	const real RHO_MAX = 100.0;
	const real RHO_MIN = 0.01;
	const real LAMBDA_MAX = 1e+6;
	const real LAMBDA_MIN = 1.0;
	const real MU_MAX = 1e+6;
	const real MU_MIN = 1.0;

	const real rho0 = 8.0;
	const real lambda0 = 12e+4;
	const real mu0 = 77e+3;

	void testTraces(const TGcmMatrices& matrix) {
		for (int i = 0; i < matrix.DIMENSIONALITY; i++) {
			ASSERT_NEAR(matrix.A(i).A.trace(), matrix.A(i).L.trace(), EQUALITY_TOLERANCE)
										<< "(" << i << ") A = " << matrix.A(i).A << "L = " << matrix.A(i).L;
		}
	};

	void testLeftEigenVectors(const TGcmMatrices& matrix) {
		for (int s = 0; s < matrix.DIMENSIONALITY; s++) {
			Matrix<TGcmMatrices::M,TGcmMatrices::M> AU1 = matrix.A(s).A * matrix.A(s).U1;
			Matrix<TGcmMatrices::M,TGcmMatrices::M> U1L = matrix.A(s).U1 * matrix.A(s).L;
			for (int i = 0; i < TGcmMatrices::M; i++) {
				for (int j = 0; j < TGcmMatrices::M; j++) {
					ASSERT_NEAR(AU1(i, j), U1L(i, j), EQUALITY_TOLERANCE) << "(" << s << ") U1 = " << matrix.A(s).U1;
				}
			}
		}
	};

	void testRightEigenVectors(const TGcmMatrices& matrix) {
		for (int s = 0; s < matrix.DIMENSIONALITY; s++) {
			Matrix<TGcmMatrices::M,TGcmMatrices::M> UA = matrix.A(s).U * matrix.A(s).A;
			Matrix<TGcmMatrices::M,TGcmMatrices::M> LU = matrix.A(s).L * matrix.A(s).U;
			for(int i = 0; i < TGcmMatrices::M; i++) {
				for(int j = 0; j < TGcmMatrices::M; j++) {
					ASSERT_NEAR(UA(i, j), LU(i, j), EQUALITY_TOLERANCE) << "(" << s << ") U = " << matrix.A(s).U;
				}
			}
		}
	};

	void testInverseMatrix(const TGcmMatrices& matrix) {
		for (int s = 0; s < matrix.DIMENSIONALITY; s++) {
			Matrix<TGcmMatrices::M,TGcmMatrices::M> UU1 = matrix.A(s).U * matrix.A(s).U1;
			for(int i = 0; i < TGcmMatrices::M; i++) {
				for(int j = 0; j < TGcmMatrices::M; j++) {
					ASSERT_NEAR(UU1(i, j), (i == j), EQUALITY_TOLERANCE) << "(" << s << ") UU1 = " << UU1;
				}
			}
		}
	};

protected:
	void testDiagonalization() {
		testTraces(TGcmMatrices(rho0, lambda0, mu0));
		testLeftEigenVectors(TGcmMatrices(rho0, lambda0, mu0));
		testRightEigenVectors(TGcmMatrices(rho0, lambda0, mu0));
		testInverseMatrix(TGcmMatrices(rho0, lambda0, mu0));

		srand(time(0));
		for (int i = 0; i < NUMBER_ITERATIONS; i++) {
			real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
			real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
			real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;
			testTraces(TGcmMatrices(rho, lambda, mu));
			testLeftEigenVectors(TGcmMatrices(rho, lambda, mu));
			testRightEigenVectors(TGcmMatrices(rho, lambda, mu));
			testInverseMatrix(TGcmMatrices(rho, lambda, mu));
		}
	}
};


/** Look at https://github.com/google/googletest/blob/master/googletest/samples/sample6_unittest.cc for explaination */
#if GTEST_HAS_TYPED_TEST_P
using testing::Types;
TYPED_TEST_CASE_P(TestGcmMatrices);

TYPED_TEST_P(TestGcmMatrices, Diagonalization) {
	this->testDiagonalization();
}

REGISTER_TYPED_TEST_CASE_P(TestGcmMatrices, Diagonalization);

// write in generics all the GcmMatrices implementations using in mpi connections
typedef Types<IdealElastic1DGcmMatrices, IdealElastic2DGcmMatrices, IdealElastic3DGcmMatrices> AllImplementations;

INSTANTIATE_TYPED_TEST_CASE_P(AllGcmMatrices, TestGcmMatrices, AllImplementations);
#endif // GTEST_HAS_TYPED_TEST_P