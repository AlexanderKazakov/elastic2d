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

	/** Check that order of values in node and matrices and order of indices in WAVE_COLUMNS is concerted */
	void testValuesOrdersConcerted(const TGcmMatrices& matrix);

protected:
	void testDiagonalization() {
		testTraces(TGcmMatrices(rho0, lambda0, mu0));
		testLeftEigenVectors(TGcmMatrices(rho0, lambda0, mu0));
		testRightEigenVectors(TGcmMatrices(rho0, lambda0, mu0));
		testInverseMatrix(TGcmMatrices(rho0, lambda0, mu0));
		testValuesOrdersConcerted(TGcmMatrices(rho0, lambda0, mu0));

		srand(time(0));
		for (int i = 0; i < NUMBER_ITERATIONS; i++) {
			real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
			real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
			real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;

			testTraces(TGcmMatrices(rho, lambda, mu));
			testLeftEigenVectors(TGcmMatrices(rho, lambda, mu));
			testRightEigenVectors(TGcmMatrices(rho, lambda, mu));
			testInverseMatrix(TGcmMatrices(rho, lambda, mu));
			testValuesOrdersConcerted(TGcmMatrices(rho, lambda, mu));
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


template<>
void TestGcmMatrices<IdealElastic1DGcmMatrices>::testValuesOrdersConcerted(const IdealElastic1DGcmMatrices &matrix) {
	// Vx 0, Sxx 0
	
	int indexForward = IdealElastic1DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::P_FORWARD);
	int indexBackward = IdealElastic1DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::P_BACKWARD);

	linal::Vector<IdealElastic1DGcmMatrices::M> pForward = matrix.A(0).U1.getColumn(indexForward);
	ASSERT_LT(pForward(0) /* Vx */ * pForward(1) /* Sxx */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(0).L(indexForward, indexForward), 0);

	linal::Vector<IdealElastic1DGcmMatrices::M> pBackward = matrix.A(0).U1.getColumn(indexBackward);
	ASSERT_GT(pBackward(0) /* Vx */ * pBackward(1) /* Sxx */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(0).L(indexBackward, indexBackward), 0);
};

template<>
void TestGcmMatrices<IdealElastic2DGcmMatrices>::testValuesOrdersConcerted(const IdealElastic2DGcmMatrices &matrix) {
	// Vx 0, Vy 1, Sxx 2, Sxy 3, Syy 4

	int indexForwardP = IdealElastic2DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::P_FORWARD);
	int indexBackwardP = IdealElastic2DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::P_BACKWARD);
	int indexForwardS = IdealElastic2DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::S1_FORWARD);
	int indexBackwardS = IdealElastic2DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::S1_BACKWARD);

	linal::Vector<IdealElastic2DGcmMatrices::M> pForwardX = matrix.A(0).U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardX(0) /* Vx */ * pForwardX(2) /* Sxx */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(0).L(indexForwardP, indexForwardP), 0);

	linal::Vector<IdealElastic2DGcmMatrices::M> pBackwardX = matrix.A(0).U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardX(0) /* Vx */ * pBackwardX(2) /* Sxx */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(0).L(indexBackwardP, indexBackwardP), 0);


	linal::Vector<IdealElastic2DGcmMatrices::M> sForwardX = matrix.A(0).U1.getColumn(indexForwardS);
	ASSERT_LT(sForwardX(1) /* Vy */ * sForwardX(3) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.A(0).L(indexForwardS, indexForwardS), 0);

	linal::Vector<IdealElastic2DGcmMatrices::M> sBackwardX = matrix.A(0).U1.getColumn(indexBackwardS);
	ASSERT_GT(sBackwardX(1) /* Vy */ * sBackwardX(3) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.A(0).L(indexBackwardS, indexBackwardS), 0);



	linal::Vector<IdealElastic2DGcmMatrices::M> pForwardY = matrix.A(1).U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardY(1) /* Vy */ * pForwardY(4) /* Syy */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(1).L(indexForwardP, indexForwardP), 0);

	linal::Vector<IdealElastic2DGcmMatrices::M> pBackwardY = matrix.A(1).U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardY(1) /* Vy */ * pBackwardY(4) /* Syy */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(1).L(indexBackwardP, indexBackwardP), 0);


	linal::Vector<IdealElastic2DGcmMatrices::M> sForwardY = matrix.A(1).U1.getColumn(indexForwardS);
	ASSERT_LT(sForwardY(0) /* Vx */ * sForwardY(3) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.A(1).L(indexForwardS, indexForwardS), 0);

	linal::Vector<IdealElastic2DGcmMatrices::M> sBackwardY = matrix.A(1).U1.getColumn(indexBackwardS);
	ASSERT_GT(sBackwardY(0) /* Vx */ * sBackwardY(3) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.A(1).L(indexBackwardS, indexBackwardS), 0);
};

template<>
void TestGcmMatrices<IdealElastic3DGcmMatrices>::testValuesOrdersConcerted(const IdealElastic3DGcmMatrices &matrix) {
	// Vx 0, Vy 1, Vz 2, Sxx 3, Sxy 4, Sxz 5, Syy 6, Syz 7, Szz 8

	int indexForwardP = IdealElastic3DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::P_FORWARD);
	int indexBackwardP = IdealElastic3DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::P_BACKWARD);
	int indexForwardS1 = IdealElastic3DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::S1_FORWARD);
	int indexBackwardS1 = IdealElastic3DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::S1_BACKWARD);
	int indexForwardS2 = IdealElastic3DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::S2_FORWARD);
	int indexBackwardS2 = IdealElastic3DGcmMatrices::WAVE_COLUMNS.at(Waves::WAVE::S2_BACKWARD);

	linal::Vector<IdealElastic3DGcmMatrices::M> pForwardX = matrix.A(0).U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardX(0) /* Vx */ * pForwardX(3) /* Sxx */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(0).L(indexForwardP, indexForwardP), 0);

	linal::Vector<IdealElastic3DGcmMatrices::M> pBackwardX = matrix.A(0).U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardX(0) /* Vx */ * pBackwardX(3) /* Sxx */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(0).L(indexBackwardP, indexBackwardP), 0);


	linal::Vector<IdealElastic3DGcmMatrices::M> s1ForwardX = matrix.A(0).U1.getColumn(indexForwardS1);
	ASSERT_LT(s1ForwardX(1) /* Vy */ * s1ForwardX(4) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.A(0).L(indexForwardS1, indexForwardS1), 0);

	linal::Vector<IdealElastic3DGcmMatrices::M> s1BackwardX = matrix.A(0).U1.getColumn(indexBackwardS1);
	ASSERT_GT(s1BackwardX(1) /* Vy */ * s1BackwardX(4) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.A(0).L(indexBackwardS1, indexBackwardS1), 0);


	linal::Vector<IdealElastic3DGcmMatrices::M> s2ForwardX = matrix.A(0).U1.getColumn(indexForwardS2);
	ASSERT_LT(s2ForwardX(2) /* Vz */ * s2ForwardX(5) /* Sxz */, 0); // shear wave
	ASSERT_GT(matrix.A(0).L(indexForwardS2, indexForwardS2), 0);

	linal::Vector<IdealElastic3DGcmMatrices::M> s2BackwardX = matrix.A(0).U1.getColumn(indexBackwardS2);
	ASSERT_GT(s2BackwardX(2) /* Vz */ * s2BackwardX(5) /* Sxz */, 0); // shear wave
	ASSERT_LT(matrix.A(0).L(indexBackwardS2, indexBackwardS2), 0);



	linal::Vector<IdealElastic3DGcmMatrices::M> pForwardY = matrix.A(1).U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardY(1) /* Vy */ * pForwardY(6) /* Syy */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(1).L(indexForwardP, indexForwardP), 0);

	linal::Vector<IdealElastic3DGcmMatrices::M> pBackwardY = matrix.A(1).U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardY(1) /* Vy */ * pBackwardY(6) /* Syy */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(1).L(indexBackwardP, indexBackwardP), 0);


	linal::Vector<IdealElastic3DGcmMatrices::M> s1ForwardY = matrix.A(1).U1.getColumn(indexForwardS1);
	ASSERT_LT(s1ForwardY(0) /* Vx */ * s1ForwardY(4) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.A(1).L(indexForwardS1, indexForwardS1), 0);

	linal::Vector<IdealElastic3DGcmMatrices::M> s1BackwardY = matrix.A(1).U1.getColumn(indexBackwardS1);
	ASSERT_GT(s1BackwardY(0) /* Vx */ * s1BackwardY(4) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.A(1).L(indexBackwardS1, indexBackwardS1), 0);


	linal::Vector<IdealElastic3DGcmMatrices::M> s2ForwardY = matrix.A(1).U1.getColumn(indexForwardS2);
	ASSERT_LT(s2ForwardY(2) /* Vz */ * s2ForwardY(7) /* Syz */, 0); // shear wave
	ASSERT_GT(matrix.A(1).L(indexForwardS2, indexForwardS2), 0);

	linal::Vector<IdealElastic3DGcmMatrices::M> s2BackwardY = matrix.A(1).U1.getColumn(indexBackwardS2);
	ASSERT_GT(s2BackwardY(2) /* Vz */ * s2BackwardY(7) /* Syz */, 0); // shear wave
	ASSERT_LT(matrix.A(1).L(indexBackwardS2, indexBackwardS2), 0);



	linal::Vector<IdealElastic3DGcmMatrices::M> pForwardZ = matrix.A(2).U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardZ(2) /* Vz */ * pForwardZ(8) /* Szz */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(2).L(indexForwardP, indexForwardP), 0);

	linal::Vector<IdealElastic3DGcmMatrices::M> pBackwardZ = matrix.A(2).U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardZ(2) /* Vz */ * pBackwardZ(8) /* Szz */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(2).L(indexBackwardP, indexBackwardP), 0);


	linal::Vector<IdealElastic3DGcmMatrices::M> s1ForwardZ = matrix.A(2).U1.getColumn(indexForwardS1);
	ASSERT_LT(s1ForwardZ(0) /* Vx */ * s1ForwardZ(5) /* Sxz */, 0); // shear wave
	ASSERT_GT(matrix.A(2).L(indexForwardS1, indexForwardS1), 0);

	linal::Vector<IdealElastic3DGcmMatrices::M> s1BackwardZ = matrix.A(2).U1.getColumn(indexBackwardS1);
	ASSERT_GT(s1BackwardZ(0) /* Vx */ * s1BackwardZ(5) /* Sxz */, 0); // shear wave
	ASSERT_LT(matrix.A(2).L(indexBackwardS1, indexBackwardS1), 0);


	linal::Vector<IdealElastic3DGcmMatrices::M> s2ForwardZ = matrix.A(2).U1.getColumn(indexForwardS2);
	ASSERT_LT(s2ForwardZ(1) /* Vy */ * s2ForwardZ(7) /* Syz */, 0); // shear wave
	ASSERT_GT(matrix.A(2).L(indexForwardS2, indexForwardS2), 0);

	linal::Vector<IdealElastic3DGcmMatrices::M> s2BackwardZ = matrix.A(2).U1.getColumn(indexBackwardS2);
	ASSERT_GT(s2BackwardZ(1) /* Vy */ * s2BackwardZ(7) /* Syz */, 0); // shear wave
	ASSERT_LT(matrix.A(2).L(indexBackwardS2, indexBackwardS2), 0);
};