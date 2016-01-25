#include <gtest/gtest.h>

#include <lib/rheology/gcm_matrices/GcmMatrices.hpp>

using namespace gcm;
using namespace gcm::linal;

template<class TGcmMatrices>
class TestGcmMatrices : public testing::Test {

	static const int NUMBER_ITERATIONS = 1000;

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
		srand((unsigned int)time(0));
		for (int i = 0; i < NUMBER_ITERATIONS; i++) {
			typename TGcmMatrices::Material material = TGcmMatrices::Material::generateRandomMaterial();

			testTraces(TGcmMatrices(material));
			testLeftEigenVectors(TGcmMatrices(material));
			testRightEigenVectors(TGcmMatrices(material));
			testInverseMatrix(TGcmMatrices(material));
			testValuesOrdersConcerted(TGcmMatrices(material));
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
typedef Types<GcmMatrices<2, 1, IsotropicMaterial>, GcmMatrices<5, 2, IsotropicMaterial>, GcmMatrices<9, 3, IsotropicMaterial>, 
        GcmMatrices<9, 3, OrthotropicMaterial>> AllImplementations;

INSTANTIATE_TYPED_TEST_CASE_P(AllGcmMatrices, TestGcmMatrices, AllImplementations);
#endif // GTEST_HAS_TYPED_TEST_P


template<>
void TestGcmMatrices<GcmMatrices<2, 1, IsotropicMaterial>>::testValuesOrdersConcerted(const GcmMatrices<2, 1, IsotropicMaterial> &matrix) {
	// Vx 0, Sxx 0
	
	int indexForward = GcmMatrices<2, 1, IsotropicMaterial>::WAVE_COLUMNS.at(Waves::T::P_FORWARD);
	int indexBackward = GcmMatrices<2, 1, IsotropicMaterial>::WAVE_COLUMNS.at(Waves::T::P_BACKWARD);

	linal::Vector<GcmMatrices<2, 1, IsotropicMaterial>::M> pForward = matrix.A(0).U1.getColumn(indexForward);
	ASSERT_LT(pForward(0) /* Vx */ * pForward(1) /* Sxx */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(0).L(indexForward), 0);

	linal::Vector<GcmMatrices<2, 1, IsotropicMaterial>::M> pBackward = matrix.A(0).U1.getColumn(indexBackward);
	ASSERT_GT(pBackward(0) /* Vx */ * pBackward(1) /* Sxx */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(0).L(indexBackward), 0);
};

template<>
void TestGcmMatrices<GcmMatrices<5, 2, IsotropicMaterial>>::testValuesOrdersConcerted(const GcmMatrices<5, 2, IsotropicMaterial> &matrix) {
	// Vx 0, Vy 1, Sxx 2, Sxy 3, Syy 4

	int indexForwardP = GcmMatrices<5, 2, IsotropicMaterial>::WAVE_COLUMNS.at(Waves::T::P_FORWARD);
	int indexBackwardP = GcmMatrices<5, 2, IsotropicMaterial>::WAVE_COLUMNS.at(Waves::T::P_BACKWARD);
	int indexForwardS = GcmMatrices<5, 2, IsotropicMaterial>::WAVE_COLUMNS.at(Waves::T::S1_FORWARD);
	int indexBackwardS = GcmMatrices<5, 2, IsotropicMaterial>::WAVE_COLUMNS.at(Waves::T::S1_BACKWARD);

	linal::Vector<GcmMatrices<5, 2, IsotropicMaterial>::M> pForwardX = matrix.A(0).U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardX(0) /* Vx */ * pForwardX(2) /* Sxx */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(0).L(indexForwardP), 0);

	linal::Vector<GcmMatrices<5, 2, IsotropicMaterial>::M> pBackwardX = matrix.A(0).U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardX(0) /* Vx */ * pBackwardX(2) /* Sxx */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(0).L(indexBackwardP), 0);


	linal::Vector<GcmMatrices<5, 2, IsotropicMaterial>::M> sForwardX = matrix.A(0).U1.getColumn(indexForwardS);
	ASSERT_LT(sForwardX(1) /* Vy */ * sForwardX(3) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.A(0).L(indexForwardS), 0);

	linal::Vector<GcmMatrices<5, 2, IsotropicMaterial>::M> sBackwardX = matrix.A(0).U1.getColumn(indexBackwardS);
	ASSERT_GT(sBackwardX(1) /* Vy */ * sBackwardX(3) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.A(0).L(indexBackwardS), 0);



	linal::Vector<GcmMatrices<5, 2, IsotropicMaterial>::M> pForwardY = matrix.A(1).U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardY(1) /* Vy */ * pForwardY(4) /* Syy */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(1).L(indexForwardP), 0);

	linal::Vector<GcmMatrices<5, 2, IsotropicMaterial>::M> pBackwardY = matrix.A(1).U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardY(1) /* Vy */ * pBackwardY(4) /* Syy */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(1).L(indexBackwardP), 0);


	linal::Vector<GcmMatrices<5, 2, IsotropicMaterial>::M> sForwardY = matrix.A(1).U1.getColumn(indexForwardS);
	ASSERT_LT(sForwardY(0) /* Vx */ * sForwardY(3) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.A(1).L(indexForwardS), 0);

	linal::Vector<GcmMatrices<5, 2, IsotropicMaterial>::M> sBackwardY = matrix.A(1).U1.getColumn(indexBackwardS);
	ASSERT_GT(sBackwardY(0) /* Vx */ * sBackwardY(3) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.A(1).L(indexBackwardS), 0);
};

template<class TGcmMatrices>
void TestGcmMatrices<TGcmMatrices>::testValuesOrdersConcerted(const TGcmMatrices &matrix) {
	static_assert(TGcmMatrices::DIMENSIONALITY == 3, "instantiate your test separately");
	static_assert(TGcmMatrices::M == 9, "instantiate your test separately");
	// Vx 0, Vy 1, Vz 2, Sxx 3, Sxy 4, Sxz 5, Syy 6, Syz 7, Szz 8

	int indexForwardP = TGcmMatrices::WAVE_COLUMNS.at(Waves::T::P_FORWARD);
	int indexBackwardP = TGcmMatrices::WAVE_COLUMNS.at(Waves::T::P_BACKWARD);
	int indexForwardS1 = TGcmMatrices::WAVE_COLUMNS.at(Waves::T::S1_FORWARD);
	int indexBackwardS1 = TGcmMatrices::WAVE_COLUMNS.at(Waves::T::S1_BACKWARD);
	int indexForwardS2 = TGcmMatrices::WAVE_COLUMNS.at(Waves::T::S2_FORWARD);
	int indexBackwardS2 = TGcmMatrices::WAVE_COLUMNS.at(Waves::T::S2_BACKWARD);

	linal::Vector<TGcmMatrices::M> pForwardX = matrix.A(0).U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardX(0) /* Vx */ * pForwardX(3) /* Sxx */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(0).L(indexForwardP), 0);

	linal::Vector<TGcmMatrices::M> pBackwardX = matrix.A(0).U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardX(0) /* Vx */ * pBackwardX(3) /* Sxx */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(0).L(indexBackwardP), 0);


	linal::Vector<TGcmMatrices::M> s1ForwardX = matrix.A(0).U1.getColumn(indexForwardS1);
	ASSERT_LT(s1ForwardX(1) /* Vy */ * s1ForwardX(4) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.A(0).L(indexForwardS1), 0);

	linal::Vector<TGcmMatrices::M> s1BackwardX = matrix.A(0).U1.getColumn(indexBackwardS1);
	ASSERT_GT(s1BackwardX(1) /* Vy */ * s1BackwardX(4) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.A(0).L(indexBackwardS1), 0);


	linal::Vector<TGcmMatrices::M> s2ForwardX = matrix.A(0).U1.getColumn(indexForwardS2);
	ASSERT_LT(s2ForwardX(2) /* Vz */ * s2ForwardX(5) /* Sxz */, 0); // shear wave
	ASSERT_GT(matrix.A(0).L(indexForwardS2), 0);

	linal::Vector<TGcmMatrices::M> s2BackwardX = matrix.A(0).U1.getColumn(indexBackwardS2);
	ASSERT_GT(s2BackwardX(2) /* Vz */ * s2BackwardX(5) /* Sxz */, 0); // shear wave
	ASSERT_LT(matrix.A(0).L(indexBackwardS2), 0);



	linal::Vector<TGcmMatrices::M> pForwardY = matrix.A(1).U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardY(1) /* Vy */ * pForwardY(6) /* Syy */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(1).L(indexForwardP), 0);

	linal::Vector<TGcmMatrices::M> pBackwardY = matrix.A(1).U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardY(1) /* Vy */ * pBackwardY(6) /* Syy */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(1).L(indexBackwardP), 0);


	linal::Vector<TGcmMatrices::M> s1ForwardY = matrix.A(1).U1.getColumn(indexForwardS1);
	ASSERT_LT(s1ForwardY(0) /* Vx */ * s1ForwardY(4) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.A(1).L(indexForwardS1), 0);

	linal::Vector<TGcmMatrices::M> s1BackwardY = matrix.A(1).U1.getColumn(indexBackwardS1);
	ASSERT_GT(s1BackwardY(0) /* Vx */ * s1BackwardY(4) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.A(1).L(indexBackwardS1), 0);


	linal::Vector<TGcmMatrices::M> s2ForwardY = matrix.A(1).U1.getColumn(indexForwardS2);
	ASSERT_LT(s2ForwardY(2) /* Vz */ * s2ForwardY(7) /* Syz */, 0); // shear wave
	ASSERT_GT(matrix.A(1).L(indexForwardS2), 0);

	linal::Vector<TGcmMatrices::M> s2BackwardY = matrix.A(1).U1.getColumn(indexBackwardS2);
	ASSERT_GT(s2BackwardY(2) /* Vz */ * s2BackwardY(7) /* Syz */, 0); // shear wave
	ASSERT_LT(matrix.A(1).L(indexBackwardS2), 0);



	linal::Vector<TGcmMatrices::M> pForwardZ = matrix.A(2).U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardZ(2) /* Vz */ * pForwardZ(8) /* Szz */, 0); // compression or depression wave
	ASSERT_GT(matrix.A(2).L(indexForwardP), 0);

	linal::Vector<TGcmMatrices::M> pBackwardZ = matrix.A(2).U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardZ(2) /* Vz */ * pBackwardZ(8) /* Szz */, 0); // compression or depression wave
	ASSERT_LT(matrix.A(2).L(indexBackwardP), 0);


	linal::Vector<TGcmMatrices::M> s1ForwardZ = matrix.A(2).U1.getColumn(indexForwardS1);
	ASSERT_LT(s1ForwardZ(0) /* Vx */ * s1ForwardZ(5) /* Sxz */, 0); // shear wave
	ASSERT_GT(matrix.A(2).L(indexForwardS1), 0);

	linal::Vector<TGcmMatrices::M> s1BackwardZ = matrix.A(2).U1.getColumn(indexBackwardS1);
	ASSERT_GT(s1BackwardZ(0) /* Vx */ * s1BackwardZ(5) /* Sxz */, 0); // shear wave
	ASSERT_LT(matrix.A(2).L(indexBackwardS1), 0);


	linal::Vector<TGcmMatrices::M> s2ForwardZ = matrix.A(2).U1.getColumn(indexForwardS2);
	ASSERT_LT(s2ForwardZ(1) /* Vy */ * s2ForwardZ(7) /* Syz */, 0); // shear wave
	ASSERT_GT(matrix.A(2).L(indexForwardS2), 0);

	linal::Vector<TGcmMatrices::M> s2BackwardZ = matrix.A(2).U1.getColumn(indexBackwardS2);
	ASSERT_GT(s2BackwardZ(1) /* Vy */ * s2BackwardZ(7) /* Syz */, 0); // shear wave
	ASSERT_LT(matrix.A(2).L(indexBackwardS2), 0);
};