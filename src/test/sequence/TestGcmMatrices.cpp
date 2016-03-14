#include <gtest/gtest.h>

#include <lib/rheology/models/Model.hpp>
#include <lib/rheology/materials/materials.hpp>

using namespace gcm;

template<class Case>
class TestGcmMatrices : public testing::Test {

	typedef typename Case::Model           Model;
	typedef typename Case::Material        Material;
	typedef typename Model::GCM_MATRICES   GCM_MATRICES;
	typedef typename Model::PdeVector      PdeVector;

	static const int DIMENSIONALITY = GCM_MATRICES::DIMENSIONALITY;
	static const int M = GCM_MATRICES::M;

	static const int NUMBER_OF_TEST_ITERATIONS = 1000;

	void testTraces(const GCM_MATRICES& matrix) {
		for (int i = 0; i < DIMENSIONALITY; i++) {
			ASSERT_NEAR(matrix.m[i].A.trace(), matrix.m[i].L.trace(), EQUALITY_TOLERANCE)
				<< "(" << i << ") A = " << matrix.m[i].A << "L = " << matrix.m[i].L;
		}
	}

	void testLeftEigenVectors(const GCM_MATRICES& matrix) {
		for (int s = 0; s < DIMENSIONALITY; s++) {
			auto AU1 = matrix.m[s].A * matrix.m[s].U1;
			auto U1L = matrix.m[s].U1 * matrix.m[s].L;
			for (int i = 0; i < M; i++) {
				for (int j = 0; j < M; j++) {
					ASSERT_NEAR(AU1(i, j), U1L(i, j), EQUALITY_TOLERANCE) << "(" << s << ") U1 = " << matrix.m[s].U1;
				}
			}
		}
	}

	void testRightEigenVectors(const GCM_MATRICES& matrix) {
		for (int s = 0; s < DIMENSIONALITY; s++) {
			auto UA = matrix.m[s].U * matrix.m[s].A;
			auto LU = matrix.m[s].L * matrix.m[s].U;
			for(int i = 0; i < M; i++) {
				for(int j = 0; j < M; j++) {
					ASSERT_NEAR(UA(i, j), LU(i, j), EQUALITY_TOLERANCE) << "(" << s << ") U = " << matrix.m[s].U;
				}
			}
		}
	}

	void testInverseMatrix(const GCM_MATRICES& matrix) {
		for (int s = 0; s < DIMENSIONALITY; s++) {
			auto UU1 = matrix.m[s].U * matrix.m[s].U1;
			for(int i = 0; i < M; i++) {
				for(int j = 0; j < M; j++) {
					ASSERT_NEAR(UU1(i, j), (i == j), EQUALITY_TOLERANCE) << "(" << s << ") UU1 = " << UU1;
				}
			}
		}
	}

	/** Check that order of values in node and matrices and order of indices in WAVE_COLUMNS is concerted */
	void testValuesOrdersConcerted(const GCM_MATRICES& matrix);

protected:
	void testDiagonalization() {
		srand((unsigned int)time(0));
		for (int i = 0; i < NUMBER_OF_TEST_ITERATIONS; i++) {
			Material material = Material::generateRandomMaterial();
			auto matrix = std::shared_ptr<GCM_MATRICES>(new GCM_MATRICES());
			Model::constructGcmMatrices(matrix, PdeVector::zeros(), material);

			testTraces(*matrix);
			testLeftEigenVectors(*matrix);
			testRightEigenVectors(*matrix);
			testInverseMatrix(*matrix);
			testValuesOrdersConcerted(*matrix);
		}
	}
};

template <typename TModel, typename TMaterial>
struct Wrap {
	typedef TModel     Model;
	typedef TMaterial  Material;
};


/** Look at https://github.com/google/googletest/blob/master/googletest/samples/sample6_unittest.cc for explaination */
#if GTEST_HAS_TYPED_TEST_P
using testing::Types;
TYPED_TEST_CASE_P(TestGcmMatrices);

TYPED_TEST_P(TestGcmMatrices, Diagonalization) {
	this->testDiagonalization();
}

REGISTER_TYPED_TEST_CASE_P(TestGcmMatrices, Diagonalization);

typedef Types<
		Wrap<Elastic1DModel, IsotropicMaterial>,
		Wrap<Elastic2DModel, IsotropicMaterial>,
		Wrap<Elastic3DModel, IsotropicMaterial>,
		Wrap<Elastic3DModel, OrthotropicMaterial>,
		Wrap<SuperDuperModel, IsotropicMaterial>,
		Wrap<SuperDuperModel, OrthotropicMaterial>
> AllImplementations;

INSTANTIATE_TYPED_TEST_CASE_P(AllGcmMatrices, TestGcmMatrices, AllImplementations);
#endif // GTEST_HAS_TYPED_TEST_P


template<>
void TestGcmMatrices<Wrap<Elastic1DModel, IsotropicMaterial>>::
testValuesOrdersConcerted(const GCM_MATRICES& matrix) {
	// Vx 0, Sxx 1
	
	int indexForward = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_FORWARD);
	int indexBackward = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_BACKWARD);

	auto pForward = matrix.m[0].U1.getColumn(indexForward);
	ASSERT_LT(pForward(0) /* Vx */ * pForward(1) /* Sxx */, 0); // compression or depression wave
	ASSERT_GT(matrix.m[0].L(indexForward), 0);

	auto pBackward = matrix.m[0].U1.getColumn(indexBackward);
	ASSERT_GT(pBackward(0) /* Vx */ * pBackward(1) /* Sxx */, 0); // compression or depression wave
	ASSERT_LT(matrix.m[0].L(indexBackward), 0);
}

template<>
void TestGcmMatrices<Wrap<Elastic2DModel, IsotropicMaterial>>::
testValuesOrdersConcerted(const GCM_MATRICES& matrix) {
	// Vx 0, Vy 1, Sxx 2, Sxy 3, Syy 4

	int indexForwardP = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_FORWARD);
	int indexBackwardP = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_BACKWARD);
	int indexForwardS = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S1_FORWARD);
	int indexBackwardS = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S1_BACKWARD);

	auto pForwardX = matrix.m[0].U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardX(0) /* Vx */ * pForwardX(2) /* Sxx */, 0); // compression or depression wave
	ASSERT_GT(matrix.m[0].L(indexForwardP), 0);

	auto pBackwardX = matrix.m[0].U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardX(0) /* Vx */ * pBackwardX(2) /* Sxx */, 0); // compression or depression wave
	ASSERT_LT(matrix.m[0].L(indexBackwardP), 0);


	auto sForwardX = matrix.m[0].U1.getColumn(indexForwardS);
	ASSERT_LT(sForwardX(1) /* Vy */ * sForwardX(3) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.m[0].L(indexForwardS), 0);

	auto sBackwardX = matrix.m[0].U1.getColumn(indexBackwardS);
	ASSERT_GT(sBackwardX(1) /* Vy */ * sBackwardX(3) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.m[0].L(indexBackwardS), 0);



	auto pForwardY = matrix.m[1].U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardY(1) /* Vy */ * pForwardY(4) /* Syy */, 0); // compression or depression wave
	ASSERT_GT(matrix.m[1].L(indexForwardP), 0);

	auto pBackwardY = matrix.m[1].U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardY(1) /* Vy */ * pBackwardY(4) /* Syy */, 0); // compression or depression wave
	ASSERT_LT(matrix.m[1].L(indexBackwardP), 0);


	auto sForwardY = matrix.m[1].U1.getColumn(indexForwardS);
	ASSERT_LT(sForwardY(0) /* Vx */ * sForwardY(3) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.m[1].L(indexForwardS), 0);

	auto sBackwardY = matrix.m[1].U1.getColumn(indexBackwardS);
	ASSERT_GT(sBackwardY(0) /* Vx */ * sBackwardY(3) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.m[1].L(indexBackwardS), 0);
}

template<class Wrap>
void TestGcmMatrices<Wrap>::testValuesOrdersConcerted(const GCM_MATRICES& matrix) {
	// Vx 0, Vy 1, Vz 2, Sxx 3, Sxy 4, Sxz 5, Syy 6, Syz 7, Szz 8

	int indexForwardP = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_FORWARD);
	int indexBackwardP = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_BACKWARD);
	int indexForwardS1 = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S1_FORWARD);
	int indexBackwardS1 = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S1_BACKWARD);
	int indexForwardS2 = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S2_FORWARD);
	int indexBackwardS2 = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S2_BACKWARD);

	auto pForwardX = matrix.m[0].U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardX(0) /* Vx */ * pForwardX(3) /* Sxx */, 0); // compression or depression wave
	ASSERT_GT(matrix.m[0].L(indexForwardP), 0);

	auto pBackwardX = matrix.m[0].U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardX(0) /* Vx */ * pBackwardX(3) /* Sxx */, 0); // compression or depression wave
	ASSERT_LT(matrix.m[0].L(indexBackwardP), 0);


	auto s1ForwardX = matrix.m[0].U1.getColumn(indexForwardS1);
	ASSERT_LT(s1ForwardX(1) /* Vy */ * s1ForwardX(4) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.m[0].L(indexForwardS1), 0);

	auto s1BackwardX = matrix.m[0].U1.getColumn(indexBackwardS1);
	ASSERT_GT(s1BackwardX(1) /* Vy */ * s1BackwardX(4) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.m[0].L(indexBackwardS1), 0);


	auto s2ForwardX = matrix.m[0].U1.getColumn(indexForwardS2);
	ASSERT_LT(s2ForwardX(2) /* Vz */ * s2ForwardX(5) /* Sxz */, 0); // shear wave
	ASSERT_GT(matrix.m[0].L(indexForwardS2), 0);

	auto s2BackwardX = matrix.m[0].U1.getColumn(indexBackwardS2);
	ASSERT_GT(s2BackwardX(2) /* Vz */ * s2BackwardX(5) /* Sxz */, 0); // shear wave
	ASSERT_LT(matrix.m[0].L(indexBackwardS2), 0);



	auto pForwardY = matrix.m[1].U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardY(1) /* Vy */ * pForwardY(6) /* Syy */, 0); // compression or depression wave
	ASSERT_GT(matrix.m[1].L(indexForwardP), 0);

	auto pBackwardY = matrix.m[1].U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardY(1) /* Vy */ * pBackwardY(6) /* Syy */, 0); // compression or depression wave
	ASSERT_LT(matrix.m[1].L(indexBackwardP), 0);


	auto s1ForwardY = matrix.m[1].U1.getColumn(indexForwardS1);
	ASSERT_LT(s1ForwardY(0) /* Vx */ * s1ForwardY(4) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.m[1].L(indexForwardS1), 0);

	auto s1BackwardY = matrix.m[1].U1.getColumn(indexBackwardS1);
	ASSERT_GT(s1BackwardY(0) /* Vx */ * s1BackwardY(4) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.m[1].L(indexBackwardS1), 0);


	auto s2ForwardY = matrix.m[1].U1.getColumn(indexForwardS2);
	ASSERT_LT(s2ForwardY(2) /* Vz */ * s2ForwardY(7) /* Syz */, 0); // shear wave
	ASSERT_GT(matrix.m[1].L(indexForwardS2), 0);

	auto s2BackwardY = matrix.m[1].U1.getColumn(indexBackwardS2);
	ASSERT_GT(s2BackwardY(2) /* Vz */ * s2BackwardY(7) /* Syz */, 0); // shear wave
	ASSERT_LT(matrix.m[1].L(indexBackwardS2), 0);



	auto pForwardZ = matrix.m[2].U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardZ(2) /* Vz */ * pForwardZ(8) /* Szz */, 0); // compression or depression wave
	ASSERT_GT(matrix.m[2].L(indexForwardP), 0);

	auto pBackwardZ = matrix.m[2].U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardZ(2) /* Vz */ * pBackwardZ(8) /* Szz */, 0); // compression or depression wave
	ASSERT_LT(matrix.m[2].L(indexBackwardP), 0);


	auto s1ForwardZ = matrix.m[2].U1.getColumn(indexForwardS1);
	ASSERT_LT(s1ForwardZ(0) /* Vx */ * s1ForwardZ(5) /* Sxz */, 0); // shear wave
	ASSERT_GT(matrix.m[2].L(indexForwardS1), 0);

	auto s1BackwardZ = matrix.m[2].U1.getColumn(indexBackwardS1);
	ASSERT_GT(s1BackwardZ(0) /* Vx */ * s1BackwardZ(5) /* Sxz */, 0); // shear wave
	ASSERT_LT(matrix.m[2].L(indexBackwardS1), 0);


	auto s2ForwardZ = matrix.m[2].U1.getColumn(indexForwardS2);
	ASSERT_LT(s2ForwardZ(1) /* Vy */ * s2ForwardZ(7) /* Syz */, 0); // shear wave
	ASSERT_GT(matrix.m[2].L(indexForwardS2), 0);

	auto s2BackwardZ = matrix.m[2].U1.getColumn(indexBackwardS2);
	ASSERT_GT(s2BackwardZ(1) /* Vy */ * s2BackwardZ(7) /* Syz */, 0); // shear wave
	ASSERT_LT(matrix.m[2].L(indexBackwardS2), 0);
}
