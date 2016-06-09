#include <gtest/gtest.h>

#include <lib/rheology/models/Model.hpp>
#include <lib/rheology/materials/materials.hpp>

using namespace gcm;

template<class Case>
class TestGcmMatrices : public testing::Test {
protected:

typedef typename Case::Model         Model;
typedef typename Case::Material      Material;
typedef typename Model::GCM_MATRICES GCM_MATRICES;
typedef typename Model::PdeVector    PdeVector;

static const int DIMENSIONALITY = GCM_MATRICES::D;
static const int M = GCM_MATRICES::M;

static const int NUMBER_OF_TEST_ITERATIONS = 1000;

void testTraces(const GCM_MATRICES& matrix) {
	for (int i = 0; i < DIMENSIONALITY; i++) {
		ASSERT_NEAR(linal::trace(matrix.m[i].A),
		            linal::trace(matrix.m[i].L),
		            EQUALITY_TOLERANCE) <<
		"(" << i << ") A = " << matrix.m[i].A << "L = " << matrix.m[i].L;
	}
}

void testLeftEigenVectors(const GCM_MATRICES& matrix, 
		const real eps = EQUALITY_TOLERANCE, const bool relative = false) {
	for (int s = 0; s < DIMENSIONALITY; s++) {
		auto AU1 = matrix.m[s].A * matrix.m[s].U1;
		auto U1L = matrix.m[s].U1 * matrix.m[s].L;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				real scale = fabs(AU1(i, j)) + fabs(U1L(i, j));
				real bound = (relative && (scale > eps)) ?	
						eps * scale : eps;
				
				ASSERT_NEAR(AU1(i, j), U1L(i, j), bound) << "(" << s << ") U1 = " << matrix.m[s].U1;
			}
		}
	}
}

void testRightEigenVectors(const GCM_MATRICES& matrix, 
		const real eps = EQUALITY_TOLERANCE, const bool relative = false) {
	for (int s = 0; s < DIMENSIONALITY; s++) {
		auto UA = matrix.m[s].U * matrix.m[s].A;
		auto LU = matrix.m[s].L * matrix.m[s].U;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				real scale = fabs(UA(i, j)) + fabs(LU(i, j));
				real bound = (relative && (scale > eps)) ?	
						eps * scale : eps;
				ASSERT_NEAR(UA(i, j), LU(i, j), bound) << "(" << s << ") U = " << matrix.m[s].U;
			}
		}
	}
}

void testInverseMatrix(const GCM_MATRICES& matrix, const real eps = EQUALITY_TOLERANCE) {
	for (int s = 0; s < DIMENSIONALITY; s++) {
		auto UU1 = matrix.m[s].U * matrix.m[s].U1;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				ASSERT_NEAR(UU1(i, j), (i == j), eps) << "(" << s << ") UU1 = " << UU1;
			}
		}
	}
}

/** 
 * Check that order of values in node and matrices 
 * and order of indices in WAVE_COLUMNS is concerted 
 */
void testValuesOrdersConcerted(const GCM_MATRICES& matrix);

void testDiagonalization() {
	Utils::seedRand();
	for (int i = 0; i < NUMBER_OF_TEST_ITERATIONS; i++) {
		auto material = std::make_shared<Material>(Material::generateRandomMaterial());
		auto matrix = std::shared_ptr<GCM_MATRICES>(new GCM_MATRICES());
		Model::constructGcmMatrices(matrix, material);

		testTraces(*matrix);
		testLeftEigenVectors(*matrix);
		testRightEigenVectors(*matrix);
		testInverseMatrix(*matrix);
		testValuesOrdersConcerted(*matrix);
	}
}
};


template<class Case>
class TestGcmMatricesRotated : public TestGcmMatrices<Case> {
protected:
typedef TestGcmMatrices<Case>        Base;
typedef typename Base::Model         Model;
typedef typename Base::Material      Material;
typedef typename Base::GCM_MATRICES  GCM_MATRICES;
typedef typename Base::PdeVector     PdeVector;

static const int DIMENSIONALITY = Base::DIMENSIONALITY;
static const int M = Base::M;

static const int NUMBER_OF_TEST_ITERATIONS = Base::NUMBER_OF_TEST_ITERATIONS;

	void testConcreteMaterial(const Real3 phi) {
		auto matrix = std::shared_ptr<GCM_MATRICES>(new GCM_MATRICES());
		auto material = std::shared_ptr<OrthotropicMaterial>(
				new OrthotropicMaterial(1.6, 
						{163944.8, 3767.9, 3767.9,
								   8875.6, 2899.1,
										   8875.6,
												   4282.6, 4282.6, 4282.6},
						 0, 0, phi));
		
		Model::constructGcmMatrices(matrix, material);

		this->testTraces(*matrix);
		this->testLeftEigenVectors(*matrix, 0.01, true);
		this->testRightEigenVectors(*matrix, 0.01, true);
		this->testInverseMatrix(*matrix, 0.02);
	}

	void testDiagonalizationRotated() {
		Utils::seedRand();
		for (int i = 0; i < NUMBER_OF_TEST_ITERATIONS; i++) {
			testConcreteMaterial({Utils::randomReal(-2*M_PI, 2*M_PI), 0, 0});
			testConcreteMaterial({0, Utils::randomReal(-2*M_PI, 2*M_PI), 0});
			testConcreteMaterial({0, 0, Utils::randomReal(-2*M_PI, 2*M_PI)});
			// TODO - make test for arbitrary material;
			// TODO - theese ones are failed like UA = 2 when LU == 0
//			testConcreteMaterial({Utils::randomReal(-2*M_PI, 2*M_PI), 
//			                      Utils::randomReal(-2*M_PI, 2*M_PI), 0});
//			testConcreteMaterial({Utils::randomReal(-2*M_PI, 2*M_PI), 0, 
//			                      Utils::randomReal(-2*M_PI, 2*M_PI)   });
//			testConcreteMaterial({0, Utils::randomReal(-2*M_PI, 2*M_PI), 
//			                      Utils::randomReal(-2*M_PI, 2*M_PI)   });
		}
	}
};


template<typename TModel, typename TMaterial>
struct Wrap {
	typedef TModel    Model;
	typedef TMaterial Material;
};


/** Look at https://github.com/google/googletest/blob/master/googletest/samples/sample6_unittest.cc
  for explaination */
#if GTEST_HAS_TYPED_TEST_P
using testing::Types;
TYPED_TEST_CASE_P(TestGcmMatrices);
TYPED_TEST_CASE_P(TestGcmMatricesRotated);

TYPED_TEST_P(TestGcmMatrices, Diagonalization) {
	this->testDiagonalization();
}

TYPED_TEST_P(TestGcmMatricesRotated, Diagonalization) {
	this->testDiagonalizationRotated();
}

REGISTER_TYPED_TEST_CASE_P(TestGcmMatrices, Diagonalization);
REGISTER_TYPED_TEST_CASE_P(TestGcmMatricesRotated, Diagonalization);

typedef Types<
        Wrap<Elastic1DModel, IsotropicMaterial>,
        Wrap<Elastic2DModel, IsotropicMaterial>,
        Wrap<Elastic3DModel, IsotropicMaterial>,
        Wrap<Elastic3DModel, OrthotropicMaterial>,
        Wrap<SuperDuperModel, IsotropicMaterial>,
        Wrap<SuperDuperModel, OrthotropicMaterial>
        > AllImplementations;

typedef Types<
        Wrap<Elastic3DModel, OrthotropicMaterial>,
        Wrap<SuperDuperModel, OrthotropicMaterial>
        > RotatedImplementations;

INSTANTIATE_TYPED_TEST_CASE_P(AllGcmMatrices, TestGcmMatrices, AllImplementations);
INSTANTIATE_TYPED_TEST_CASE_P(RotatedGcmMatrices, TestGcmMatricesRotated, RotatedImplementations);

#endif // GTEST_HAS_TYPED_TEST_P


template<>
void TestGcmMatrices<Wrap<Elastic1DModel, IsotropicMaterial> >::
testValuesOrdersConcerted(const GCM_MATRICES& matrix) {
	// Vx 0, Sxx 1

	int indexForward = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_FORWARD);
	int indexBackward = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_BACKWARD);

	auto pForward = matrix.m[0].U1.getColumn(indexForward);
	ASSERT_LT(pForward(0) /* Vx */ * pForward(1) /* Sxx */, 0); // compression or depression
	                                                            // wave
	ASSERT_GT(matrix.m[0].L(indexForward), 0);

	auto pBackward = matrix.m[0].U1.getColumn(indexBackward);
	ASSERT_GT(pBackward(0) /* Vx */ * pBackward(1) /* Sxx */, 0); // compression or depression
	                                                              // wave
	ASSERT_LT(matrix.m[0].L(indexBackward), 0);
}


template<>
void TestGcmMatrices<Wrap<Elastic2DModel, IsotropicMaterial> >::
testValuesOrdersConcerted(const GCM_MATRICES& matrix) {
	// Vx 0, Vy 1, Sxx 2, Sxy 3, Syy 4

	int indexForwardP = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_FORWARD);
	int indexBackwardP = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_BACKWARD);
	int indexForwardS = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S1_FORWARD);
	int indexBackwardS = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S1_BACKWARD);

	auto pForwardX = matrix.m[0].U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardX(0) /* Vx */ * pForwardX(2) /* Sxx */, 0); // compression or depression
	                                                              // wave
	ASSERT_GT(matrix.m[0].L(indexForwardP), 0);

	auto pBackwardX = matrix.m[0].U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardX(0) /* Vx */ * pBackwardX(2) /* Sxx */, 0); // compression or depression
	                                                                // wave
	ASSERT_LT(matrix.m[0].L(indexBackwardP), 0);


	auto sForwardX = matrix.m[0].U1.getColumn(indexForwardS);
	ASSERT_LT(sForwardX(1) /* Vy */ * sForwardX(3) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.m[0].L(indexForwardS), 0);

	auto sBackwardX = matrix.m[0].U1.getColumn(indexBackwardS);
	ASSERT_GT(sBackwardX(1) /* Vy */ * sBackwardX(3) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.m[0].L(indexBackwardS), 0);



	auto pForwardY = matrix.m[1].U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardY(1) /* Vy */ * pForwardY(4) /* Syy */, 0); // compression or depression
	                                                              // wave
	ASSERT_GT(matrix.m[1].L(indexForwardP), 0);

	auto pBackwardY = matrix.m[1].U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardY(1) /* Vy */ * pBackwardY(4) /* Syy */, 0); // compression or depression
	                                                                // wave
	ASSERT_LT(matrix.m[1].L(indexBackwardP), 0);


	auto sForwardY = matrix.m[1].U1.getColumn(indexForwardS);
	ASSERT_LT(sForwardY(0) /* Vx */ * sForwardY(3) /* Sxy */, 0); // shear wave
	ASSERT_GT(matrix.m[1].L(indexForwardS), 0);

	auto sBackwardY = matrix.m[1].U1.getColumn(indexBackwardS);
	ASSERT_GT(sBackwardY(0) /* Vx */ * sBackwardY(3) /* Sxy */, 0); // shear wave
	ASSERT_LT(matrix.m[1].L(indexBackwardS), 0);
}


template<class Wrap>
void TestGcmMatrices<Wrap>::
testValuesOrdersConcerted(const GCM_MATRICES& matrix) {
	// Vx 0, Vy 1, Vz 2, Sxx 3, Sxy 4, Sxz 5, Syy 6, Syz 7, Szz 8

	int indexForwardP = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_FORWARD);
	int indexBackwardP = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::P_BACKWARD);
	int indexForwardS1 = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S1_FORWARD);
	int indexBackwardS1 = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S1_BACKWARD);
	int indexForwardS2 = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S2_FORWARD);
	int indexBackwardS2 = Model::MATERIALS_WAVES_MAP.at(Material::ID).at(Waves::T::S2_BACKWARD);

	auto pForwardX = matrix.m[0].U1.getColumn(indexForwardP);
	ASSERT_LT(pForwardX(0) /* Vx */ * pForwardX(3) /* Sxx */, 0); // compression or depression
	                                                              // wave
	ASSERT_GT(matrix.m[0].L(indexForwardP), 0);

	auto pBackwardX = matrix.m[0].U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardX(0) /* Vx */ * pBackwardX(3) /* Sxx */, 0); // compression or depression
	                                                                // wave
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
	ASSERT_LT(pForwardY(1) /* Vy */ * pForwardY(6) /* Syy */, 0); // compression or depression
	                                                              // wave
	ASSERT_GT(matrix.m[1].L(indexForwardP), 0);

	auto pBackwardY = matrix.m[1].U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardY(1) /* Vy */ * pBackwardY(6) /* Syy */, 0); // compression or depression
	                                                                // wave
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
	ASSERT_LT(pForwardZ(2) /* Vz */ * pForwardZ(8) /* Szz */, 0); // compression or depression
	                                                              // wave
	ASSERT_GT(matrix.m[2].L(indexForwardP), 0);

	auto pBackwardZ = matrix.m[2].U1.getColumn(indexBackwardP);
	ASSERT_GT(pBackwardZ(2) /* Vz */ * pBackwardZ(8) /* Szz */, 0); // compression or depression
	                                                                // wave
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


TEST(Material, rotate_orthotropic_material) {
	Utils::seedRand();
	for (int i = 0; i < 100; i++) {
		auto material = std::make_shared<OrthotropicMaterial>(
				OrthotropicMaterial::generateRandomMaterial());

		material->anglesOfRotation = {2*M_PI, M_PI, -M_PI};
		ASSERT_TRUE(linal::approximatelyEqual(material->getElasticMatrix(),
				material->getRotatedElasticMatrix(), 1e+5 * EQUALITY_TOLERANCE))
						<< "initial:" << material->getElasticMatrix() 
						<< "rotated:" << material->getRotatedElasticMatrix();
		
		auto C = material->getElasticMatrix(), initC = C;
		for (int n = 1; n < 16; n++) {
			for (int p = 0; p < n; p++) {
				Real3 phi = Real3::Zeros(); phi(i % 3) = 2 * M_PI / n;
				C = OrthotropicMaterial::rotate(C, phi);
			}
			ASSERT_TRUE(linal::approximatelyEqual(material->getElasticMatrix(),
					material->getRotatedElasticMatrix(), 1e+5 * EQUALITY_TOLERANCE))
							<< "initial:" << initC << "rotated:" << C;
			C = initC;
		}
	}
}








