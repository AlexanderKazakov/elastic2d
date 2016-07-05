#include <gtest/gtest.h>

#include <lib/rheology/models/models.hpp>
#include <lib/rheology/materials/materials.hpp>


using namespace gcm;


template<typename TGcmMatrix>
inline void testTrace(const TGcmMatrix& m) {
	ASSERT_NEAR(linal::trace(m.A), linal::trace(m.L), EQUALITY_TOLERANCE)
			<< "A:" << m.A << "L:" << m.L;
}


template<typename TGcmMatrix>
inline void testEigenvectors(const TGcmMatrix& m,
		const real eps = EQUALITY_TOLERANCE, const bool relative = false) {
	auto AU1 = m.A * m.U1;
	auto U1L = m.U1 * m.L;
	for (int i = 0; i < AU1.M; i++) {
		for (int j = 0; j < AU1.M; j++) {
			real scale = fabs(AU1(i, j)) + fabs(U1L(i, j));
			real bound = (relative && (scale > eps)) ? eps * scale : eps;
			ASSERT_NEAR(AU1(i, j), U1L(i, j), bound) 
					<< "i = " << i << ", j = " << j 
					<< "A:" << m.A << "U1:" << m.U1 << "L:" << m.L;
		}
	}
}


template<typename TGcmMatrix>
inline void testEigenstrings(const TGcmMatrix& m,
		const real eps = EQUALITY_TOLERANCE, const bool relative = false) {
	auto UA = m.U * m.A;
	auto LU = m.L * m.U;
	for (int i = 0; i < UA.M; i++) {
		for (int j = 0; j < UA.M; j++) {
			real scale = fabs(UA(i, j)) + fabs(LU(i, j));
			real bound = (relative && (scale > eps)) ? eps * scale : eps;
			ASSERT_NEAR(UA(i, j), LU(i, j), bound) 
					<< "i = " << i << ", j = " << j 
					<< "A:" << m.A << "U:" << m.U << "L:" << m.L;
		}
	}
}


template<typename TGcmMatrix>
inline void testInverse(const TGcmMatrix& m, const real eps = EQUALITY_TOLERANCE) {
	auto UU1 = m.U * m.U1;
	for (int i = 0; i < UU1.M; i++) {
		for (int j = 0; j < UU1.M; j++) {
			ASSERT_NEAR(UU1(i, j), (i == j), eps) 
					<< "i = " << i << ", j = " << j 
					<< "\nUU1:" << UU1;
		}
	}
}


template<template<int> class TModel, int Dimensionality>
void testIsotropicRotatedGcmMatrix() {
	typedef TModel<Dimensionality>      Model;
	typedef typename Model::MatrixDD    Basis;
	typedef typename Model::RealD       RealD;
	typedef typename Model::GcmMatrix   GcmMatrix;
	
	Utils::seedRand();
	for (int i = 0; i < 10000; i++) {
		auto material = std::make_shared<const IsotropicMaterial>(
				IsotropicMaterial::generateRandomMaterial());
		const Basis b = linal::randomBasis(Basis());
		const real l = Utils::randomReal(0.1, 10);
		
		GcmMatrix m;
		Model::constructGcmMatrix(m, material, b, l);
		testTrace(m);
		testEigenstrings(m, 1000 * EQUALITY_TOLERANCE, true);
		testEigenvectors(m, 1000 * EQUALITY_TOLERANCE, true);
		testInverse(m);
		
		for (int j = 0; j < Dimensionality; j++) {
			RealD n = RealD::Zeros();
			n(j) = 1; //< axis aligned basis
			Model::constructGcmMatrix(m, material, linal::createLocalBasis(n), l);
			testTrace(m);
			testEigenstrings(m, 1000 * EQUALITY_TOLERANCE, true);
			testEigenvectors(m, 1000 * EQUALITY_TOLERANCE, true);
			testInverse(m);
		}
	}
}


TEST(IsotropicGcmMatrix, RotatedBasis) {
	testIsotropicRotatedGcmMatrix<ElasticModel, 1>();
	testIsotropicRotatedGcmMatrix<ElasticModel, 2>();
	testIsotropicRotatedGcmMatrix<ElasticModel, 3>();
	
	testIsotropicRotatedGcmMatrix<AcousticModel, 1>();
	testIsotropicRotatedGcmMatrix<AcousticModel, 2>();
	testIsotropicRotatedGcmMatrix<AcousticModel, 3>();
}


template<typename Model, typename Material>
void testGcmMatrices() {
	typedef typename Model::GCM_MATRICES GCM_MATRICES;
	
	Utils::seedRand();
	for (int i = 0; i < 1000; i++) {
		auto material = std::make_shared<Material>(Material::generateRandomMaterial());
		auto matrix = std::shared_ptr<GCM_MATRICES>(new GCM_MATRICES());
		Model::constructGcmMatrices(matrix, material);
		
		for (int s = 0; s < Model::DIMENSIONALITY; s++) {
			testTrace((*matrix)(s));
			testEigenstrings((*matrix)(s));
			testEigenvectors((*matrix)(s));
			testInverse((*matrix)(s));
		}
	}
}


TEST(OrthotropicGcmMatrix, NotRotatedBasis) {
	testGcmMatrices<ElasticModel<2>, OrthotropicMaterial>();
	testGcmMatrices<ElasticModel<3>, OrthotropicMaterial>();
}


TEST(GcmMatrices, RotatedOrthotropicMaterial) {
	typedef ElasticModel<3>::GCM_MATRICES GCM_MATRICES;
	
	auto matrix = std::shared_ptr<GCM_MATRICES>(new GCM_MATRICES());
	auto material = std::shared_ptr<OrthotropicMaterial>(
			new OrthotropicMaterial(1.6, 
					{163944.8, 3767.9, 3767.9,
					           8875.6, 2899.1,
					                   8875.6,
					                          4282.6, 4282.6, 4282.6}));
	auto test = [&](const Real3 phi) {
		material->anglesOfRotation = phi;
		ElasticModel<3>::constructGcmMatrices(matrix, material);
		for (int s = 0; s < 3; s++) {
			testTrace((*matrix)(s));
			testEigenstrings((*matrix)(s), 0.01, true);
			testEigenvectors((*matrix)(s), 0.01, true);
			testInverse((*matrix)(s), 0.02);
		}
	};
	
	Utils::seedRand();
	for (int i = 0; i < 1000; i++) {	
		test({Utils::randomReal(-2*M_PI, 2*M_PI), 0, 0});
		test({0, Utils::randomReal(-2*M_PI, 2*M_PI), 0});
		test({0, 0, Utils::randomReal(-2*M_PI, 2*M_PI)});
		// TODO - theese ones are failed:
		// test({Utils::randomReal(-2*M_PI, 2*M_PI), 
		//       Utils::randomReal(-2*M_PI, 2*M_PI), 0});
	}
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


