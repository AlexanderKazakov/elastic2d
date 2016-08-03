#include <gtest/gtest.h>

#include <lib/Engine.hpp>
#include <lib/util/Area.hpp>
#include <lib/rheology/models/models.hpp>

#include <lib/mesh/DefaultMesh.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>


using namespace gcm;


struct Wrapper {
	typedef DefaultMesh<ElasticModel<2>, CubicGrid<2>, IsotropicMaterial> Mesh;
	static const Mesh* getMesh(const Engine& engine) {
		const AbstractGrid* grid = engine.getSolver()->getActualMesh();
		const Mesh* mesh = dynamic_cast<const Mesh*>(grid);
		assert_true(mesh);
		return mesh;
	}
};



TEST(Engine, runStatement) {
	Task task;
	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::CUBIC;
	
	task.bodies = {
		{Materials::T::ISOTROPIC, Models::T::ELASTIC, 0},
	};
	
	Statement statement;
	task.cubicGrid.borderSize = 5;
	statement.globalSettings.CourantNumber = 4.5;
	statement.materialConditions.byAreas.defaultMaterial = 
			std::make_shared<IsotropicMaterial>(4, 2, 0.5);
	task.cubicGrid.sizes = {20, 40};
	task.cubicGrid.lengths = {7, 3};
	statement.globalSettings.numberOfSnaps = 9;
	statement.globalSettings.requiredTime = 100.0;
	
	Statement::InitialCondition::Wave wave;
	wave.waveType = Waves::T::S1_FORWARD;
	wave.direction = 1; // along y
	wave.quantity = PhysicalQuantities::T::Vx;
	wave.quantityValue = 1;
	Real3 min({ -1, 0.1125, -1});
	Real3 max({ 8, 0.6375, 1});
	wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
	statement.initialCondition.waves.push_back(wave);
	
	task.statements.push_back(statement);
	
	Engine engine(task);
	engine.beforeStatement(statement);
	
	// s-wave
	auto expected = Wrapper::getMesh(engine)->pde({task.cubicGrid.sizes.at(0) / 2, 3});
	engine.runStatement();
	auto actual = Wrapper::getMesh(engine)->pde({task.cubicGrid.sizes.at(0) / 2, 22});
	
	ASSERT_TRUE(linal::approximatelyEqual(expected, actual)) << expected << actual;
}


TEST(Engine, TwoLayersDifferentRho) {
	real rho2rho0Initial = 0.25;
	for (int i = 0; i < 5; i++) {
		
		Task task;
		task.globalSettings.dimensionality = 2;
		task.globalSettings.gridId = Grids::T::CUBIC;
		
		task.bodies = {
			{Materials::T::ISOTROPIC, Models::T::ELASTIC, 0},
		};
		
		Statement statement;
		task.cubicGrid.borderSize = 3;
		statement.globalSettings.CourantNumber = 1.5;
		
		real rho0 = 1, lambda0 = 2, mu0 = 0.8;
		real rho2rho0 = rho2rho0Initial * pow(2, i), lambda2lambda0 = 1, mu2mu0 = 1;
		real rho = rho2rho0 * rho0, lambda = lambda2lambda0 * lambda0, mu = mu2mu0 * mu0;
		
		statement.materialConditions.byAreas.defaultMaterial = 
				std::make_shared<IsotropicMaterial>(rho0, lambda0, mu0);
		Statement::MaterialCondition::ByAreas::Inhomogenity newMaterial;
		newMaterial.area = std::make_shared<AxisAlignedBoxArea>(
				Real3({-10, 0.5 - 1e-5,-10}), Real3({10, 10, 10}));
		newMaterial.material = std::make_shared<IsotropicMaterial>(rho, lambda, mu);
		statement.materialConditions.byAreas.materials.push_back(newMaterial);
		
		task.cubicGrid.sizes = {50, 100};
		task.cubicGrid.lengths = {2, 1};
		statement.globalSettings.numberOfSnaps = 0;
		statement.globalSettings.requiredTime = 0.24;
		
		Statement::InitialCondition::Wave wave;
		wave.waveType = Waves::T::P_FORWARD;
		wave.direction = 1; // along y
		wave.quantity = PhysicalQuantities::T::Vy;
		wave.quantityValue = -2;
		Real3 min({ -1, 0.015, -1});
		Real3 max({ 4, 0.455, 1});
		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
		statement.initialCondition.waves.push_back(wave);
		
		task.statements.push_back(statement);
		
		Engine engine(task);
		engine.beforeStatement(statement);
		
		int leftNodeIndex = (int) (task.cubicGrid.sizes.at(1) * 0.25);
		auto init = Wrapper::getMesh(engine)->pdeVars(
				{task.cubicGrid.sizes.at(0) / 2, leftNodeIndex});
		
		engine.runStatement();
		
//		int rightNodeIndex = (int) (task.cubicGrid.sizes(1) * 0.7);
		auto reflect = Wrapper::getMesh(engine)->pdeVars(
				{task.cubicGrid.sizes.at(0) / 2, leftNodeIndex});
				
//		auto transfer = Wrapper::getMesh(engine)->pdeVars(
//				{task.cubicGrid.sizes(0) / 2, rightNodeIndex});
		
		real E0 = mu0 * (3 * lambda0 + 2 * mu0) / (lambda0 + mu0); // Young's modulus
		real Z0 = sqrt(E0 * rho0);                                 // acoustic impedance
		real E = mu * (3 * lambda + 2 * mu) / (lambda + mu);       // Young's modulus
		real Z = sqrt(E * rho);                                    // acoustic impedance
		
		ASSERT_NEAR(reflect.sigma(1, 1) / init.sigma(1, 1),
		            (Z - Z0) / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(reflect.velocity(1) / init.velocity(1),
		            (Z0 - Z) / (Z + Z0),
		            1e-2);
//		TODO - complete stable beautiful tests
//		ASSERT_NEAR(transfer.sigma(1, 1) / init.sigma(1, 1),
//		            2 * Z / (Z + Z0),
//		            1e-2);
//		ASSERT_NEAR(transfer.velocity(1) / init.velocity(1),
//		            2 * Z0 / (Z + Z0),
//		            1e-2);
		
		
	}
}


TEST(Engine, TwoLayersDifferentE) {
	real E2E0Initial = 0.25;
	for (int i = 0; i < 5; i++) {
		
		Task task;
		task.globalSettings.dimensionality = 2;
		task.globalSettings.gridId = Grids::T::CUBIC;
		
		task.bodies = {
			{Materials::T::ISOTROPIC, Models::T::ELASTIC, 0},
		};
		
		Statement statement;
		task.cubicGrid.borderSize = 3;
		statement.globalSettings.CourantNumber = 1.5;
		
		real rho0 = 1, lambda0 = 2, mu0 = 0.8;
		real rho2rho0 = 1, lambda2lambda0 = E2E0Initial * pow(2, i), mu2mu0 = E2E0Initial * pow(2, i);
		real rho = rho2rho0 * rho0, lambda = lambda2lambda0 * lambda0, mu = mu2mu0 * mu0;
		
		statement.materialConditions.byAreas.defaultMaterial = 
				std::make_shared<IsotropicMaterial>(rho0, lambda0, mu0);
		Statement::MaterialCondition::ByAreas::Inhomogenity newMaterial;
		newMaterial.area = std::make_shared<AxisAlignedBoxArea>(
				Real3({-10, 0.5 - 1e-5, -10}), Real3({10, 10, 10}));
		newMaterial.material = std::make_shared<IsotropicMaterial>(rho, lambda, mu);
		statement.materialConditions.byAreas.materials.push_back(newMaterial);
		
		task.cubicGrid.sizes = {50, 100};
		task.cubicGrid.lengths = {2, 1};
		
		statement.globalSettings.numberOfSnaps = 0;
		statement.globalSettings.requiredTime = 0.24;
		
		Statement::InitialCondition::Wave wave;
		wave.waveType = Waves::T::P_FORWARD;
		wave.direction = 1; // along y
		wave.quantity = PhysicalQuantities::T::Vy;
		wave.quantityValue = -2;
		Real3 min({ -1, 0.015, -1});
		Real3 max({ 4, 0.455, 1});
		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
		statement.initialCondition.waves.push_back(wave);
		
		task.statements.push_back(statement);
		
		Engine engine(task);
		engine.beforeStatement(statement);

		int leftNodeIndex = (int) (task.cubicGrid.sizes.at(1) * 0.25);
		auto init = Wrapper::getMesh(engine)->pdeVars(
				{task.cubicGrid.sizes.at(0) / 2, leftNodeIndex});

		engine.runStatement();

//		int rightNodeIndex = (int) (task.cubicGrid.sizes(1) * 0.7);
		auto reflect = Wrapper::getMesh(engine)->pdeVars(
				{task.cubicGrid.sizes.at(0) / 2, leftNodeIndex});
//		auto transfer = Wrapper::getMesh(engine)->pdeVars(
//		{task.cubicGrid.sizes(0) / 2, rightNodeIndex});

		real E0 = mu0 * (3 * lambda0 + 2 * mu0) / (lambda0 + mu0); // Young's modulus
		real Z0 = sqrt(E0 * rho0);                                 // acoustic impedance
		real E = mu * (3 * lambda + 2 * mu) / (lambda + mu);       // Young's modulus
		real Z = sqrt(E * rho);                                    // acoustic impedance

		ASSERT_NEAR(reflect.sigma(1, 1) / init.sigma(1, 1),
		            (Z - Z0) / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(reflect.velocity(1) / init.velocity(1),
		            (Z0 - Z) / (Z + Z0),
		            1e-2);

//		ASSERT_NEAR(transfer.sigma(1, 1) / init.sigma(1, 1),
//		            2 * Z / (Z + Z0),
//		            1e-2);
//		ASSERT_NEAR(transfer.velocity(1) / init.velocity(1),
//		            2 * Z0 / (Z + Z0),
//		            1e-2);


	}
}
