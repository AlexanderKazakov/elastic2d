#include <gtest/gtest.h>

#include <lib/util/areas/areas.hpp>

#include <test/wrappers/Wrappers.hpp>
#include <lib/rheology/models/Model.hpp>

#include <lib/util/snapshot/VtkSnapshotter.hpp>

using namespace gcm;

TEST(Engine, runStatementForTest)
{
	Task task;
	Statement statement;
	task.dimensionality = 2;
	task.borderSize = 5;
	statement.CourantNumber = 4.5;
	statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(4, 2, 0.5);
	task.sizes(0) = 20;
	task.sizes(1) = 40;
	task.lengthes = {7, 3, 1};
	statement.numberOfSnaps = 9;
	statement.T = 100.0;

	Statement::InitialCondition::Wave wave;
	wave.waveType = Waves::T::S1_FORWARD;
	wave.direction = 1; // along y
	wave.quantity = PhysicalQuantities::T::Vx;
	wave.quantityValue = 1;
	linal::Vector3 min({ -1, 0.1125, -1});
	linal::Vector3 max({ 8, 0.6375, 1});
	wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
	statement.initialCondition.waves.push_back(wave);
	
	task.statements.push_back(statement);

	EngineWrapper<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>> engine;
	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>());
	engine.initialize(task);
	engine.beforeStatementForTest(statement);
	auto sWave = engine.getSolverForTest()->getMesh()->pde({task.sizes(0) / 2, 3, 0});
	engine.runStatementForTest();
	ASSERT_EQ(sWave, engine.getSolverForTest()->getMesh()->pde({task.sizes(0) / 2, 22, 0}));
}


TEST(Engine, TwoLayersDifferentRho)
{
	real rho2rho0Initial = 0.25;
	for (int i = 0; i < 5; i++) {

		Task task;
		Statement statement;
		task.borderSize = 3;
		task.dimensionality = 2;
		statement.CourantNumber = 1.5;

		real rho0 = 1, lambda0 = 2, mu0 = 0.8;
		real rho2rho0 = rho2rho0Initial * pow(2, i), lambda2lambda0 = 1, mu2mu0 = 1;
		real rho = rho2rho0 * rho0, lambda = lambda2lambda0 * lambda0, mu = mu2mu0 * mu0;

		statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(rho0, lambda0, mu0);
		Statement::MaterialCondition::Inhomogenity newMaterial;
		newMaterial.area = std::make_shared<AxisAlignedBoxArea>(linal::Vector3({-10, 0.5 - 1e-5, -10}), linal::Vector3({10, 10, 10}));
		newMaterial.material = std::make_shared<IsotropicMaterial>(rho, lambda, mu);
		statement.materialConditions.materials.push_back(newMaterial);

		task.sizes(0) = 50;
		task.sizes(1) = 100;
		task.lengthes = {2, 1, 1};
		statement.numberOfSnaps = 0;
		statement.T = 0.24;

		Statement::InitialCondition::Wave wave;
		wave.waveType = Waves::T::P_FORWARD;
		wave.direction = 1; // along y
		wave.quantity = PhysicalQuantities::T::Vy;
		wave.quantityValue = -2;
		linal::Vector3 min({ -1, 0.015, -1});
		linal::Vector3 max({ 4, 0.455, 1});
		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
		statement.initialCondition.waves.push_back(wave);

		task.statements.push_back(statement);

		EngineWrapper<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>> engine;
		engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>());
		engine.initialize(task);
		engine.beforeStatementForTest(statement);

		int leftNodeIndex = (int) (task.sizes(1) * 0.25);
		auto init = engine.getSolverForTest()->getMesh()->pde({task.sizes(0) / 2, leftNodeIndex, 0});

		engine.runStatementForTest();

//		int rightNodeIndex = (int) (task.sizes(1) * 0.7);
		auto reflect = engine.getSolverForTest()->getMesh()->pde({task.sizes(0) / 2, leftNodeIndex, 0});
//		auto transfer = engine.getSolverForTest()->getMesh()->pde({task.sizes(0) / 2, rightNodeIndex, 0});

		real E0 = mu0 * (3 * lambda0 + 2 * mu0) / (lambda0 + mu0); // Young's modulus
		real Z0 = sqrt(E0 * rho0); // acoustic impedance
		real E = mu * (3 * lambda + 2 * mu) / (lambda + mu); // Young's modulus
		real Z = sqrt(E * rho); // acoustic impedance

		ASSERT_NEAR(reflect.sigma(1, 1) / init.sigma(1, 1),
		            (Z - Z0) / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(reflect.V[1] / init.V[1],
		            (Z0 - Z) / (Z + Z0),
		            1e-2);
//		TODO - complete stable beautiful tests
//		ASSERT_NEAR(transfer.sigma(1, 1) / init.sigma(1, 1),
//		            2 * Z / (Z + Z0),
//		            1e-2);
//		ASSERT_NEAR(transfer.V[1] / init.V[1],
//		            2 * Z0 / (Z + Z0),
//		            1e-2);


	}
}


TEST(Engine, TwoLayersDifferentE)
{
	real E2E0Initial = 0.25;
	for (int i = 0; i < 5; i++) {

		Task task;
		Statement statement;
		task.borderSize = 3;
		task.dimensionality = 2;
		statement.CourantNumber = 1.5;

		real rho0 = 1, lambda0 = 2, mu0 = 0.8;
		real rho2rho0 = 1, lambda2lambda0 = E2E0Initial * pow(2, i), mu2mu0 = E2E0Initial * pow(2, i);
		real rho = rho2rho0 * rho0, lambda = lambda2lambda0 * lambda0, mu = mu2mu0 * mu0;

		statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(rho0, lambda0, mu0);
		Statement::MaterialCondition::Inhomogenity newMaterial;
		newMaterial.area = std::make_shared<AxisAlignedBoxArea>(linal::Vector3({-10, 0.5 - 1e-5, -10}), linal::Vector3({10, 10, 10}));
		newMaterial.material = std::make_shared<IsotropicMaterial>(rho, lambda, mu);
		statement.materialConditions.materials.push_back(newMaterial);

		task.sizes(0) = 50;
		task.sizes(1) = 100;
		task.lengthes = {2, 1, 1};

		statement.numberOfSnaps = 0;
		statement.T = 0.24;
		task.enableSnapshotting = true;

		Statement::InitialCondition::Wave wave;
		wave.waveType = Waves::T::P_FORWARD;
		wave.direction = 1; // along y
		wave.quantity = PhysicalQuantities::T::Vy;
		wave.quantityValue = -2;
		linal::Vector3 min({ -1, 0.015, -1});
		linal::Vector3 max({ 4, 0.455, 1});
		wave.area = std::make_shared<AxisAlignedBoxArea>(min, max);
		statement.initialCondition.waves.push_back(wave);

		task.statements.push_back(statement);
		
		EngineWrapper<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>> engine;
		engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>());
		engine.addSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>());
		engine.initialize(task);
		engine.beforeStatementForTest(statement);

		int leftNodeIndex = (int) (task.sizes(1) * 0.25);
		auto init = engine.getSolverForTest()->getMesh()->pde({task.sizes(0) / 2, leftNodeIndex, 0});

		engine.runStatementForTest();

//		int rightNodeIndex = (int) (task.sizes(1) * 0.7);
		auto reflect = engine.getSolverForTest()->getMesh()->pde({task.sizes(0) / 2, leftNodeIndex, 0});
//		auto transfer = engine.getSolverForTest()->getMesh()->pde({task.sizes(0) / 2, rightNodeIndex, 0});

		real E0 = mu0 * (3 * lambda0 + 2 * mu0) / (lambda0 + mu0); // Young's modulus
		real Z0 = sqrt(E0 * rho0); // acoustic impedance
		real E = mu * (3 * lambda + 2 * mu) / (lambda + mu); // Young's modulus
		real Z = sqrt(E * rho); // acoustic impedance

		ASSERT_NEAR(reflect.sigma(1, 1) / init.sigma(1, 1),
		            (Z - Z0) / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(reflect.V[1] / init.V[1],
		            (Z0 - Z) / (Z + Z0),
		            1e-2);

//		ASSERT_NEAR(transfer.sigma(1, 1) / init.sigma(1, 1),
//		            2 * Z / (Z + Z0),
//		            1e-2);
//		ASSERT_NEAR(transfer.V[1] / init.V[1],
//		            2 * Z0 / (Z + Z0),
//		            1e-2);


	}
}
