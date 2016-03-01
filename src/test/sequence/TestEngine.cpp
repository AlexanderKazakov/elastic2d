#include <gtest/gtest.h>

#include <lib/util/areas/areas.hpp>

#include <test/wrappers/Wrappers.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

TEST(Engine, runStatementForTest)
{
	Task task;
	Statement statement;
	task.accuracyOrder = 5;
	statement.CourantNumber = 4.5;
	statement.isotropicMaterial = IsotropicMaterial(4.0, 2.0, 0.5);
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

	EngineWrapper<DefaultMesh<Elastic2DModel, CubicGrid>> engine;
	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid>>());
	engine.initialize(task);
	engine.beforeStatementForTest(statement);
	auto sWave = engine.getSolverForTest()->getMesh()->getPde(task.sizes(0) / 2, 3, 0);
	engine.runStatementForTest();
	ASSERT_EQ(sWave, engine.getSolverForTest()->getMesh()->getPde(task.sizes(0) / 2, 22, 0));
}


TEST(Engine, TwoLayersDifferentRho)
{
	real rho2rho0Initial = 0.25;
	int numberOfSnapsInitial = 30;
	for (int i = 0; i < 5; i++) {

		Task task;
		Statement statement;
		task.accuracyOrder = 3;
		statement.CourantNumber = 1.5;
		statement.isotropicMaterial = IsotropicMaterial(1.0, 2.0, 0.8);
		task.sizes(0) = 50;
		task.sizes(1) = 100;
		task.lengthes = {2, 1, 1};
		statement.numberOfSnaps = numberOfSnapsInitial + 2 * i; // in order to catch the impulses

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

		EngineWrapper<DefaultMesh<Elastic2DModel, CubicGrid>> engine;
		engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid>>());
		engine.initialize(task);
		engine.beforeStatementForTest(statement);

		real rho2rho0 = rho2rho0Initial * pow(2, i);
		real lambda2lambda0 = 1;
		real mu2mu0 = 1;
		engine.getSolverForTest()->getMesh()->changeRheology(rho2rho0, lambda2lambda0, mu2mu0);

		int leftNodeIndex = (int) (task.sizes(1) * 0.25);
		auto init = engine.getSolverForTest()->getMesh()->getPde(task.sizes(0) / 2, leftNodeIndex, 0);

		engine.runStatementForTest();

		int rightNodeIndex = (int) (task.sizes(1) * 0.7);
		auto reflect = engine.getSolverForTest()->getMesh()->getPde(task.sizes(0) / 2, leftNodeIndex, 0);
		auto transfer = engine.getSolverForTest()->getMesh()->getPde(task.sizes(0) / 2, rightNodeIndex, 0);

		real rho0 = statement.isotropicMaterial.rho;
		real lambda0 = statement.isotropicMaterial.lambda;
		real mu0 = statement.isotropicMaterial.mu;
		real E0 = mu0 * (3 * lambda0 + 2 * mu0) / (lambda0 + mu0); // Young's modulus
		real Z0 = sqrt(E0 * rho0); // acoustic impedance


		real rho = rho2rho0 * rho0;
		real lambda = lambda2lambda0 * lambda0;
		real mu = mu2mu0 * mu0;
		real E = mu * (3 * lambda + 2 * mu) / (lambda + mu); // Young's modulus
		real Z = sqrt(E * rho); // acoustic impedance

		ASSERT_NEAR(reflect.sigma(1, 1) / init.sigma(1, 1),
		            (Z - Z0) / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(reflect.V[1] / init.V[1],
		            (Z0 - Z) / (Z + Z0),
		            1e-2);

		ASSERT_NEAR(transfer.sigma(1, 1) / init.sigma(1, 1),
		            2 * Z / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(transfer.V[1] / init.V[1],
		            2 * Z0 / (Z + Z0),
		            1e-2);


	}
}


TEST(Engine, TwoLayersDifferentE)
{
	real E2E0Initial = 0.25;
	int numberOfSnapsInitial = 40;
	for (int i = 0; i < 5; i++) {

		Task task;
		Statement statement;
		task.accuracyOrder = 3;
		statement.CourantNumber = 1.5;
		statement.isotropicMaterial = IsotropicMaterial(1.0, 2.0, 0.8);
		task.sizes(0) = 50;
		task.sizes(1) = 100;
		task.lengthes = {2, 1, 1};
		statement.numberOfSnaps = numberOfSnapsInitial - 2 * i; // in order to catch the impulses

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
		
		EngineWrapper<DefaultMesh<Elastic2DModel, CubicGrid>> engine;
		engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid>>());
		engine.initialize(task);
		engine.beforeStatementForTest(statement);

		real rho2rho0 = 1;
		real lambda2lambda0 = E2E0Initial * pow(2, i);
		real mu2mu0 = E2E0Initial * pow(2, i);
		engine.getSolverForTest()->getMesh()->changeRheology(rho2rho0, lambda2lambda0, mu2mu0);

		int leftNodeIndex = (int) (task.sizes(1) * 0.25);
		auto init = engine.getSolverForTest()->getMesh()->getPde(task.sizes(0) / 2, leftNodeIndex, 0);

		engine.runStatementForTest();

		int rightNodeIndex = (int) (task.sizes(1) * 0.7);
		auto reflect = engine.getSolverForTest()->getMesh()->getPde(task.sizes(0) / 2, leftNodeIndex, 0);
		auto transfer = engine.getSolverForTest()->getMesh()->getPde(task.sizes(0) / 2, rightNodeIndex, 0);

		real rho0 = statement.isotropicMaterial.rho;
		real lambda0 = statement.isotropicMaterial.lambda;
		real mu0 = statement.isotropicMaterial.mu;
		real E0 = mu0 * (3 * lambda0 + 2 * mu0) / (lambda0 + mu0); // Young's modulus
		real Z0 = sqrt(E0 * rho0); // acoustic impedance


		real rho = rho2rho0 * rho0;
		real lambda = lambda2lambda0 * lambda0;
		real mu = mu2mu0 * mu0;
		real E = mu * (3 * lambda + 2 * mu) / (lambda + mu); // Young's modulus
		real Z = sqrt(E * rho); // acoustic impedance

		ASSERT_NEAR(reflect.sigma(1, 1) / init.sigma(1, 1),
		            (Z - Z0) / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(reflect.V[1] / init.V[1],
		            (Z0 - Z) / (Z + Z0),
		            1e-2);

		ASSERT_NEAR(transfer.sigma(1, 1) / init.sigma(1, 1),
		            2 * Z / (Z + Z0),
		            1e-2);
		ASSERT_NEAR(transfer.V[1] / init.V[1],
		            2 * Z0 / (Z + Z0),
		            1e-2);


	}
}