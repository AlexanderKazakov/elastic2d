#include <chrono>

#include <lib/Engine.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/util/snapshot/VtkSnapshotter.hpp>
#include <lib/util/areas/areas.hpp>


using namespace gcm;

Task parseTaskCagi();
Task parseTaskCgal();
Task parseTaskDemo();


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	USE_AND_INIT_LOGGER("gcm.main");

	Engine engine;
	/*engine.setSolver(new DefaultSolver<DefaultMesh<SuperDuperModel, CubicGrid>>());
	engine.setSnapshotter(new VtkSnapshotter<DefaultMesh<SuperDuperModel, CubicGrid>>());*/
	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid>>());
	engine.setSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic2DModel, CubicGrid>>());
//	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, Cgal2DGrid>>());
//	engine.setSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic2DModel, Cgal2DGrid>>());

	try {
		engine.initialize(parseTaskCagi());

		auto t1 = std::chrono::high_resolution_clock::now();
		engine.run();
		auto t2 = std::chrono::high_resolution_clock::now();

		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		std::cout << "Time of calculation, microseconds = " << duration << std::endl;
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}


Task parseTaskCagi() {
	Task task;

	task.accuracyOrder = 2;

	real X = 0.016 / 2, Y = 0.004;
	task.lengthes = {X, Y, 1};
	task.sizes = {81, 41, 1};

	real rho = 1e+3;
	real lambda = 3e+10;
	real mu = 2e+10;
	task.isotropicMaterial = IsotropicMaterial(rho, lambda, mu);

	task.CourantNumber = 1.0;

	task.enableSnapshotting = true;
	task.numberOfSnaps = 121;
	task.stepsPerSnap = 1;

	// border conditions
	Task::BorderCondition borderCondition;
	// x right free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({X - 1e-5, -10, -10}), linal::Vector3({10, 10, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxx, [](real){return 0;}},
		{PhysicalQuantities::T::Sxy, [](real){return 0;}}
	};
	task.borderConditions.push_back(borderCondition);
	
	// y up sigma
	real frequency = 10e+6, T = 1.0 / frequency;
	real A = - 1e+6;
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({-10, Y - 1e-5, -10}), linal::Vector3({10, 10, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxy, [A, T](real t){return (t < T) ? A : 0;}},
		{PhysicalQuantities::T::Syy, [](real){return 0;}}
	};
	task.borderConditions.push_back(borderCondition);
	// y bottom fixed border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({-10, -10, -10}), linal::Vector3({10, 1e-5, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Vx, [](real){return 0;}},
		{PhysicalQuantities::T::Vy, [](real){return 0;}}
	};
	task.borderConditions.push_back(borderCondition);

	
//	Task::InitialCondition::Wave wave;
//	wave.waveType = Waves::T::P_FORWARD;
//	wave.direction = 0;
//	wave.quantity = PhysicalQuantities::T::PRESSURE;
//	wave.quantityValue = 10.0;
//	wave.area = std::make_shared<AxisAlignedBoxArea>(linal::Vector3({0.2*X, -10, -10}), linal::Vector3({0.4*X, 10, 10}));
//	task.initialCondition.waves.push_back(wave);
	
	// quantities to snapshot
	task.quantitiesToWrite = {PhysicalQuantities::T::PRESSURE,
	                          PhysicalQuantities::T::Sxx,
	                          PhysicalQuantities::T::Sxy,
	                          PhysicalQuantities::T::Syy};

	return task;
}

Task parseTaskCgal() {
	Task task;
	
	task.spatialStep = 0.4;

	real rho = 4;
	real lambda = 2;
	real mu = 1;
	task.isotropicMaterial = IsotropicMaterial(rho, lambda, mu, 1, 1);
	task.orthotropicMaterial = OrthotropicMaterial(rho, {360, 70, 70, 180, 70, 90, 10, 10, 10}, 1, 1);

	task.CourantNumber = 1.0;

	task.enableSnapshotting = true;
	task.numberOfSnaps = 21;
	task.stepsPerSnap = 1;

	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.4, linal::Vector3({0.5, 0.5, 0}));
	task.initialCondition.quantities.push_back(pressure);
	
	task.quantitiesToWrite = {PhysicalQuantities::T::PRESSURE,
	                          PhysicalQuantities::T::Sxx,
	                          PhysicalQuantities::T::Sxy,
	                          PhysicalQuantities::T::Syy};

	return task;
}

Task parseTaskDemo() {
	Task task;

	task.accuracyOrder = 2;

	task.lengthes = {4, 2, 1};
	task.sizes = {100, 50, 25};

	real rho = 4;
	real lambda = 2;
	real mu = 1;
	task.isotropicMaterial = IsotropicMaterial(rho, lambda, mu, 1, 1);
	task.orthotropicMaterial = OrthotropicMaterial(rho, {360, 70, 70, 180, 70, 90, 10, 10, 10}, 1, 1);

	task.CourantNumber = 0.9;

	task.enableSnapshotting = true;
	task.numberOfSnaps = 20;
	task.stepsPerSnap = 1;

	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.2, linal::Vector3({2, 1, 0.5}));
	task.initialCondition.quantities.push_back(pressure);

	/*Task::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_FORWARD;
	wave.direction = 0;
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = 10.0;
	wave.area = std::make_shared<AxisAlignedBoxArea>(linal::Vector3({0.2, -1, -1}), linal::Vector3({0.5, 3, 3}));
	task.initialCondition.waves.push_back(wave);*/

	task.quantitiesToWrite = {PhysicalQuantities::T::PRESSURE};

	return task;
}
