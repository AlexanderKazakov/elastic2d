#include <chrono>

#include <lib/Engine.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/util/snapshot/VtkSnapshotter.hpp>
#include <lib/util/areas/areas.hpp>


using namespace gcm;

Task parseTaskCgal();
Task parseTask2d();
Task parseTaskDemo();


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	USE_AND_INIT_LOGGER("gcm.main");

	Engine engine;
//	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic3DModel, CubicGrid>>());
//	engine.addSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic3DModel, CubicGrid>>());
	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid>>());
	engine.addSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic2DModel, CubicGrid>>());
//	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, Cgal2DGrid>>());
//	engine.addSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic2DModel, Cgal2DGrid>>());

	try {
		engine.initialize(parseTask2d());

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
	
	task.quantitiesToVtk = {PhysicalQuantities::T::PRESSURE,
	                          PhysicalQuantities::T::Sxx,
	                          PhysicalQuantities::T::Sxy,
	                          PhysicalQuantities::T::Syy};

	return task;
}

Task parseTask2d() {
	Task task;

	task.accuracyOrder = 2;

	task.lengthes = {4, 2, 1};
	task.sizes = {100, 50, 1};

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
	pressure.area = std::make_shared<SphereArea>(0.2, linal::Vector3({2, 1, 0}));
	task.initialCondition.quantities.push_back(pressure);

	task.quantitiesToVtk = {PhysicalQuantities::T::PRESSURE};

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

	task.quantitiesToVtk = {PhysicalQuantities::T::PRESSURE};

	return task;
}
