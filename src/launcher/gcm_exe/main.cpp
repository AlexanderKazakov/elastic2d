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
//	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic3DModel, CubicGrid, IsotropicMaterial>>());
//	engine.addSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic3DModel, CubicGrid, IsotropicMaterial>>());
	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>());
	engine.addSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>());
//	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, Cgal2DGrid, IsotropicMaterial>>());
//	engine.addSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic2DModel, Cgal2DGrid, IsotropicMaterial>>());

	try {
		engine.initialize(parseTask2d());

		auto t1 = std::chrono::high_resolution_clock::now();
		engine.run();
		auto t2 = std::chrono::high_resolution_clock::now();

		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		std::cout << "Time of calculation, microseconds = " << duration << std::endl;
	} catch (Exception e) {
		std::cout << e.what() << std::endl;
//		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}

Task parseTaskCgal() {
	Task task;	
	task.spatialStep = 0.4;
	task.enableSnapshotting = true;

	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);

	statement.CourantNumber = 1.0;

	statement.numberOfSnaps = 21;
	statement.stepsPerSnap = 1;

	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.4, linal::Vector3({0.5, 0.5, 0}));
	statement.initialCondition.quantities.push_back(pressure);
	
	statement.quantitiesToVtk = {PhysicalQuantities::T::PRESSURE,
	                          PhysicalQuantities::T::Sxx,
	                          PhysicalQuantities::T::Sxy,
	                          PhysicalQuantities::T::Syy};

	task.statements.push_back(statement);
	return task;
}

Task parseTask2d() {
	Task task;
	task.enableSnapshotting = true;

	task.borderSize = 2;
	task.dimensionality = 2;
	task.lengthes = {4, 2, 1};
	task.sizes = {100, 50, 1};

	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);

	statement.CourantNumber = 0.9;

	statement.numberOfSnaps = 20;
	statement.stepsPerSnap = 1;

	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.2, linal::Vector3({2, 1, 0}));
	statement.initialCondition.quantities.push_back(pressure);

	statement.quantitiesToVtk = {PhysicalQuantities::T::PRESSURE};

	task.statements.push_back(statement);
	return task;
}

Task parseTaskDemo() {
	Task task;
	task.enableSnapshotting = true;

	task.borderSize = 2;
	task.dimensionality = 3;
	task.lengthes = {4, 2, 1};
	task.sizes = {100, 50, 25};

	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);
//	statement.materialConditions.defaultMaterial = std::make_shared<OrthotropicMaterial>
//			(rho, {360, 70, 70, 180, 70, 90, 10, 10, 10}, 1, 1);

	statement.CourantNumber = 0.9;

	statement.numberOfSnaps = 20;
	statement.stepsPerSnap = 1;

	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.2, linal::Vector3({2, 1, 0.5}));
	statement.initialCondition.quantities.push_back(pressure);

	statement.quantitiesToVtk = {PhysicalQuantities::T::PRESSURE};

	task.statements.push_back(statement);
	return task;
}
