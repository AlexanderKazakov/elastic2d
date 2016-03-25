#include <chrono>

#include <lib/Engine.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/util/snapshot/VtkSnapshotter.hpp>
#include <lib/util/areas/areas.hpp>


using namespace gcm;

Task parseTaskCgal2d();
Task parseTask2d();
Task parseTask3d();


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	USE_AND_INIT_LOGGER("gcm.main");

	try {
		auto t1 = std::chrono::high_resolution_clock::now();
		Engine(parseTaskCgal2d()).run();
		auto t2 = std::chrono::high_resolution_clock::now();

		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		LOG_INFO("Time of calculation, microseconds = " << duration);
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}

Task parseTaskCgal2d() {
	Task task;
	
	task.modelId = Models::T::ELASTIC2D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CGAL2D;
	task.snapshottersId = {Snapshotters::T::VTK};
	
	task.cgal2DGrid.spatialStep = 0.4;

	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);

	statement.globalSettings.CourantNumber = 1.0;

	statement.globalSettings.numberOfSnaps = 21;
	statement.globalSettings.stepsPerSnap = 1;

	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.4, Real3({0.5, 0.5, 0}));
	statement.initialCondition.quantities.push_back(pressure);
	
	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {
			PhysicalQuantities::T::PRESSURE,
			PhysicalQuantities::T::Sxx,
			PhysicalQuantities::T::Sxy,
			PhysicalQuantities::T::Syy
	};

	task.statements.push_back(statement);
	return task;
}

Task parseTask2d() {
	Task task;
	
	task.modelId = Models::T::ELASTIC2D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {Snapshotters::T::VTK};

	task.cubicGrid.borderSize = 2;
	task.cubicGrid.dimensionality = 2;
	task.cubicGrid.lengths = {4, 2, 1};
	task.cubicGrid.sizes = {100, 50, 1};

	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);

	statement.globalSettings.CourantNumber = 0.9;

	statement.globalSettings.numberOfSnaps = 20;
	statement.globalSettings.stepsPerSnap = 1;

	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.2, Real3({2, 1, 0}));
	statement.initialCondition.quantities.push_back(pressure);

	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {PhysicalQuantities::T::PRESSURE};

	task.statements.push_back(statement);
	return task;
}

Task parseTask3d() {
	Task task;
	
	task.modelId = Models::T::ELASTIC3D;
	task.materialId = Materials::T::ISOTROPIC;
	task.gridId = Grids::T::CUBIC;
	task.snapshottersId = {Snapshotters::T::VTK};
	
	task.cubicGrid.borderSize = 2;
	task.cubicGrid.dimensionality = 3;
	task.cubicGrid.lengths = {4, 2, 1};
	task.cubicGrid.sizes = {100, 50, 25};

	Statement statement;
	real rho = 4;
	real lambda = 2;
	real mu = 1;
	statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);
//	statement.materialConditions.defaultMaterial = std::make_shared<OrthotropicMaterial>
//			(rho, {360, 70, 70, 180, 70, 90, 10, 10, 10}, 1, 1);

	statement.globalSettings.CourantNumber = 0.9;

	statement.globalSettings.numberOfSnaps = 20;
	statement.globalSettings.stepsPerSnap = 1;

	Statement::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 10.0;
	pressure.area = std::make_shared<SphereArea>(0.2, Real3({2, 1, 0.5}));
	statement.initialCondition.quantities.push_back(pressure);

	statement.vtkSnapshotter.enableSnapshotting = true;
	statement.vtkSnapshotter.quantitiesToSnap = {PhysicalQuantities::T::PRESSURE};

	task.statements.push_back(statement);
	return task;
}
