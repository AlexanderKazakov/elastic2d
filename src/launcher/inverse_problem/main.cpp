#include <lib/Engine.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/util/snapshot/VtkSnapshotter.hpp>
#include <lib/util/snapshot/Detector.hpp>
#include <lib/util/areas/areas.hpp>


using namespace gcm;

Task parseTaskCagi2d();
Task parseTaskCagi3d();


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	USE_AND_INIT_LOGGER("gcm.main");
	Engine engine;
	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid>>());
	engine.addSnapshotter(new VtkSnapshotter<DefaultMesh<Elastic2DModel, CubicGrid>>());
	engine.addSnapshotter(new Detector<DefaultMesh<Elastic2DModel, CubicGrid>>());

	const int numberOfStatements = 1;
	try {
		for (int i = 0; i < numberOfStatements; i++) {
			engine.initialize(parseTaskCagi2d());
			engine.run();
		}
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}


Task parseTaskCagi2d() {
	Task task;
	task.accuracyOrder = 2;
	task.forceSequence = true;

	real X = 0.016, Y = 0.004;
	task.lengthes = {X, Y, 1};
	task.sizes = {201, 51, 1};

	real rho = 1e+3;
	real lambda = 3e+10;
	real mu = 2e+10;
	task.isotropicMaterial = IsotropicMaterial(rho, lambda, mu);

	task.CourantNumber = 1.0;

	task.enableSnapshotting = true;
	task.numberOfSnaps = 101;
	task.stepsPerSnap = 1;

	// border conditions
	auto area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({X/2-0.001*X, Y - 1e-5, -10}), linal::Vector3({X/2+0.001*X, 10, 10}));
	Task::BorderCondition borderCondition;	
	// y up
	real frequency = 10e+6, T = 1.0 / frequency;
	real A = - 1e+6;
	borderCondition.area = area;
	borderCondition.values = {
		{PhysicalQuantities::T::Syy, [A, T](real t){return (t < T) ? A : 0;}},
		{PhysicalQuantities::T::Sxy, [](real){return 0;}}
	};
	task.borderConditions.push_back(borderCondition);
	// y bottom
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({-10, -10, -10}), linal::Vector3({10, 1e-5, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxy, [](real){return 0;}},
		{PhysicalQuantities::T::Syy, [](real){return 0;}}
	};
	task.borderConditions.push_back(borderCondition);
	

	Task::Fracture fracture;
	fracture.direction = 1;
	fracture.coordinate = Y/2;
	fracture.area = std::make_shared<SphereArea>(Y/2, linal::Vector3({X/2, Y/2, 0}));
	fracture.values = {
		{PhysicalQuantities::T::Sxy, [](real){return 0;}},
		{PhysicalQuantities::T::Syy, [](real){return 0;}}
	};
	task.fractures.push_back(fracture); 
	
	// quantities to snapshot
	task.quantitiesToVtk = {PhysicalQuantities::T::PRESSURE,
	                        PhysicalQuantities::T::Sxx,
	                        PhysicalQuantities::T::Sxy,
	                        PhysicalQuantities::T::Syy};
	
	task.detector.quantities = {PhysicalQuantities::T::Syy};
	task.detector.area = area;
	
	return task;
}


Task parseTaskCagi3d() {
	Task task;
	task.accuracyOrder = 2;
	task.forceSequence = true;

	real X = 0.016, Y = 0.016, Z = 0.004;
	task.lengthes = {X, Y, Z};
	task.sizes = {81, 81, 41};

	real rho = 1e+3;
	real lambda = 3e+10;
	real mu = 2e+10;
	task.isotropicMaterial = IsotropicMaterial(rho, lambda, mu);

	task.CourantNumber = 1.0;

	task.enableSnapshotting = true;
	task.numberOfSnaps = 61;
	task.stepsPerSnap = 1;

	// border conditions
	Task::BorderCondition borderCondition;	
	// z up
	real frequency = 10e+6, T = 1.0 / frequency;
	real A = - 1e+6;
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({-100*0.2*X, -100*0.2*Y, Z - 1e-5}), linal::Vector3({100*0.6*X, 100*0.6*Y, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxz, [](real){return 0;}},
		{PhysicalQuantities::T::Syz, [](real){return 0;}},
		{PhysicalQuantities::T::Szz, [A, T](real t){return (t < T) ? A : 0;}}
	};
	task.borderConditions.push_back(borderCondition);
	// z bottom
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({-10, -10, -10}), linal::Vector3({10, 10, 1e-5}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxz, [](real){return 0;}},
		{PhysicalQuantities::T::Syz, [](real){return 0;}},
		{PhysicalQuantities::T::Szz, [](real){return 0;}}
	};
	task.borderConditions.push_back(borderCondition);
		
	Task::Fracture fracture;
	fracture.direction = 2;
	fracture.coordinate = Z/2;
	fracture.area = std::make_shared<SphereArea>(Z/2, linal::Vector3({X/2, Y/2, Z/2}));
	fracture.values = {
		{PhysicalQuantities::T::Sxz, [](real){return 0;}},
		{PhysicalQuantities::T::Syz, [](real){return 0;}},
		{PhysicalQuantities::T::Szz, [](real){return 0;}}
	};
	task.fractures.push_back(fracture); 
	
	// quantities to snapshot
	task.quantitiesToVtk = {PhysicalQuantities::T::PRESSURE,
	                          PhysicalQuantities::T::Sxx,
	                          PhysicalQuantities::T::Sxy,
	                          PhysicalQuantities::T::Syy};

	
	return task;
}