#include <lib/Engine.hpp>
#include <lib/util/areas/areas.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/util/snapshot/Binary2DSeismograph.hpp>

using namespace gcm;

Task parseTask();


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	USE_AND_INIT_LOGGER("gcm.seismo");
	Engine engine;
	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>());
	engine.addSnapshotter(new Binary2DSeismograph<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>());

	try {
		engine.initialize(parseTask());
		engine.run();
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}


Task parseTask() {
	Task task;
	task.forceSequence = true;

	task.borderSize = 3;
	task.dimensionality = 2;
	task.lengthes = {1, 1, 1};
	task.sizes = {50, 50, 1};
	task.enableSnapshotting = true;	
	
	Statement statement;
	statement.CourantNumber = 1.2; // number from Courant–Friedrichs–Lewy condition

	real rho = 4; // default density
	real lambda = 2; // default Lame parameter
	real mu = 1; // default Lame parameter
	statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(rho, lambda, mu, 1, 1);

	statement.numberOfSnaps = 100;
	statement.stepsPerSnap = 1;

	Statement::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_FORWARD;
	wave.direction = 1;
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = MPI::COMM_WORLD.Get_rank() + 1;
	wave.area = std::make_shared<AxisAlignedBoxArea>(linal::Vector3({-1, 0.2, -1}), linal::Vector3({3, 0.5, 3}));
	statement.initialCondition.waves.push_back(wave);

	statement.id = "0000";
	task.statements.push_back(statement);
	
	Statement::BorderCondition borderCondition;	
	// x right free border
	borderCondition.area = std::make_shared<AxisAlignedBoxArea>
		(linal::Vector3({-10, 0.99, -10}), linal::Vector3({10, 10, 10}));
	borderCondition.values = {
		{PhysicalQuantities::T::Sxy, [](real){return 0;}},
		{PhysicalQuantities::T::Syy, [](real){return 0;}}
	};
	statement.borderConditions.push_back(borderCondition);
	
	statement.id = "0001";
	task.statements.push_back(statement);
	
	return task;
}
