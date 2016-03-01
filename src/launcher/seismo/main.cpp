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
	engine.setSolver(new DefaultSolver<DefaultMesh<Elastic2DModel, CubicGrid>>());
	engine.addSnapshotter(new Binary2DSeismograph<DefaultMesh<Elastic2DModel, CubicGrid>>());

	const int numberOfStatements = 1;
	try {
		for (int i = 0; i < numberOfStatements; i++) {
			engine.initialize(parseTask());
			engine.run();
		}
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}


Task parseTask() {
	Task task;
	task.forceSequence = true;
	task.accuracyOrder = 3;

	task.lengthes = {1, 1, 1};
	task.sizes = {50, 50, 1};
	task.CourantNumber = 1.2; // number from Courant–Friedrichs–Lewy condition

	real rho0 = 4; // default density
	real lambda0 = 2; // default Lame parameter
	real mu0 = 1; // default Lame parameter
	real yieldStrength0 = 1;
	task.isotropicMaterial = IsotropicMaterial(rho0, lambda0, mu0, yieldStrength0);

	task.enableSnapshotting = true;
	task.numberOfSnaps = 100;
	task.stepsPerSnap = 1;


	Task::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_FORWARD;
	wave.direction = 1;
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = MPI::COMM_WORLD.Get_rank() + 1;
	wave.area = std::make_shared<AxisAlignedBoxArea>(linal::Vector3({-1, 0.2, -1}), linal::Vector3({3, 0.5, 3}));
	task.initialCondition.waves.push_back(wave);

	return task;
}
