#include <lib/util/DataBus.hpp>
#include <lib/numeric/gcmethod/MpiStructuredSolver.hpp>

using namespace gcm;

Task parseTask();


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();
	USE_AND_INIT_LOGGER("gcm.main");

	MpiStructuredSolver<IdealElastic3DNode> solver;

	try {
		solver.initialize(parseTask());
		solver.calculate();
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}


Task parseTask() {
	Task task;

	task.accuracyOrder = 2;

	task.lengthes = {2, 2, 5};
	task.sizes = {21, 21, 51};

	real rho0 = 8.0; // default density
	real lambda0 = 12e+4; // default Lame parameter
	real mu0 = 77e+3; // default Lame parameter
	task.material = IsotropicMaterial(rho0, lambda0, mu0);

	task.CourantNumber = 1.2; // number from Courant–Friedrichs–Lewy condition
	task.numberOfSnaps = 20;

	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 2.0;
	pressure.area = std::make_shared<SphereArea>(0.2, linal::Vector3({1, 1, 2}));
	task.initialCondition.quantities.push_back(pressure);

	task.borderConditions.at(CUBIC_BORDERS::X_LEFT) = BorderCondition::T::FREE_BORDER;

	task.enableSnapshotting = true;

	return task;
}
