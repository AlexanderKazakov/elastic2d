#include <lib/util/snapshot/Binary2DSeismograph.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/util/areas/AxisAlignedBoxArea.hpp>

using namespace gcm;

Task nextStatement(int number);


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	USE_AND_INIT_LOGGER("inverse_problem.main");
	DefaultSolver<StructuredGrid<Elastic2DModel>> solver;
	Binary2DSeismograph<StructuredGrid<Elastic2DModel>> seismograph;
	const int numberOfStatements = 2;

	try {
		for (int i = 0; i < numberOfStatements; i++) {
			solver.initialize(nextStatement(i + 1));
			solver.nextTimeStep();
		}
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}


Task nextStatement(int number) {
	Task task;
	task.forceSequence = true;
	task.accuracyOrder = 2;

	task.lengthes = {1, 1, 0};
	task.sizes = {50, 50, 1};

	real rho0 = 8.0; // default density
	real lambda0 = 12e+4; // default Lame parameter
	real mu0 = 77e+3; // default Lame parameter
	task.material = IsotropicMaterial(rho0, lambda0, mu0);

	task.CourantNumber = 1.2; // number from Courant–Friedrichs–Lewy condition
	task.numberOfSnaps = 100;

	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = (MPI::COMM_WORLD.Get_rank() % 2 ? 1 : -1) * number;
	pressure.area = std::make_shared<SphereArea>(0.1, linal::Vector3({0.5, 0.5, 0}));
	task.initialCondition.quantities.push_back(pressure);

	return task;
}
