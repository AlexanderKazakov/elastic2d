#include <lib/util/DataBus.hpp>
#include <lib/Engine.hpp>
#include <lib/util/areas/AxisAlignedBoxArea.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/util/snapshot/VtkTextStructuredSnapshotter.hpp>

using namespace gcm;

Task parseTask();


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();
	USE_AND_INIT_LOGGER("gcm.main");

	Engine engine;
	engine.setSolver(new DefaultSolver<StructuredGrid<SuperDuperModel>>());
	engine.setSnapshotter(new VtkTextStructuredSnapshotter<StructuredGrid<SuperDuperModel>>());

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

	task.accuracyOrder = 2;

	task.lengthes = {4, 2, 1};
	task.sizes = {100, 50, 25};

	real rho = 4; // default density
	real lambda = 2; // default Lame parameter
	real mu = 1; // default Lame parameter
	task.yieldStrength = 1;
	task.continualDamageParameter = 1;
	task.isotropicMaterial = IsotropicMaterial(rho, lambda, mu, task.yieldStrength, task.continualDamageParameter);
	task.orthotropicMaterial = OrthotropicMaterial(rho, {360, 70, 70, 180, 70, 90, 10, 10, 10},
	                                               task.yieldStrength, task.continualDamageParameter);

	task.CourantNumber = 0.9; // number from Courant–Friedrichs–Lewy condition

	task.enableSnapshotting = true;
	task.numberOfSnaps = 20;
	task.stepsPerSnap = 2;

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

	task.borderConditions.at(CUBIC_BORDERS::X_LEFT) = BorderCondition::T::FREE_BORDER;

	return task;
}
