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
	engine.setSolver(new DefaultSolver<StructuredGrid<IdealPlastic2DModel>>());
	engine.setSnapshotter(new VtkTextStructuredSnapshotter<StructuredGrid<IdealPlastic2DModel>>());

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

	task.accuracyOrder = 3;

	task.lengthes = {4, 2, 2};
	task.sizes = {40, 10, 1};

	real rho0 = 4; // default density
	real lambda0 = 2; // default Lame parameter
	real mu0 = 1; // default Lame parameter
	real yieldStrength0 = 1;
	real continualDamageParameter0 = 1;
	task.material = IsotropicMaterial(rho0, lambda0, mu0, yieldStrength0);
	task.material.continualDamageParameter = continualDamageParameter0;
	task.plasticityFlowCorrector = true;

	task.CourantNumber = 1.9; // number from Courant–Friedrichs–Lewy condition

	task.enableSnapshotting = true;
	task.numberOfSnaps = 10;
	task.stepsPerSnap = 3;

	/*Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 2.0;
	pressure.area = std::make_shared<SphereArea>(0.2, linal::Vector3({1, 1, 2}));
	task.initialCondition.quantities.push_back(pressure);*/

	Task::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_FORWARD;
	wave.direction = 0;
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = 10.0;
	wave.area = std::make_shared<AxisAlignedBoxArea>(linal::Vector3({0.2, -1, -1}), linal::Vector3({0.5, 3, 3}));
	task.initialCondition.waves.push_back(wave);

	task.borderConditions.at(CUBIC_BORDERS::X_LEFT) = BorderCondition::T::FREE_BORDER;

	return task;
}
