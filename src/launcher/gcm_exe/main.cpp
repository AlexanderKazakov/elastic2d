#include <lib/util/DataBus.hpp>
#include <lib/Engine.hpp>
#include <lib/util/areas/AxisAlignedBoxArea.hpp>

using namespace gcm;

Task parseTask();


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();
	USE_AND_INIT_LOGGER("gcm.main");

	Engine engine;

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

	task.lengthes = {2, 2, 2};
	task.sizes = {50, 50, 50};

	real rho0 = 8.0; // default density
	real lambda0 = 12e+4; // default Lame parameter
	real mu0 = 77e+3; // default Lame parameter
	task.material = IsotropicMaterial(rho0, lambda0, mu0);

	task.CourantNumber = 1.2; // number from Courant–Friedrichs–Lewy condition
	task.numberOfSnaps = 20;

	/*Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 2.0;
	pressure.area = std::make_shared<SphereArea>(0.2, linal::Vector3({1, 1, 2}));
	task.initialCondition.quantities.push_back(pressure);*/

	Task::InitialCondition::Wave wave;
	wave.waveType = Waves::T::P_BACKWARD;
	wave.direction = 0;
	wave.quantity = PhysicalQuantities::T::PRESSURE;
	wave.quantityValue = 2.0;
	wave.area = std::make_shared<AxisAlignedBoxArea>(linal::Vector3({0.2, -1, -1}), linal::Vector3({0.5, 3, 3}));
	task.initialCondition.waves.push_back(wave);

	task.borderConditions.at(CUBIC_BORDERS::X_LEFT) = BorderCondition::T::FREE_BORDER;

	task.enableSnapshotting = true;

	return task;
}
