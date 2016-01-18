#include "lib/model/IdealElastic1DModel.hpp"
#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/model/IdealElastic3DModel.hpp"
#include "lib/util/DataBus.hpp"
#include "lib/solver/MpiStructuredSolver.hpp"

using namespace gcm;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();
	USE_AND_INIT_LOGGER("gcm.main");

	Task task;

	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 2.0;
	pressure.area = std::make_shared<SphereArea>(0.2, linal::Vector3({1, 1, 0}));
	task.initialCondition.quantities.push_back(pressure);

	task.borderConditions.at(CUBIC_BORDERS::X_LEFT) = BorderCondition::T::FREE_BORDER;

	task.enableSnapshotting = true;
	MpiStructuredSolver<IdealElastic2DModel> solver;

	try {
		solver.initialize(task);
		solver.calculate();
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}
