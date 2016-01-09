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
