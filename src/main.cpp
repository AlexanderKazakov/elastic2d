#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/grid/StructuredGrid.hpp"
#include "lib/util/DataBus.hpp"
#include "lib/util/Logging.hpp"
#include "lib/solver/MpiStructuredSolver.hpp"

using namespace gcm;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();
	USE_AND_INIT_LOGGER("gcm.main");

	StructuredGrid<IdealElastic2DModel> mesh1;
	StructuredGrid<IdealElastic2DModel> mesh2;
	Task task;
	task.enableSnapshotting = true;
	MpiStructuredSolver<IdealElastic2DModel> solver;

	try {
		solver.initialize(task, &mesh1, &mesh2);
		solver.calculate();
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}
