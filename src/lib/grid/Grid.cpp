#include <fstream>
#include <algorithm>

#include "lib/model/IdealElastic1DModel.hpp"
#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/model/IdealElastic3DModel.hpp"
#include "lib/grid/StructuredGrid.hpp"

using namespace gcm;

void Grid::initialize(const Task &task) {

	rank = MPI::COMM_WORLD.Get_rank();
	numberOfWorkers = MPI::COMM_WORLD.Get_size();

	if (task.forceSequence) {
		rank = 0;
		numberOfWorkers = 1;
	}

	initializeImpl(task);

	/* ------------------ Properties and conditions ------------------ */

	initialConditions = task.initialConditions;
	borderConditions = task.borderConditions;

	/* ------------------ Properties and conditions (end) ------------------ */

	applyInitialConditions();
}
