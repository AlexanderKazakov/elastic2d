#include <lib/grid/Grid.hpp>


using namespace gcm;

void Grid::initialize(const Task &task) {
	rank = MPI::COMM_WORLD.Get_rank();
	numberOfWorkers = MPI::COMM_WORLD.Get_size();

	if (task.forceSequence) {
		rank = 0;
		numberOfWorkers = 1;
	}

	initializeImpl(task);
	applyInitialConditions(task);
}
