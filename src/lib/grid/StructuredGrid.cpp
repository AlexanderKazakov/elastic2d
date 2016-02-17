#include <fstream>
#include <algorithm>

#include <lib/grid/StructuredGrid.hpp>
#include <lib/util/task/InitialCondition.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;


void StructuredGrid::initializeImpl(const Task &task) {
	LOG_INFO("Start initialization");
	accuracyOrder = task.accuracyOrder; // order of accuracy of spatial interpolation
	assert_ge(accuracyOrder, 1);

	h = linal::plainDivision(task.lengthes, task.sizes - linal::VectorInt<3>({1, 1, 1}));
	this->minimalSpatialStep = fmin(h(0), fmin(h(1), h(2)));
	sizes = task.sizes;
	startR = task.startR;

	// MPI - we divide the grid among processes equally along x-axis
	int numberOfNodesAlongXPerOneCore = (int) std::round((real) task.sizes(0) / this->numberOfWorkers);
	sizes(0) = numberOfNodesAlongXPerOneCore; // number of nodes along x direction on this mesh
	// in order to keep specified in task number of nodes
	if (this->rank == this->numberOfWorkers - 1) sizes(0) = task.sizes(0) - numberOfNodesAlongXPerOneCore * (this->numberOfWorkers - 1);
	startR(0) += (this->rank * numberOfNodesAlongXPerOneCore) * h(0); // mpi parallel along X axis
	// MPI - (end)

	for (int j = 0; j < 3; j++) {
		if (sizes(j) != 1) assert_ge(sizes(j), 2 * accuracyOrder);
		assert_gt(h(j), 0.0);
		assert_eq(h(j), h(j)); // this is supposed to catch NaN
		assert_eq(startR(j), startR(j));
	}
//	borderConditions.initialize(task);

	initializeImplImpl(task);
}


void StructuredGrid::beforeStageImpl() {
	exchangeNodesWithNeighbors();
//	borderConditions.applyBorderConditions(this);
}


void StructuredGrid::exchangeNodesWithNeighbors() {
	// TODO - move it to DataBus
	/*LOG_DEBUG("Start data exchange with neighbor cores");
	if (this->numberOfWorkers == 1) return;

	int bufferSize = this->accuracyOrder * (sizes(1) + 2 * this->accuracyOrder) * (sizes(2) + 2 * this->accuracyOrder);
	size_t size = this->pdeVectors.size();

	if (this->rank == 0) {
		MPI_Sendrecv(&(this->pdeVectors[size - 2 * bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, this->rank + 1, 1,
		             &(this->pdeVectors[size - bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, this->rank + 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);

	} else if (this->rank == this->numberOfWorkers - 1) {
		MPI_Sendrecv(&(this->pdeVectors[(size_t)bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, this->rank - 1, 1,
		             &(this->pdeVectors[0]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, this->rank - 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);

	} else {
		MPI_Sendrecv(&(this->pdeVectors[size - 2 * bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, this->rank + 1, 1,
		             &(this->pdeVectors[size - bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, this->rank + 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(&(this->pdeVectors[(size_t)bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, this->rank - 1, 1,
		             &(this->pdeVectors[0]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, this->rank - 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
	}*/
}