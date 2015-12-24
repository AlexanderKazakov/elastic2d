#include "lib/MPISolver.hpp"

#include "lib/DataBus.hpp"

using namespace gcm;

MPISolver::MPISolver(Mesh *mesh, Mesh *newMesh) :
		mesh(mesh), newMesh(newMesh) {}


void MPISolver::calculate() {
	exchangeNodesWithNeighbors();
	real currentTime = 0.0; int step = 0;
	if (makeSnapshots) mesh->snapshot(step);
	while(currentTime < mesh->T) {
	if (splittingSecondOrder) {
		stage(0, mesh->tau / 2);
		stage(1, mesh->tau);
		stage(0, mesh->tau / 2);
	} else {
		stage(0, mesh->tau);
		stage(1, mesh->tau);
	}
		currentTime += mesh->tau; step += 1;
		if (makeSnapshots) mesh->snapshot(step);
	}
}


void MPISolver::stage(const int s, const real& timeStep) {

	mesh->applyBorderConditions();

	for (int y = 0; y < mesh->Y; y++) {
		for (int x = 0; x < mesh->X; x++) {

			// points to interpolate values on previous time layer
			Vector dx = (*mesh)(y, x).matrix->A(s).L.getDiagonalMultipliedBy( - timeStep);

			/* new values = U1 * Riemann solvers */
			(*newMesh)(y, x).u = (*mesh)(y, x).matrix->A(s).U1 *
			                     /* Riemann solvers = U * old values */
			                     (*mesh)(y, x).matrix->A(s).U.diagonalMultiply
					                     /* old values are in columns of the matrix */
					                     (mesh->interpolateValuesAround(s, y, x, dx));
		}
	}

	std::swap(mesh, newMesh);
	exchangeNodesWithNeighbors();
}


void MPISolver::exchangeNodesWithNeighbors() {

	int rank = mesh->rank;
	int numberOfWorkers = mesh->numberOfWorkers;
	if (numberOfWorkers == 1) return;

	int sizeOfBuffer = mesh->accuracyOrder * (mesh->X + 2 * mesh->accuracyOrder);
	int nodesSize = (mesh->X + 2 * mesh->accuracyOrder) * (mesh->Y + 2 * mesh->accuracyOrder);

	if (rank == 0) {
		MPI_Sendrecv(&(mesh->nodes[nodesSize - 2 * sizeOfBuffer]), sizeOfBuffer, DataBus::MPI_NODE, rank + 1, 1,
		             &(mesh->nodes[nodesSize - sizeOfBuffer]), sizeOfBuffer, DataBus::MPI_NODE, rank + 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
	} else if (rank == numberOfWorkers - 1) {
		MPI_Sendrecv(&(mesh->nodes[sizeOfBuffer]), sizeOfBuffer, DataBus::MPI_NODE, rank - 1, 1,
		             mesh->nodes, sizeOfBuffer, DataBus::MPI_NODE, rank - 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
	} else {
		MPI_Sendrecv(&(mesh->nodes[nodesSize - 2 * sizeOfBuffer]), sizeOfBuffer, DataBus::MPI_NODE, rank + 1, 1,
			     	 &(mesh->nodes[nodesSize - sizeOfBuffer]), sizeOfBuffer, DataBus::MPI_NODE, rank + 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(&(mesh->nodes[sizeOfBuffer]), sizeOfBuffer, DataBus::MPI_NODE, rank - 1, 1,
		             mesh->nodes, sizeOfBuffer, DataBus::MPI_NODE, rank - 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
	}

}