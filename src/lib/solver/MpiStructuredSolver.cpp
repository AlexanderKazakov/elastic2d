#include "lib/model/IdealElastic1DModel.hpp"
#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/model/IdealElastic3DModel.hpp"
#include "MpiStructuredSolver.hpp"
#include "lib/util/DataBus.hpp"

using namespace gcm;

template<class TModel>
void MpiStructuredSolver<TModel>::initialize(const Task &task, StructuredGrid<TModel> *mesh, StructuredGrid<TModel> *newMesh) {

	this->mesh = mesh;
	this->mesh->initialize(task);
	this->newMesh = newMesh;
	this->newMesh->initialize(task);
	CourantNumber = task.CourantNumber;
	tau = CourantNumber * mesh->getMinimalSpatialStep() / mesh->getMaximalLambda(); // time step
	T = task.numberOfSnaps * tau; // required time
	if (task.numberOfSnaps == 0) T = task.T;
	splittingSecondOrder = task.splittingSecondOrder;
	snapshotter = new VtkTextStructuredSnapshotter<TModel>();
	snapshotter->initialize(mesh, task.enableSnapshotting);
}

template<class TModel>
void MpiStructuredSolver<TModel>::calculate() {
	tau = CourantNumber * mesh->getMinimalSpatialStep() / mesh->getMaximalLambda(); // time step
	exchangeNodesWithNeighbors();
	real currentTime = 0.0; int step = 0;
	snapshotter->snapshot(step);
	while(currentTime < T) {
		if (splittingSecondOrder) {
			switch (TModel::DIMENSIONALITY) {
				case 1:
					stage(0, tau);
					break;
				case 2:
					stage(0, tau / 2);
					stage(1, tau);
					stage(0, tau / 2);
					break;
				case 3:
					THROW_UNSUPPORTED("TODO splitting second order in 3D");
					break;
			}
		} else {
			for (int s = 0; s < TModel::DIMENSIONALITY; s++) {
				stage(s, tau);
			}
		}
		currentTime += tau; step += 1;
		snapshotter->snapshot(step);
	}
}

template<class TModel>
void MpiStructuredSolver<TModel>::stage(const int s, const real& timeStep) {

	mesh->applyBorderConditions();

	for (int y = 0; y < mesh->Y; y++) {
		for (int x = 0; x < mesh->X; x++) {

			// points to interpolate values on previous time layer
			auto dx = ((*mesh)(y, x)).matrix->A(s).L.getDiagonalMultipliedBy( - timeStep);

			/* new values = U1 * Riemann solvers */
			(*newMesh)(y, x) = (*mesh)(y, x).matrix->A(s).U1 *
			                     /* Riemann solvers = U * old values */
			                     (*mesh)(y, x).matrix->A(s).U.diagonalMultiply
					                     /* old values are in columns of the matrix */
					                     (mesh->interpolateValuesAround(s, y, x, dx));
		}
	}

	std::swap(mesh, newMesh);
	exchangeNodesWithNeighbors();
}

template<class TModel>
void MpiStructuredSolver<TModel>::exchangeNodesWithNeighbors() {

	int rank = mesh->getRank();
	int numberOfWorkers = mesh->getNumberOfWorkers();
	if (numberOfWorkers == 1) return;

	int sizeOfBuffer = mesh->accuracyOrder * (mesh->X + 2 * mesh->accuracyOrder);
	int nodesSize = (mesh->X + 2 * mesh->accuracyOrder) * (mesh->Y + 2 * mesh->accuracyOrder);

	if (rank == 0) {
		MPI_Sendrecv(&(mesh->nodes[nodesSize - 2 * sizeOfBuffer]), sizeOfBuffer, TModel::Node::MPI_NODE_TYPE, rank + 1, 1,
		             &(mesh->nodes[nodesSize - sizeOfBuffer]), sizeOfBuffer, TModel::Node::MPI_NODE_TYPE, rank + 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);

	} else if (rank == numberOfWorkers - 1) {
		MPI_Sendrecv(&(mesh->nodes[sizeOfBuffer]), sizeOfBuffer, TModel::Node::MPI_NODE_TYPE, rank - 1, 1,
		             mesh->nodes, sizeOfBuffer, TModel::Node::MPI_NODE_TYPE, rank - 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);

	} else {
		MPI_Sendrecv(&(mesh->nodes[nodesSize - 2 * sizeOfBuffer]), sizeOfBuffer, TModel::Node::MPI_NODE_TYPE, rank + 1, 1,
			     	 &(mesh->nodes[nodesSize - sizeOfBuffer]), sizeOfBuffer, TModel::Node::MPI_NODE_TYPE, rank + 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(&(mesh->nodes[sizeOfBuffer]), sizeOfBuffer, TModel::Node::MPI_NODE_TYPE, rank - 1, 1,
		             mesh->nodes, sizeOfBuffer, TModel::Node::MPI_NODE_TYPE, rank - 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

//template class MpiStructuredSolver<IdealElastic1DModel>;
template class MpiStructuredSolver<IdealElastic2DModel>;
template class MpiStructuredSolver<IdealElastic3DModel>;