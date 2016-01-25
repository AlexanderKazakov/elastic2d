#include <lib/numeric/gcmethod/MpiStructuredSolver.hpp>
#include <lib/util/DataBus.hpp>

using namespace gcm;

template<class TNode>
void MpiStructuredSolver<TNode>::initialize(const Task &task) {
	LOG_INFO("Start initialization");
	this->mesh = new StructuredGrid<TNode>();
	this->newMesh = new StructuredGrid<TNode>();
	this->mesh->initialize(task);
	this->newMesh->initialize(task);

	CourantNumber = task.CourantNumber;
	tau = CourantNumber * mesh->getMinimalSpatialStep() / mesh->getMaximalLambda(); // time step
	T = task.numberOfSnaps * tau; // required time
	if (task.numberOfSnaps == 0) T = task.T;
	splittingSecondOrder = task.splittingSecondOrder;

	snapshotter = new VtkTextStructuredSnapshotter<TNode>();
	snapshotter->initialize(task);
}

template<class TNode>
MpiStructuredSolver<TNode>::~MpiStructuredSolver() {
	delete this->mesh;
	delete this->newMesh;
	delete this->snapshotter;
}

template<class TNode>
void MpiStructuredSolver<TNode>::calculate() {
	LOG_INFO("Start calculation");
	tau = CourantNumber * mesh->getMinimalSpatialStep() / mesh->getMaximalLambda(); // time step
	exchangeNodesWithNeighbors();
	real currentTime = 0.0; int step = 0;
	snapshotter->snapshot(mesh, step);
	while(currentTime < T) {
		LOG_INFO("Start time step " << step);
		if (splittingSecondOrder) {
			switch (TNode::DIMENSIONALITY) {
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
			for (int s = 0; s < TNode::DIMENSIONALITY; s++) {
				stage(s, tau);
			}
		}
		currentTime += tau; step += 1;
		snapshotter->snapshot(mesh, step);
	}
}

template<class TNode>
void MpiStructuredSolver<TNode>::stage(const int s, const real& timeStep) {
	LOG_DEBUG("Start stage " << s << " tau = " << timeStep);
	mesh->applyBorderConditions();

	for (int x = 0; x < mesh->sizes(0); x++) {
		for (int y = 0; y < mesh->sizes(1); y++) {
			for (int z = 0; z < mesh->sizes(2); z++) {

				// points to interpolate values on previous time layer
				auto dx = - timeStep * linal::diag((*mesh)(x, y, z).matrix->A(s).L);

				/* new values = U1 * Riemann solvers */
				(*newMesh)(x, y, z).u = (*mesh)(x, y, z).matrix->A(s).U1 *
				                      /* Riemann solvers = U * old values */
				                      (*mesh)(x, y, z).matrix->A(s).U.diagonalMultiply
						                     /* old values are in columns of the matrix */
						                     (mesh->interpolateValuesAround(s, x, y, z, dx));
			}
		}
	}

	std::swap(mesh, newMesh);
	exchangeNodesWithNeighbors();
}

template<class TNode>
void MpiStructuredSolver<TNode>::exchangeNodesWithNeighbors() {
	LOG_DEBUG("Start data exchange with neighbor cores");

	int rank = mesh->getRank();
	int numberOfWorkers = mesh->getNumberOfWorkers();
	if (numberOfWorkers == 1) return;

	int sizeOfBuffer = mesh->accuracyOrder * (mesh->sizes(1) + 2 * mesh->accuracyOrder) * (mesh->sizes(2) + 2 * mesh->accuracyOrder);
	unsigned long nodesSize = mesh->nodes.size();

	if (rank == 0) {
		MPI_Sendrecv(&(mesh->nodes[nodesSize - 2 * sizeOfBuffer]), sizeOfBuffer, TNode::MPI_NODE_TYPE, rank + 1, 1,
		             &(mesh->nodes[nodesSize - sizeOfBuffer]), sizeOfBuffer, TNode::MPI_NODE_TYPE, rank + 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);

	} else if (rank == numberOfWorkers - 1) {
		MPI_Sendrecv(&(mesh->nodes[(unsigned long)sizeOfBuffer]), sizeOfBuffer, TNode::MPI_NODE_TYPE, rank - 1, 1,
		             &(mesh->nodes[0]), sizeOfBuffer, TNode::MPI_NODE_TYPE, rank - 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);

	} else {
		MPI_Sendrecv(&(mesh->nodes[nodesSize - 2 * sizeOfBuffer]), sizeOfBuffer, TNode::MPI_NODE_TYPE, rank + 1, 1,
			     	 &(mesh->nodes[nodesSize - sizeOfBuffer]), sizeOfBuffer, TNode::MPI_NODE_TYPE, rank + 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(&(mesh->nodes[(unsigned long)sizeOfBuffer]), sizeOfBuffer, TNode::MPI_NODE_TYPE, rank - 1, 1,
		             &(mesh->nodes[0]), sizeOfBuffer, TNode::MPI_NODE_TYPE, rank - 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

template class MpiStructuredSolver<IdealElastic1DNode>;
template class MpiStructuredSolver<IdealElastic2DNode>;
template class MpiStructuredSolver<IdealElastic3DNode>;