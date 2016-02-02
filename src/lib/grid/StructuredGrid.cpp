#include <fstream>
#include <algorithm>

#include <lib/grid/StructuredGrid.hpp>
#include <lib/util/task/InitialCondition.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<typename TModel, template<typename> class GcmMatricesStorage>
void StructuredGrid<TModel, GcmMatricesStorage>::initializeImpl(const Task &task) {
	LOG_INFO("Start initialization");
	accuracyOrder = task.accuracyOrder; // order of accuracy of spatial interpolation
	assert_ge(accuracyOrder, 1);

	h = linal::plainDivision(task.lengthes, task.sizes - linal::VectorInt<3>({1, 1, 1}));
	sizes = task.sizes;
	startR = task.startR;

	// MPI - we divide the grid among processes equally along x-axis
	int numberOfNodesAlongXPerOneCore = (int) std::round((real) task.sizes(0) / this->numberOfWorkers);
	sizes(0) = numberOfNodesAlongXPerOneCore; // number of nodes along x direction on this mesh
	// in order to keep specified in task number of nodes
	if (this->rank == this->numberOfWorkers - 1) sizes(0) = task.sizes(0) - numberOfNodesAlongXPerOneCore * (this->numberOfWorkers - 1);
	startR(0) += (this->rank * numberOfNodesAlongXPerOneCore) * h(0); // mpi parallel along X axis
	globalStartXindex = this->rank * numberOfNodesAlongXPerOneCore; // for testing
	// MPI - (end)

	for (int j = 0; j < 3; j++) {
		if (sizes(j) != 1) assert_ge(sizes(j), 2 * accuracyOrder);
		assert_gt(h(j), 0.0);
		assert_eq(h(j), h(j)); // this is supposed to catch NaN
		assert_eq(startR(j), startR(j));
	}
	
	borderConditions.initialize(task);

	this->nodes.resize( (unsigned long) (sizes(0) + 2 * accuracyOrder) * (sizes(1) + 2 * accuracyOrder) * (sizes(2) + 2 * accuracyOrder) );

	typename TModel::Material material;
	material.initialize(task);
	auto gcmMatricesPtr = std::make_shared<GCM_MATRICES>(material);
	for (auto& node : nodes) {
		node.matrix = gcmMatricesPtr;
	}
	maximalLambda = gcmMatricesPtr->getMaximalEigenvalue();
	minimalSpatialStep = fmin(h(0), fmin(h(1), h(2)));
}

template<typename TModel, template<typename> class GcmMatricesStorage>
void StructuredGrid<TModel, GcmMatricesStorage>::applyInitialConditions(const Task& task) {
	InitialCondition<TModel> initialCondition;
	initialCondition.initialize(task);

	for (int x = 0; x < sizes(0); x++) {
		for (int y = 0; y < sizes(1); y++) {
			for (int z = 0; z < sizes(2); z++) {
				initialCondition.apply((*this)(x, y, z).u, getCoordinates(x, y, z));
			}
		}
	}
}

template<typename TModel, template<typename> class GcmMatricesStorage>
typename StructuredGrid<TModel, GcmMatricesStorage>::Matrix StructuredGrid<TModel, GcmMatricesStorage>::interpolateValuesAround
		(const int stage, const int x, const int y, const int z, const Vector& dx) const {

	Matrix ans;
	std::vector<Vector> src( (unsigned long) (accuracyOrder + 1) );
	Vector res;
	for (int k = 0; k < Vector::M; k++) {
		findSourcesForInterpolation(stage, x, y, z, dx(k), src);
		interpolator.minMaxInterpolate(res, src, fabs(dx(k)) / h(stage));
		ans.setColumn(k, res);
	}

	return ans;
}

template<typename TModel, template<typename> class GcmMatricesStorage>
void StructuredGrid<TModel, GcmMatricesStorage>::findSourcesForInterpolation(const int stage, const int x, const int y, const int z,
                                                         const real &dx, std::vector<Vector>& src) const {

	const int alongX = (stage == 0) * ( (dx > 0) ? 1 : -1 );
	const int alongY = (stage == 1) * ( (dx > 0) ? 1 : -1 );
	const int alongZ = (stage == 2) * ( (dx > 0) ? 1 : -1 );
	for (int k = 0; k < src.size(); k++) {
		src[(unsigned long)k] = get(x + alongX * k, y + alongY * k, z + alongZ * k).u;
	}
}

template<typename TModel, template<typename> class GcmMatricesStorage>
void StructuredGrid<TModel, GcmMatricesStorage>::beforeStageImpl() {
	exchangeNodesWithNeighbors();
	applyBorderConditions();
}

template<typename TModel, template<typename> class GcmMatricesStorage>
void StructuredGrid<TModel, GcmMatricesStorage>::exchangeNodesWithNeighbors() {
	LOG_DEBUG("Start data exchange with neighbor cores");
	if (this->numberOfWorkers == 1) return;

	int sizeOfBuffer = this->accuracyOrder * (sizes(1) + 2 * this->accuracyOrder) * (sizes(2) + 2 * this->accuracyOrder);
	unsigned long nodesSize = this->nodes.size();

	if (this->rank == 0) {
		MPI_Sendrecv(&(this->nodes[nodesSize - 2 * sizeOfBuffer]), sizeOfBuffer, NODE::MPI_NODE_TYPE, this->rank + 1, 1,
		             &(this->nodes[nodesSize - sizeOfBuffer]), sizeOfBuffer, NODE::MPI_NODE_TYPE, this->rank + 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);

	} else if (this->rank == this->numberOfWorkers - 1) {
		MPI_Sendrecv(&(this->nodes[(unsigned long)sizeOfBuffer]), sizeOfBuffer, NODE::MPI_NODE_TYPE, this->rank - 1, 1,
		             &(this->nodes[0]), sizeOfBuffer, NODE::MPI_NODE_TYPE, this->rank - 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);

	} else {
		MPI_Sendrecv(&(this->nodes[nodesSize - 2 * sizeOfBuffer]), sizeOfBuffer, NODE::MPI_NODE_TYPE, this->rank + 1, 1,
		             &(this->nodes[nodesSize - sizeOfBuffer]), sizeOfBuffer, NODE::MPI_NODE_TYPE, this->rank + 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(&(this->nodes[(unsigned long)sizeOfBuffer]), sizeOfBuffer, NODE::MPI_NODE_TYPE, this->rank - 1, 1,
		             &(this->nodes[0]), sizeOfBuffer, NODE::MPI_NODE_TYPE, this->rank - 1, 1,
		             MPI::COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

template<typename TModel, template<typename> class GcmMatricesStorage>
void StructuredGrid<TModel, GcmMatricesStorage>::applyBorderConditions() {
	borderConditions.applyBorderConditions(this);
}

template class StructuredGrid<Elastic1DModel>;

template class StructuredGrid<Elastic2DModel>;
template class StructuredGrid<ContinualDamageElastic2DModel>;
template class StructuredGrid<IdealPlastic2DModel>;

template class StructuredGrid<Elastic3DModel>;
template class StructuredGrid<OrthotropicElastic3DModel>;

template class StructuredGrid<SuperDuperModel>;

