#include "lib/grid/Grid.hpp"

#include "lib/model/IdealElastic1DModel.hpp"
#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/model/IdealElastic3DModel.hpp"

using namespace gcm;

template<class TModel>
void Grid<TModel>::initialize(const Task &task) {

	rank = MPI::COMM_WORLD.Get_rank();
	numberOfWorkers = MPI::COMM_WORLD.Get_size();

	if (task.forceSequence) {
		rank = 0;
		numberOfWorkers = 1;
	}

	initializeImpl(task);

	defaultMatrix = std::make_shared<typename TModel::GcmMatrices>(task.rho0, task.lambda0, task.mu0);
	assert_true(defaultMatrix);
	maximalLambda = defaultMatrix->getMaximalEigenvalue();

	for (auto& node : nodes) {
		linal::clear(node.u);
		node.matrix = defaultMatrix;
	}

	applyInitialConditions(task);
}

template class Grid<IdealElastic1DModel>;
template class Grid<IdealElastic2DModel>;
template class Grid<IdealElastic3DModel>;