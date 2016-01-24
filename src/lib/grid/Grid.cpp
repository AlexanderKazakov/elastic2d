#include <lib/grid/Grid.hpp>
#include <lib/grid/nodes/Node.hpp>


using namespace gcm;

template<class TNode>
void Grid<TNode>::initialize(const Task &task) {

	rank = MPI::COMM_WORLD.Get_rank();
	numberOfWorkers = MPI::COMM_WORLD.Get_size();

	if (task.forceSequence) {
		rank = 0;
		numberOfWorkers = 1;
	}

	initializeImpl(task);

	auto gcmMatricesPtr = std::make_shared<typename TNode::GcmMatrices>(task.material);
	for (auto& node : nodes) {
		node.matrix = gcmMatricesPtr;
	}
	maximalLambda = gcmMatricesPtr->getMaximalEigenvalue();

	applyInitialConditions(task);
}

template class Grid<IdealElastic1DNode>;
template class Grid<IdealElastic2DNode>;
template class Grid<IdealElastic3DNode>;