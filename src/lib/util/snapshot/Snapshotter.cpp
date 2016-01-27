#include <lib/util/snapshot/Snapshotter.hpp>
#include <lib/grid/nodes/Node.hpp>

using namespace gcm;

template <class TNode>
void Snapshotter<TNode>::snapshot(Grid<TNode>* _grid, const int step) {
	if (enableSnapshotting) {
		this->grid = _grid;
		snapshotImpl(step);
	}
}

template <class TNode>
void Snapshotter<TNode>::initialize(const Task& task) {
		this->enableSnapshotting = task.enableSnapshotting;
}

template class Snapshotter<IdealElastic1DNode>;
template class Snapshotter<IdealElastic2DNode>;
template class Snapshotter<IdealElastic3DNode>;