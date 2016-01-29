#include <lib/util/snapshot/Snapshotter.hpp>

using namespace gcm;

void Snapshotter::snapshot(const Grid* _grid, const int step) {
	if (enableSnapshotting) {
		this->grid = _grid;
		snapshotImpl(step);
	}
}

void Snapshotter::initialize(const Task& task) {
		this->enableSnapshotting = task.enableSnapshotting;
}
