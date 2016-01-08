#include "lib/snapshot/Snapshotter.hpp"

using namespace gcm;


void Snapshotter::snapshot(const int step) {
	if (enableSnapshotting) {
		snapshotImpl(makeFileNameForSnapshot(step));
	}
}

void Snapshotter::initialize(Grid *grid, bool enableSnapshotting) {
		this->grid = grid;
		this->enableSnapshotting = enableSnapshotting;
}

std::string Snapshotter::makeFileNameForSnapshot(const int step) {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d%s", "snaps/core", MPI::COMM_WORLD.Get_rank(), "_snapshot", step, ".vtk");
	return std::string(buffer);
}