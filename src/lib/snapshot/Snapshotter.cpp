#include <lib/snapshot/Snapshotter.hpp>
#include <lib/model/IdealElastic1DModel.hpp>
#include <lib/model/IdealElastic3DModel.hpp>
#include <lib/model/IdealElastic2DModel.hpp>

using namespace gcm;

template <class TModel>
void Snapshotter<TModel>::snapshot(const int step) {
	if (enableSnapshotting) {
		snapshotImpl(makeFileNameForSnapshot(step));
	}
}

template <class TModel>
void Snapshotter<TModel>::initialize(Grid<TModel> *grid, bool enableSnapshotting) {
		this->grid = grid;
		this->enableSnapshotting = enableSnapshotting;
}

template <class TModel>
std::string Snapshotter<TModel>::makeFileNameForSnapshot(const int step) {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d%s", "snaps/core", grid->getRank(), "_snapshot", step, ".vtk");
	return std::string(buffer);
}


template class Snapshotter<IdealElastic1DModel>;
template class Snapshotter<IdealElastic2DModel>;
template class Snapshotter<IdealElastic3DModel>;