#include <lib/util/snapshot/VtkTextStructuredSnapshotter.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TGrid>
void VtkTextStructuredSnapshotter<TGrid>::snapshotImpl(const Grid* _grid, const int step) {
	LOG_DEBUG("Start snapshot writing to " << makeFileNameForSnapshot(step));
	grid = static_cast<const TGrid*>(_grid);
	openSnapshotFileStream(makeFileNameForSnapshot(step));

	snapshotFileStream << "# vtk DataFile Version 3.0" << std::endl;
	snapshotFileStream << "U data" << std::endl;
	snapshotFileStream << "ASCII" << std::endl;
	snapshotFileStream << "DATASET STRUCTURED_POINTS" << std::endl;
	snapshotFileStream << "DIMENSIONS " << grid->sizes(0) << " " << grid->sizes(1) << " " << grid->sizes(2) << std::endl;
	snapshotFileStream << "SPACING " << grid->h(0) << " " << grid->h(1) << " " << grid->h(2) << std::endl;
	snapshotFileStream << "ORIGIN " << grid->startR(0) << " " << grid->startR(1) << " " << grid->startR(2) << std::endl;
	snapshotFileStream << "POINT_DATA " << grid->sizes(0) * grid->sizes(1) * grid->sizes(2) << std::endl;

	for (auto& quantity : TGrid::Vector::QUANTITIES) {
		writeQuantity(PhysicalQuantities::NAME.at(quantity.first), quantity.second.Get);
	}

	for (auto& quantity : TGrid::Model::InternalOde::QUANTITIES) {
		writeQuantity(PhysicalQuantities::NAME.at(quantity.first), quantity.second.Get);
	}

	closeSnapshotFileStream();
}

template<class TGrid>
void VtkTextStructuredSnapshotter<TGrid>::writeQuantity(const std::string name,
                      const typename GetSetter<typename TGrid::NODE::Variables>::Getter Get) {
	snapshotFileStream << "SCALARS " << name << " double" << std::endl;
	snapshotFileStream << "LOOKUP_TABLE default" << std::endl;
	for (int z = 0; z < grid->sizes(2); z++) {
		for (int y = 0; y < grid->sizes(1); y++) {
			for (int x = 0; x < grid->sizes(0); x++) {
				snapshotFileStream << Get(grid->get(x, y, z).u) << std::endl;
			}
		}
	}
}

template<class TGrid>
void VtkTextStructuredSnapshotter<TGrid>::writeQuantity(const std::string name,
                                                        const typename GetSetter<typename TGrid::Model::InternalOde::Variables>::Getter Get) {
	snapshotFileStream << "SCALARS " << name << " double" << std::endl;
	snapshotFileStream << "LOOKUP_TABLE default" << std::endl;
	for (int z = 0; z < grid->sizes(2); z++) {
		for (int y = 0; y < grid->sizes(1); y++) {
			for (int x = 0; x < grid->sizes(0); x++) {
				snapshotFileStream << Get(grid->get(x, y, z)) << std::endl;
			}
		}
	}
}

template<class TGrid>
void VtkTextStructuredSnapshotter<TGrid>::openSnapshotFileStream(const std::string& fileName) {
	snapshotFileStream.open(fileName, std::ios::out);
	assert_true(snapshotFileStream.is_open());
}

template<class TGrid>
void VtkTextStructuredSnapshotter<TGrid>::closeSnapshotFileStream() {
	snapshotFileStream.close();
}

template <class TGrid>
std::string VtkTextStructuredSnapshotter<TGrid>::makeFileNameForSnapshot(const int step) {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d%s", "snaps/core", MPI::COMM_WORLD.Get_rank(), "_snapshot", step, ".vtk");
	return std::string(buffer);
}


template class VtkTextStructuredSnapshotter<StructuredGrid<Elastic1DModel>>;
template class VtkTextStructuredSnapshotter<StructuredGrid<Elastic2DModel>>;
template class VtkTextStructuredSnapshotter<StructuredGrid<Elastic3DModel>>;
template class VtkTextStructuredSnapshotter<StructuredGrid<OrthotropicElastic3DModel>>;
template class VtkTextStructuredSnapshotter<StructuredGrid<ContinualDamageElastic2DModel>>;
template class VtkTextStructuredSnapshotter<StructuredGrid<IdealPlastic2DModel>>;

template class VtkTextStructuredSnapshotter<StructuredGrid<SuperDuperModel>>;
