#include <lib/model/IdealElastic1DModel.hpp>
#include <lib/model/IdealElastic3DModel.hpp>
#include <lib/model/IdealElastic2DModel.hpp>
#include <lib/snapshot/VtkTextStructuredSnapshotter.hpp>

using namespace gcm;

template<class TModel>
void VtkTextStructuredSnapshotter<TModel>::snapshotImpl(const std::string &fileName) {
	LOG_DEBUG("Start snapshot writing to " << fileName);
	openSnapshotFileStream(fileName);
	sGrid = static_cast<StructuredGrid<TModel>*>(this->grid);
	
	snapshotFileStream << "# vtk DataFile Version 3.0" << std::endl;
	snapshotFileStream << "U data" << std::endl;
	snapshotFileStream << "ASCII" << std::endl;
	snapshotFileStream << "DATASET STRUCTURED_POINTS" << std::endl;
	snapshotFileStream << "DIMENSIONS " << sGrid->X << " " << sGrid->Y << " " << sGrid->Z << std::endl;
	snapshotFileStream << "SPACING " << sGrid->h[0] << " " << sGrid->h[1] << " " << sGrid->h[2] << std::endl;
	snapshotFileStream << "ORIGIN " << sGrid->startX << " " << sGrid->startY << " " << sGrid->startZ << std::endl;
	snapshotFileStream << "POINT_DATA " << sGrid->X * sGrid->Y * sGrid->Z << std::endl;

	for (auto& quantity : TModel::QUANTITIES) {
		writeQuantity(PhysicalQuantities::NAME.at(quantity.first), quantity.second.Get);
	}

	closeSnapshotFileStream();
}

template<class TModel>
void VtkTextStructuredSnapshotter<TModel>::writeQuantity(const std::string name, const typename TModel::Getter Get) {
	snapshotFileStream << "SCALARS " << name << " double" << std::endl;
	snapshotFileStream << "LOOKUP_TABLE default" << std::endl;
	for (int z = 0; z < sGrid->Z; z++) {
		for (int y = 0; y < sGrid->Y; y++) {
			for (int x = 0; x < sGrid->X; x++) {
				snapshotFileStream << Get(sGrid->get(x, y, z)) << std::endl;
			}
		}
	}
}

template<class TModel>
void VtkTextStructuredSnapshotter<TModel>::openSnapshotFileStream(const std::string& fileName) {
	snapshotFileStream.open(fileName, std::ios::out);
	assert_true(snapshotFileStream.is_open());
}

template<class TModel>
void VtkTextStructuredSnapshotter<TModel>::closeSnapshotFileStream() {
	snapshotFileStream.close();
}

template class VtkTextStructuredSnapshotter<IdealElastic1DModel>;
template class VtkTextStructuredSnapshotter<IdealElastic2DModel>;
template class VtkTextStructuredSnapshotter<IdealElastic3DModel>;