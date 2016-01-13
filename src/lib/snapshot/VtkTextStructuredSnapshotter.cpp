#include "lib/model/IdealElastic1DModel.hpp"
#include "lib/model/IdealElastic3DModel.hpp"
#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/snapshot/VtkTextStructuredSnapshotter.hpp"

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

	for (auto& vec : TModel::VECTORS) {
		writeVector(vec.first, vec.second.first, vec.second.second);
	}

	for (auto& scalar : TModel::SCALARS) {
		writeScalar(scalar.first, scalar.second);
	}

	// TODO - how to write a pressure, for example ?

	closeSnapshotFileStream();
}

template<class TModel>
void VtkTextStructuredSnapshotter<TModel>::writeScalar(const std::string name, const int index) {
	// TODO - will it work with floats?
	snapshotFileStream << "SCALARS " << name << " double" << std::endl;
	snapshotFileStream << "LOOKUP_TABLE default" << std::endl;
	for (int z = 0; z < sGrid->Z; z++) {
		for (int y = 0; y < sGrid->Y; y++) {
			for (int x = 0; x < sGrid->X; x++) {
				snapshotFileStream << sGrid->get(x, y, z).u(index) << std::endl;
			}
		}
	}
}

template<class TModel>
void VtkTextStructuredSnapshotter<TModel>::writeVector(const std::string& name, const int index, const int size) {
	// TODO - will it work with floats?
	snapshotFileStream << "VECTORS " << name << " double" << std::endl;
	for (int z = 0; z < sGrid->Z; z++) {
		for (int y = 0; y < sGrid->Y; y++) {
			for (int x = 0; x < sGrid->X; x++) {
				for (int i = index; i < index + size; i++) {
					snapshotFileStream << sGrid->get(x, y, z).u(i) << " ";
				}
				for (int i = index + size; i < index + VTK_VECTOR_SIZE; i++) {
					snapshotFileStream << "0 "; // vtk recognize only vectors of size 3
				}
				snapshotFileStream << std::endl;
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