#include "lib/model/IdealElastic1DModel.hpp"
#include "lib/model/IdealElastic3DModel.hpp"
#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/snapshot/VtkTextStructuredSnapshotter.hpp"

using namespace gcm;

template<class TModel>
void VtkTextStructuredSnapshotter<TModel>::snapshotImpl(const std::string &fileName) {
	LOG_DEBUG("Start snapshot writing to " << fileName);
	openSnapshotFileStream(fileName);
	structuredGrid = static_cast<StructuredGrid<TModel>*>(grid);
	
	snapshotFileStream << "# vtk DataFile Version 3.0" << std::endl;
	snapshotFileStream << "U data" << std::endl;
	snapshotFileStream << "ASCII" << std::endl;
	snapshotFileStream << "DATASET STRUCTURED_POINTS" << std::endl;
	snapshotFileStream << "DIMENSIONS " << structuredGrid->X << " " << structuredGrid->Y << " " << structuredGrid->Z << std::endl;
	snapshotFileStream << "SPACING " << structuredGrid->h[0] << " " << structuredGrid->h[1] << " " << structuredGrid->h[2] << std::endl;
	snapshotFileStream << "ORIGIN " << structuredGrid->startX * structuredGrid->h[0] << " 0 0" << std::endl;
	snapshotFileStream << "POINT_DATA " << structuredGrid->X * structuredGrid->Y * structuredGrid->Z << std::endl;

	for (auto& vec : TModel::Node::VECTORS) {
		writeVector(vec.first, vec.second.first, vec.second.second);
	}

	for (auto& scalar : TModel::Node::SCALARS) {
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
	for (int z = 0; z < structuredGrid->Z; z++) {
		for (int y = 0; y < structuredGrid->Y; y++) {
			for (int x = 0; x < structuredGrid->X; x++) {
				snapshotFileStream << structuredGrid->get(x, y, z)(index) << std::endl;
			}
		}
	}
}

template<class TModel>
void VtkTextStructuredSnapshotter<TModel>::writeVector(const std::string& name, const int index, const int size) {
	// TODO - will it work with floats?
	snapshotFileStream << "VECTORS " << name << " double" << std::endl;
	for (int z = 0; z < structuredGrid->Z; z++) {
		for (int y = 0; y < structuredGrid->Y; y++) {
			for (int x = 0; x < structuredGrid->X; x++) {
				for (int i = index; i < index + size; i++) {
					snapshotFileStream << structuredGrid->get(x, y, z)(i) << " ";
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