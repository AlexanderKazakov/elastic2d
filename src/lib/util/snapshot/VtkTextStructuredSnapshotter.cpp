#include <lib/util/snapshot/VtkTextStructuredSnapshotter.hpp>

using namespace gcm;

template<class TNode>
void VtkTextStructuredSnapshotter<TNode>::snapshotImpl(const int step) {
	LOG_DEBUG("Start snapshot writing to " << fileName);
	openSnapshotFileStream(makeFileNameForSnapshot(step));
	sGrid = static_cast<StructuredGrid<TNode>*>(this->grid);
	
	snapshotFileStream << "# vtk DataFile Version 3.0" << std::endl;
	snapshotFileStream << "U data" << std::endl;
	snapshotFileStream << "ASCII" << std::endl;
	snapshotFileStream << "DATASET STRUCTURED_POINTS" << std::endl;
	snapshotFileStream << "DIMENSIONS " << sGrid->sizes(0) << " " << sGrid->sizes(1) << " " << sGrid->sizes(2) << std::endl;
	snapshotFileStream << "SPACING " << sGrid->h(0) << " " << sGrid->h(1) << " " << sGrid->h(2) << std::endl;
	snapshotFileStream << "ORIGIN " << sGrid->startR(0) << " " << sGrid->startR(1) << " " << sGrid->startR(2) << std::endl;
	snapshotFileStream << "POINT_DATA " << sGrid->sizes(0) * sGrid->sizes(1) * sGrid->sizes(2) << std::endl;

	for (auto& quantity : TNode::Vector::QUANTITIES) {
		writeQuantity(PhysicalQuantities::NAME.at(quantity.first), quantity.second.Get);
	}

	closeSnapshotFileStream();
}

template<class TNode>
void VtkTextStructuredSnapshotter<TNode>::writeQuantity(const std::string name, const typename GetSetter<typename TNode::Variables>::Getter Get) {
	snapshotFileStream << "SCALARS " << name << " double" << std::endl;
	snapshotFileStream << "LOOKUP_TABLE default" << std::endl;
	for (int z = 0; z < sGrid->sizes(2); z++) {
		for (int y = 0; y < sGrid->sizes(1); y++) {
			for (int x = 0; x < sGrid->sizes(0); x++) {
				snapshotFileStream << Get(sGrid->get(x, y, z).u) << std::endl;
			}
		}
	}
}

template<class TNode>
void VtkTextStructuredSnapshotter<TNode>::openSnapshotFileStream(const std::string& fileName) {
	snapshotFileStream.open(fileName, std::ios::out);
	assert_true(snapshotFileStream.is_open());
}

template<class TNode>
void VtkTextStructuredSnapshotter<TNode>::closeSnapshotFileStream() {
	snapshotFileStream.close();
}

template <class TNode>
std::string VtkTextStructuredSnapshotter<TNode>::makeFileNameForSnapshot(const int step) {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d%s", "snaps/core", this->grid->getRank(), "_snapshot", step, ".vtk");
	return std::string(buffer);
}

template class VtkTextStructuredSnapshotter<IdealElastic1DNode>;
template class VtkTextStructuredSnapshotter<IdealElastic2DNode>;
template class VtkTextStructuredSnapshotter<IdealElastic3DNode>;