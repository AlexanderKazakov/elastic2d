#include "Binary2DSeismograph.hpp"

using namespace gcm;

template<class TNode>
void Binary2DSeismograph<TNode>::startSeismo(StructuredGrid<TNode>* grid) {
	static_assert(TNode::DIMENSIONALITY == 2, "This is seismograph for 2D");
	LOG_DEBUG("Start seismo writing to " << makeFileNameForSeismo());
	openSeismoFileStream(makeFileNameForSeismo(grid));
	int sizeY = grid->sizes(1);
	surface = new real[sizeY];
	seismoNumber++;
}

template<class TNode>
void Binary2DSeismograph<TNode>::finishSeismo() {
	closeSeismoFileStream();
	delete [] surface;
}

template<class TNode>
void Binary2DSeismograph<TNode>::writeNextTimeStep(StructuredGrid<TNode>* grid) {
	int sizeY = grid->sizes(1);
	for (int y = 0; y < sizeY; y++) {
		surface[y] = grid->get(0, y, 0).getPressure();
	}

	size_t bufferSize = sizeY * sizeof(real);
	auto previousNumberOfBytes = seismoFileStream.tellp();
	assert_true(seismoFileStream.write(reinterpret_cast<char*>(surface), bufferSize));
	auto currentNumberOfBytes = seismoFileStream.tellp();
	assert_eq(bufferSize, currentNumberOfBytes - previousNumberOfBytes);
}

template<class TNode>
void Binary2DSeismograph<TNode>::openSeismoFileStream(const std::string& fileName) {
	seismoFileStream.open(fileName, std::ios::binary);
	assert_true(seismoFileStream.is_open());
}

template<class TNode>
void Binary2DSeismograph<TNode>::closeSeismoFileStream() {
	seismoFileStream.close();
}

template <class TNode>
std::string Binary2DSeismograph<TNode>::makeFileNameForSeismo(StructuredGrid<TNode>* grid) {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d%s", "seismo/core", grid->getRank(), "_seismo", seismoNumber, ".vtk");
	return std::string(buffer);
}