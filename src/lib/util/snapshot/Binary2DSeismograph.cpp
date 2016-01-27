#include <lib/util/snapshot/Binary2DSeismograph.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TGrid>
void Binary2DSeismograph<TGrid>::startSeismo(TGrid* grid) {
	static_assert(TGrid::DIMENSIONALITY == 2, "This is seismograph for 2D");
	LOG_DEBUG("Start seismo writing to " << makeFileNameForSeismo());
	openSeismoFileStream(makeFileNameForSeismo(grid));
	int sizeY = grid->sizes(1);
	surface = new real[sizeY];
	seismoNumber++;
}

template<class TGrid>
void Binary2DSeismograph<TGrid>::finishSeismo() {
	closeSeismoFileStream();
	delete [] surface;
}

template<class TGrid>
void Binary2DSeismograph<TGrid>::writeNextTimeStep(TGrid* grid) {
	int sizeY = grid->sizes(1);
	for (int y = 0; y < sizeY; y++) {
		surface[y] = grid->get(0, y, 0).u.getPressure();
	}

	auto bufferSize = (std::streamsize) (sizeY * sizeof(real));
	auto previousNumberOfBytes = seismoFileStream.tellp();
	assert_true(seismoFileStream.write(reinterpret_cast<char*>(surface), bufferSize));
	auto currentNumberOfBytes = seismoFileStream.tellp();
	assert_eq(bufferSize, currentNumberOfBytes - previousNumberOfBytes);
}

template<class TGrid>
void Binary2DSeismograph<TGrid>::openSeismoFileStream(const std::string& fileName) {
	seismoFileStream.open(fileName, std::ios::binary);
	assert_true(seismoFileStream.is_open());
}

template<class TGrid>
void Binary2DSeismograph<TGrid>::closeSeismoFileStream() {
	seismoFileStream.close();
}

template <class TGrid>
std::string Binary2DSeismograph<TGrid>::makeFileNameForSeismo(TGrid* grid) {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d%s", "seismo/core", grid->getRank(), "_seismo", seismoNumber, ".vtk");
	return std::string(buffer);
}



template class Binary2DSeismograph<StructuredGrid<Elastic2DModel>>;
template class Binary2DSeismograph<StructuredGrid<PlasticFlow2DModel>>;
