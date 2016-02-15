#include <lib/util/snapshot/Binary2DSeismograph.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TGrid>
void Binary2DSeismograph<TGrid>::startSeismo(const Task& task) {
	static_assert(TGrid::DIMENSIONALITY == 2, "This is seismograph for 2D");
	LOG_DEBUG("Start seismo writing to " << makeFileNameForSeismo());
	openSeismoFileStream(makeFileNameForSeismo());

	sizeY = task.sizes(1);
	hY = 1; // TODO
	tau = 1; //TODO
	surface = new output_precision[sizeY + 1];  // plus one for auxiliary gnuplot data
	writeHeadOfTable();
	seismoNumber++;
}

template<class TGrid>
void Binary2DSeismograph<TGrid>::finishSeismo() {
	closeSeismoFileStream();
	delete [] surface;
}

template<class TGrid>
void Binary2DSeismograph<TGrid>::snapshotImpl(const Grid* _grid, const int step) {
	const TGrid* grid = static_cast<const TGrid*>(_grid);
	assert_eq(grid->sizes(1), sizeY);
	surface[0] = step * tau;
	for (int y = 0; y < sizeY; y++) {
		surface[y + 1] = (output_precision) grid->pde(0, y, 0).getPressure();
	}
	writeSurfaceToBuffer();
}

template<class TGrid>
void Binary2DSeismograph<TGrid>::openSeismoFileStream(const std::string& fileName) {
	seismoFileStream.open(fileName, std::ios::binary);
	assert_true(seismoFileStream.is_open());
	seismoFileStream.clear();
}

template<class TGrid>
void Binary2DSeismograph<TGrid>::writeSurfaceToBuffer() {
	auto bufferSize = (std::streamsize) ((sizeY + 1) * sizeof(output_precision));
	auto previousNumberOfBytes = seismoFileStream.tellp();
	assert_true(seismoFileStream.write(reinterpret_cast<char*>(surface), bufferSize));
	auto currentNumberOfBytes = seismoFileStream.tellp();
	assert_eq(bufferSize, currentNumberOfBytes - previousNumberOfBytes);
}

template<class TGrid>
void Binary2DSeismograph<TGrid>::closeSeismoFileStream() {
	assert_true(seismoFileStream.good());
	seismoFileStream.close();
}

template <class TGrid>
std::string Binary2DSeismograph<TGrid>::makeFileNameForSeismo() {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d%s", "seismo/core", MPI::COMM_WORLD.Get_rank(), "_seismo", seismoNumber, ".bin");
	return std::string(buffer);
}

template <class TGrid>
void Binary2DSeismograph<TGrid>::writeHeadOfTable() {
	surface[0] = (output_precision) sizeY;
	for (int i = 0; i < sizeY; i++) {
		surface[i + 1] = i * hY;
	}
	writeSurfaceToBuffer();
}

template <class TGrid>
void Binary2DSeismograph<TGrid>::initializeImpl(const Task &task) {
	SUPPRESS_WUNUSED(task);
}


template class Binary2DSeismograph<StructuredGrid<Elastic2DModel>>;
