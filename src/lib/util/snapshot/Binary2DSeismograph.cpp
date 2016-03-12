#include <lib/util/snapshot/Binary2DSeismograph.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>

#include "lib/numeric/solvers/Solver.hpp"

using namespace gcm;

template<class TMesh>
void Binary2DSeismograph<TMesh>::initializeImpl(const Task& task) {
	static_assert(TMesh::DIMENSIONALITY == 2, "This is seismograph for 2D");
	sizeY = task.sizes(1);
}

template<class TMesh>
void Binary2DSeismograph<TMesh>::beforeStatementImpl(const Statement&) {
	FileUtils::openBinaryFileStream(fileStream, makeFileNameForSnapshot
			(-1, FILE_EXTENSION, FOLDER_NAME));
	surface = new precision[sizeY + 1];  // plus one for auxiliary gnuplot data
	writeHeadOfTable();
}

template<class TMesh>
void Binary2DSeismograph<TMesh>::afterStatementImpl() {
	FileUtils::closeFileStream(fileStream);
	delete [] surface;
}

template<class TMesh>
void Binary2DSeismograph<TMesh>::snapshotImpl(const AbstractGrid* mesh_, const int step) {
	const TMesh* mesh = dynamic_cast<const TMesh*>(mesh_);
	assert_eq(mesh->sizes(1), sizeY);
	surface[0] = step * tau;
	for (int y = 0; y < sizeY; y++) {
		surface[y + 1] = (precision) mesh->pde(linal::VectorInt<3>({0, y, 0})).getPressure();
	}
	FileUtils::writeArrayToBinaryFileStream(fileStream, surface, (size_t)sizeY + 1);
}

template <class TMesh>
void Binary2DSeismograph<TMesh>::writeHeadOfTable() {
	surface[0] = (precision)sizeY;
	for (int i = 0; i < sizeY; i++) {
		surface[i + 1] = i * hY;
	}
	FileUtils::writeArrayToBinaryFileStream(fileStream, surface, (size_t)sizeY + 1);
}


template class Binary2DSeismograph<DefaultMesh<Elastic2DModel, CubicGrid>>;
