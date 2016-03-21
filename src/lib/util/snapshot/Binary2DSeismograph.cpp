#include <lib/util/snapshot/Binary2DSeismograph.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>

#include "lib/numeric/solvers/Solver.hpp"

using namespace gcm;

template<class TMesh>
void Binary2DSeismograph<TMesh>::initializeImpl(const Task& task) {
	assert_eq(task.cubicGrid.dimensionality, 2);
	sizeY = (size_t) task.cubicGrid.sizes(1);
	hY = task.cubicGrid.h(1);
	surface.resize(sizeY + 1); // plus one for auxiliary gnuplot data
}

template<class TMesh>
void Binary2DSeismograph<TMesh>::beforeStatementImpl(const Statement& statement) {
	auto quantityToWrite = statement.binary2DSeismograph.quantityToWrite;
	valuesGetter = PdeVariables::QUANTITIES.at(quantityToWrite).Get;
	FileUtils::openBinaryFileStream(fileStream, makeFileNameForSnapshot
			(-1, FILE_EXTENSION, FOLDER_NAME));
	writeHeadOfTable();
}

template<class TMesh>
void Binary2DSeismograph<TMesh>::afterStatement() {
	FileUtils::closeFileStream(fileStream);
}

template<class TMesh>
void Binary2DSeismograph<TMesh>::snapshotImpl(const AbstractGrid* mesh_, const int) {
	const TMesh* mesh = dynamic_cast<const TMesh*>(mesh_);
	assert_eq(mesh->sizes(1), sizeY); // TODO slice iterator
	surface[0] = (precision) Engine::Instance().getCurrentTime();
	for (size_t y = 0; y < sizeY; y++) {
		surface[y + 1] = (precision) valuesGetter(mesh->pde(linal::VectorInt<3>({0, (int)y, 0})));
	}
	FileUtils::writeStdVectorToBinaryFileStream(fileStream, surface);
}

template <class TMesh>
void Binary2DSeismograph<TMesh>::writeHeadOfTable() {
	surface[0] = (precision) sizeY;
	for (size_t y = 0; y < sizeY; y++) {
		surface[y + 1] = (precision) (y * hY);
	}
	FileUtils::writeStdVectorToBinaryFileStream(fileStream, surface);
}


template class Binary2DSeismograph<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>>;
