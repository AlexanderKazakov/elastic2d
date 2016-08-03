#include <lib/util/snapshot/Binary2DSeismograph.hpp>
#include <lib/rheology/models/models.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/mesh/DefaultMesh.hpp>

#include "lib/numeric/solvers/Solver.hpp"

using namespace gcm;


template<class TMesh>
void Binary2DSeismograph<TMesh>::
beforeStatementImpl(const Statement& statement) {
	auto quantityToWrite = statement.binary2DSeismograph.quantityToWrite;
	valuesGetter = PdeVariables::QUANTITIES.at(quantityToWrite).Get;
	FileUtils::openBinaryFileStream(fileStream, makeFileNameForSnapshot
			("?", -1, FILE_EXTENSION, FOLDER_NAME));
	firstCallInTheStatement = true;
}


template<class TMesh>
void Binary2DSeismograph<TMesh>::
afterStatement() {
	FileUtils::closeFileStream(fileStream);
}


template<class TMesh>
void Binary2DSeismograph<TMesh>::
snapshotImpl(const AbstractGrid* mesh_, const int) {
	const TMesh* mesh = dynamic_cast<const TMesh*>(mesh_);
	assert_true(mesh);
	
	if (firstCallInTheStatement) {
		firstCallInTheStatement = false;
		surface.resize((size_t)mesh->sizes(1) + 1); // plus one for auxiliary data
		writeHeadOfTable(mesh);
	}
	
    assert_eq(mesh->sizes(1), (int)surface.size() - 1);
	surface[0] = (precision) Clock::Time();
	for (size_t y = 0; y < surface.size() - 1; y++) {
		surface[y + 1] = (precision) valuesGetter(mesh->pde({0, (int)y}));
	}
	FileUtils::writeStdVectorToBinaryFileStream(fileStream, surface);
}


template<class TMesh>
void Binary2DSeismograph<TMesh>::
writeHeadOfTable(const TMesh* mesh) {
	surface[0] = (precision) surface.size() - 1;
	for (size_t y = 0; y < surface.size() - 1; y++) {
		surface[y + 1] = (precision) (mesh->coords({0, (int)y})(1));
	}
	FileUtils::writeStdVectorToBinaryFileStream(fileStream, surface);
}


template class Binary2DSeismograph<DefaultMesh<ElasticModel<2>, CubicGrid<2>, IsotropicMaterial> >;
template class Binary2DSeismograph<DefaultMesh<ElasticModel<2>, CubicGrid<2>, OrthotropicMaterial> >;
template class Binary2DSeismograph<DefaultMesh<AcousticModel<2>, CubicGrid<2>, IsotropicMaterial> >;
