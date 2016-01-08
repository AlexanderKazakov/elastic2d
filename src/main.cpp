#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/grid/StructuredGrid.hpp"
#include "lib/util/DataBus.hpp"
#include "lib/solver/MpiStructuredSolver.hpp"

using namespace gcm;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();

	StructuredGrid<IdealElastic2DModel> mesh1;
	StructuredGrid<IdealElastic2DModel> mesh2;
	Task task;
	MpiStructuredSolver<IdealElastic2DModel> solver;
	solver.initialize(task, &mesh1, &mesh2);

	solver.calculate();

	MPI_Finalize();
	return 0;
}
