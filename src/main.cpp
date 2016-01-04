#include <lib/model/IdealElastic2DModel.hpp>
#include "lib/grid/StructuredGrid.hpp"
#include "lib/util/DataBus.hpp"
#include "lib/solver/MPISolver.hpp"

using namespace gcm;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();

	StructuredGrid<IdealElastic2DModel> mesh1;
	StructuredGrid<IdealElastic2DModel> mesh2;
	Task task;
	mesh1.initialize(task);
	mesh1.changeRheology2(4, 1, 1);
	mesh2.initialize(task);
	mesh2.changeRheology2(4, 1, 1);
	MPISolver<IdealElastic2DModel> solver(&mesh1, &mesh2);
	solver.splittingSecondOrder = true;
	solver.makeSnapshots = true;

	solver.calculate();

	MPI_Finalize();
	return 0;
}
