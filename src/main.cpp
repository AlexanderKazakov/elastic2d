#include <lib/model/IdealElastic2DModel.hpp>
#include "lib/mesh/Mesh.hpp"
#include "lib/DataBus.hpp"
#include "lib/MPISolver.hpp"

using namespace gcm;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();

	Mesh<IdealElastic2DModel> mesh1;
	Mesh<IdealElastic2DModel> mesh2;
	Task task;
	mesh1.initialize(task);
	mesh1.changeRheology2(4, 1, 1);
	mesh2.initialize(task);
	mesh2.changeRheology2(4, 1, 1);
	MPISolver solver(&mesh1, &mesh2);
	solver.splittingSecondOrder = true;
	solver.makeSnapshots = true;

	solver.calculate();

	MPI_Finalize();
	return 0;
}
