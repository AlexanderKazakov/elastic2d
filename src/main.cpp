#include <lib/DataBus.hpp>
#include "lib/MPISolver.hpp"


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();

	Mesh mesh1;
	Mesh mesh2;
	Task task;
	mesh1.initialize(task);
	mesh1.changeRheology2(4, 1, 1);
	mesh2.initialize(task);
	mesh2.changeRheology2(4, 1, 1);
	MPISolver solver(&mesh1, &mesh2);
	solver.splittingSecondOrder = true;
//	solver.makeSnapshots = true;

	solver.calculate();

	MPI_Finalize();
	return 0;
}
