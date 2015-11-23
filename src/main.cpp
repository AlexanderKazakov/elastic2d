#include <iostream>

using namespace std;

#include "lib/SequenceSolver.hpp"
#include "lib/Mesh.hpp"
#include "lib/Task.hpp"

int main() {
	Mesh mesh1;
	Mesh mesh2;
	Task task;
	mesh1.initialize(task);
	mesh2.initialize(task);
	SequenceSolver sequenceSolver(&mesh1, &mesh2);
	sequenceSolver.makeSnapshots = true;

	sequenceSolver.calculate();

	return 0;
}