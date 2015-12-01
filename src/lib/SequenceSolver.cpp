#include "lib/SequenceSolver.hpp"


SequenceSolver::SequenceSolver(Mesh *mesh, Mesh *newMesh) :
		mesh(mesh), newMesh(newMesh) {}


void SequenceSolver::calculate() {
	real currentTime = 0.0; uint step = 0;
	if (makeSnapshots) mesh->snapshot(step);
	while(currentTime < mesh->T) {
		stage(0, mesh->tau);// / 2);
		stage(1, mesh->tau);
//		stage(0, mesh->tau / 2);
		currentTime += mesh->tau; step += 1;
		if (makeSnapshots) mesh->snapshot(step);
	}
}


void SequenceSolver::stage(const uint s, const real& timeStep) {


	for (uint y = 0; y < mesh->Y; y++) {
		// TODO - make all this linal operations faster
		for (uint x = 0; x < mesh->X; x++) {

			// points to interpolate values on previous time layer
			Vector dx = (*mesh)(y, x).matrix->A(s).L.getDiagonalMultipliedBy(timeStep);

			/* new values = U1 * Riemann solvers */
			(*newMesh)(y, x).u = (*mesh)(y, x).matrix->A(s).U1 *
			                     /* Riemann solvers = U * old values */
			                     (*mesh)(y, x).matrix->A(s).U.diagonalMultiply
					                     /* old values are in columns of the matrix */
					                     (mesh->interpolateValuesAround(s, y, x, dx));
		}
	}

	std::swap(mesh, newMesh);

}
