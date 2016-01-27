#include <lib/numeric/gcmethod/GridCharacteristicMethod.hpp>

using namespace gcm;

template<class TNode>
void GridCharacteristicMethod<TNode>::stage(const int s, const real& timeStep,
                                            const StructuredGrid<TNode>* mesh, StructuredGrid<TNode>* newMesh) {
	LOG_DEBUG("Start stage " << s << " timeStep = " << timeStep);

	for (int x = 0; x < mesh->sizes(0); x++) {
		for (int y = 0; y < mesh->sizes(1); y++) {
			for (int z = 0; z < mesh->sizes(2); z++) {

				// points to interpolate values on previous time layer
				auto dx = - timeStep * linal::diag(mesh->get(x, y, z).matrix->A(s).L);

				/* new values = U1 * Riemann solvers */
				(*newMesh)(x, y, z).u = mesh->get(x, y, z).matrix->A(s).U1 *
				                      /* Riemann solvers = U * old values */
				                      mesh->get(x, y, z).matrix->A(s).U.diagonalMultiply
						                     /* old values are in columns of the matrix */
						                     (mesh->interpolateValuesAround(s, x, y, z, dx));
			}
		}
	}
}


template class GridCharacteristicMethod<IdealElastic1DNode>;
template class GridCharacteristicMethod<IdealElastic2DNode>;
template class GridCharacteristicMethod<IdealElastic3DNode>;