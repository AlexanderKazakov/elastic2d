#include <lib/numeric/gcmethod/GridCharacteristicMethod.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TGrid>
void GridCharacteristicMethod<TGrid>::stage(const int s, const real& timeStep,
                                            const TGrid* mesh, TGrid* newMesh) {
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


template class GridCharacteristicMethod<DefaultStructuredGrid<Elastic1DModel>>;
template class GridCharacteristicMethod<DefaultStructuredGrid<Elastic2DModel>>;
template class GridCharacteristicMethod<DefaultStructuredGrid<Elastic3DModel>>;
template class GridCharacteristicMethod<DefaultStructuredGrid<OrthotropicElastic3DModel>>;
