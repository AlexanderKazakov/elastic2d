#ifndef LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHOD_HPP

#include <lib/linal/linal.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/Types.hpp>

#include <lib/grid/Cgal2DGrid.hpp>
#include <lib/rheology/models/Model.hpp>

namespace gcm {
	template<class TGrid>
	class GridCharacteristicMethod {
	public:
		/**
		 * Do grid-characteristic stage of splitting method
		 * @param s direction
		 * @param timeStep time step
		 * @param mesh grid on current time layer
		 * @param newMesh grid on next time layer
		 */
		void stage(const int s, const real &timeStep, TGrid* mesh) {
			LOG_DEBUG("Start stage " << s << " timeStep = " << timeStep);

			for (int x = 0; x < mesh->sizes(0); x++) {
				for (int y = 0; y < mesh->sizes(1); y++) {
					for (int z = 0; z < mesh->sizes(2); z++) {

						// points to interpolate values on previous time layer
						auto dx = - timeStep * linal::diag(mesh->getMatrix(x, y, z)->A(s).L);

						mesh->pdeVectorsNew[mesh->getIndex(x, y, z)] =
								/* new values = U1 * Riemann solvers */
								mesh->getMatrix(x, y, z)->A(s).U1 *
								/* Riemann solvers = U * old values */
								mesh->getMatrix(x, y, z)->A(s).U.diagonalMultiply
										/* old values are in columns of the matrix */
										(mesh->interpolateValuesAround(s, x, y, z, dx));
					}
				}
			}
		};

		USE_AND_INIT_LOGGER("gcm.MpiStructuredSolver");
	};

	template<>
	class GridCharacteristicMethod<Cgal2DGrid<Elastic2DModel>> {
	public:
		void stage(const int s, const real &timeStep, Cgal2DGrid<Elastic2DModel>* mesh) {};
	};
}

#endif // LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
