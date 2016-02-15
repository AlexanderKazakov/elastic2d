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
		 * @param grid grid to perform calculation
		 */
		void stage(const int s, const real &timeStep, TGrid* grid) {
			LOG_DEBUG("Start stage " << s << " timeStep = " << timeStep);

			for (auto it : *grid) {
				// points to interpolate values on previous time layer
				auto dx = -timeStep * linal::diag(grid->matrix(it)->A(s).L);

				grid->_pdeNew(it) =
						/* new values = U1 * Riemann solvers */
						grid->matrix(it)->A(s).U1 *
						/* Riemann solvers = U * old values */
						grid->matrix(it)->A(s).U.diagonalMultiply
								/* old values are in columns of the matrix */
								(grid->interpolateValuesAround(s, it, dx));
			}
		};

		USE_AND_INIT_LOGGER("gcm.MpiStructuredSolver");
	};

}

#endif // LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
