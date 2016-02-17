#ifndef LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHOD_HPP

#include <lib/linal/linal.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/Types.hpp>
#include <lib/numeric/gcmethod/GcmHandler.hpp>


namespace gcm {
	template<class TGrid>
	class GridCharacteristicMethod {
		typedef typename TGrid::GCM_HANDLER GCM_HANDLER;
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
								(GCM_HANDLER::interpolateValuesAround(*grid, s, it, dx));
			}
		};

		USE_AND_INIT_LOGGER("gcm.GridCharacteristicMethod");
	};
}

#endif // LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
