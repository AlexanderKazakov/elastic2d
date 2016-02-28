#ifndef LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHOD_HPP

#include <lib/linal/linal.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/Types.hpp>
#include <lib/numeric/gcm/GcmHandler.hpp>


namespace gcm {
	template<class TMesh>
	class GridCharacteristicMethod {
		typedef typename TMesh::GCM_HANDLER GCM_HANDLER;
	public:
		/**
		 * Do grid-characteristic stage of splitting method
		 * @param s direction
		 * @param timeStep time step
		 * @param mesh mesh to perform calculation
		 */
		void stage(const int s, const real &timeStep, TMesh* mesh) {
			LOG_DEBUG("Start stage " << s << " timeStep = " << timeStep);

			for (auto it : *mesh) {
				// points to interpolate values on previous time layer
				auto dx = -timeStep * linal::diag(mesh->matrix(it)->A(s).L);
				
				mesh->_pdeNew(it) =
						/* new values = U1 * Riemann solvers */
						mesh->matrix(it)->A(s).U1 *
						/* Riemann solvers = U * old values */
						mesh->matrix(it)->A(s).U.diagonalMultiply
								/* old values are in columns of the matrix */
								(GCM_HANDLER::interpolateValuesAround(*mesh, s, it, dx));
			}
		};

		USE_AND_INIT_LOGGER("gcm.GridCharacteristicMethod");
	};
}

#endif // LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
