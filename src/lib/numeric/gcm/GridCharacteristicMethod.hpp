#ifndef LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHOD_HPP

#include <lib/linal/linal.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/Types.hpp>
#include <lib/numeric/gcm/GcmHandler.hpp>


namespace gcm {
template<class TMesh>
class GridCharacteristicMethod {
public:
	typedef typename TMesh::GCM_HANDLER GCM_HANDLER;
	
	/**
	 * Do grid-characteristic stage of splitting method
	 * @param s direction aka stage
	 * @param timeStep time step
	 * @param mesh mesh to perform calculation
	 */
	void stage(const int s, const real& timeStep, TMesh& mesh) {
		gcmHandler.beforeStage(mesh);

		for (auto it : mesh) {
			// points to interpolate values on previous time layer
			auto dx = -timeStep * linal::diag(mesh.matrices(it)->m[s].L);

			mesh._pdeNew(it) =
				/* new values = U1 * Riemann solvers */
				mesh.matrices(it)->m[s].U1 *
					linal::diagonalMultiply(
						/* Riemann solvers = U * old values */
						mesh.matrices(it)->m[s].U,
						/* old values are in columns of the matrix */
						gcmHandler.interpolateValuesAround(mesh, s, it, dx));
			
			gcmHandler.borderCorrector(mesh, s, it);
		}
	}

private:
	GCM_HANDLER gcmHandler; ///< helper

	USE_AND_INIT_LOGGER("gcm.GridCharacteristicMethod")
};


}

#endif // LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
