#ifndef LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHOD_HPP

#include <lib/util/Logging.hpp>
#include <lib/util/Types.hpp>

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
		void stage(const int s, const real &timeStep,
		           const TGrid* mesh, TGrid* newMesh);

		USE_AND_INIT_LOGGER("gcm.MpiStructuredSolver");
	};
}

#endif // LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
