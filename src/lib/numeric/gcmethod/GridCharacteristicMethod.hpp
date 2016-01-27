#ifndef LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
#define LIBGCM_GRIDCHARACTERISTICMETHOD_HPP

#include <lib/grid/StructuredGrid.hpp>
#include <lib/util/Logging.hpp>

namespace gcm {
	template<class TNode>
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
		           const StructuredGrid<TNode>* mesh, StructuredGrid<TNode>* newMesh);

		USE_AND_INIT_LOGGER("gcm.MpiStructuredSolver");

	};
}

#endif // LIBGCM_GRIDCHARACTERISTICMETHOD_HPP
