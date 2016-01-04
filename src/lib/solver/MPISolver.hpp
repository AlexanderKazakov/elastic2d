#ifndef LIBGCM_MPISOLVER_HPP
#define LIBGCM_MPISOLVER_HPP


#include "lib/grid/StructuredGrid.hpp"

namespace gcm {
	template<class TModel>
	class MPISolver {
	public:
		MPISolver(StructuredGrid<TModel> *mesh, StructuredGrid<TModel> *newMesh);

		/**
		 * Perform calculation of the task
		 */
		void calculate();

		/**
		 * Do next stage of splitting method
		 * @param s direction
		 * @param timeStep time step
		 */
		void stage(const int s, const real &timeStep);

		bool makeSnapshots = false;
		bool splittingSecondOrder = false;

	private:

		StructuredGrid<TModel> *mesh;
		StructuredGrid<TModel> *newMesh;

		void exchangeNodesWithNeighbors();

	};
}

#endif //LIBGCM_MPISOLVER_HPP
