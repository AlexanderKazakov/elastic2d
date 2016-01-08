#ifndef LIBGCM_MPISOLVER_HPP
#define LIBGCM_MPISOLVER_HPP


#include "lib/snapshot/VtkTextStructuredSnapshotter.hpp"
#include "lib/grid/StructuredGrid.hpp"

namespace gcm {
	template<class TModel>
	class MpiStructuredSolver {
	public:

		void initialize(const Task& task, StructuredGrid<TModel> *mesh, StructuredGrid<TModel> *newMesh);

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

	private:
		real tau = 0.0; // time step
		real T = 0.0; // required time
		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition
		bool splittingSecondOrder = false;

		StructuredGrid<TModel>* mesh;
		StructuredGrid<TModel>* newMesh;

		Snapshotter* snapshotter = nullptr;

		void exchangeNodesWithNeighbors();

		/* ---------------- For testing purposes ---------------- */
	public:
		const real &getTauForTest() const { return tau; };

		const real &getTForTest() const { return T; };

		/* ---------------- For testing purposes (end) ---------------- */

	};
}

#endif //LIBGCM_MPISOLVER_HPP
