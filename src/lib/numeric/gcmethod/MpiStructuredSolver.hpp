#ifndef LIBGCM_MPISOLVER_HPP
#define LIBGCM_MPISOLVER_HPP


#include <lib/util/snapshot/VtkTextStructuredSnapshotter.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/util/Logging.hpp>

namespace gcm {
	template<class TNode>
	class MpiStructuredSolver {
	public:
		void initialize(const Task& task);
		~MpiStructuredSolver();

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

		StructuredGrid<TNode>* mesh = nullptr;
		StructuredGrid<TNode>* newMesh = nullptr;

		Snapshotter<TNode>* snapshotter = nullptr;

		USE_AND_INIT_LOGGER("gcm.MpiStructuredSolver");

		void exchangeNodesWithNeighbors();

		/* ---------------- For testing purposes ---------------- */
	public:
		const real &getTauForTest() const { return tau; };

		const real &getTForTest() const { return T; };

		StructuredGrid<TNode> *getMesh() const { return mesh; };

		StructuredGrid<TNode> *getNewMesh() const { return newMesh; };

		/* ---------------- For testing purposes (end) ---------------- */

	};
}

#endif //LIBGCM_MPISOLVER_HPP
