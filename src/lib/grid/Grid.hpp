#ifndef LIBGCM_GRID_HPP
#define LIBGCM_GRID_HPP

#include <mpi.h>

#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>

namespace gcm {
	/**
	 * Base class for all grids
	 */
	class Grid {
	public:
		/** @param task properties and initial conditions etc */
		void initialize(const Task &task);
		virtual ~Grid() { };

		int getRank() const {
			assert_ge(rank, 0);
			assert_lt(rank, numberOfWorkers);
			return rank;
		};
		int getNumberOfWorkers() const {
			assert_gt(numberOfWorkers, 0);
			return numberOfWorkers;
		};
		real getMaximalLambda() const {
			assert_gt(maximalLambda, 0.0);
			return maximalLambda;
		};
		real getMinimalSpatialStep() const {
			real minH = this->getMinimalSpatialStepImpl();
			assert_gt(minH, 0.0);
			return minH;
		};

	protected:
		int rank = -1; // index of core
		int numberOfWorkers = -1; // number of cores

		real maximalLambda = 0.0; // maximal in modulus eigenvalue among all nodes all GcmMatrices of the mesh

		virtual void initializeImpl(const Task &task) = 0;
		virtual void applyInitialConditions(const Task &task) = 0;
		virtual real getMinimalSpatialStepImpl() const = 0;
		virtual void applyBorderConditions() = 0;
	};
}

#endif // LIBGCM_GRID_HPP
