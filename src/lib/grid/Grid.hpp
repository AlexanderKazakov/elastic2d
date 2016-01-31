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
		void initialize(const Task &task) {
			rank = MPI::COMM_WORLD.Get_rank();
			numberOfWorkers = MPI::COMM_WORLD.Get_size();

			if (task.forceSequence) {
				rank = 0;
				numberOfWorkers = 1;
			}

			initializeImpl(task);
			applyInitialConditions(task);
		};
		virtual ~Grid() { };

		/** @warning this is not always just MPI-rank (see forceSequence) */
		int getRank() const {
			assert_ge(rank, 0);
			assert_lt(rank, numberOfWorkers);
			return rank;
		};
		/** @warning this is not always just MPI-numberOfWorkers (see forceSequence) */
		int getNumberOfWorkers() const {
			assert_gt(numberOfWorkers, 0);
			return numberOfWorkers;
		};

		/** @return maximal in modulus eigenvalue among all nodes all GcmMatrices of the mesh */
		real getMaximalLambda() const {
			assert_gt(maximalLambda, 0.0);
			return maximalLambda;
		};
		/** @return minimal spatial step for Courant condition */
		real getMinimalSpatialStep() const {
			assert_gt(minimalSpatialStep, 0.0);
			return minimalSpatialStep;
		};

		void beforeStage() { beforeStageImpl(); };
		void afterStage() { afterStageImpl(); };
		void beforeStep() {
			recalculateMinimalSpatialStep();
			recalculateMaximalLambda();
			beforeStepImpl();
		};
		void afterStep() { afterStepImpl(); };

	protected:
		int rank = -1; // index of core
		int numberOfWorkers = -1; // number of cores

		real maximalLambda = 0.0; // maximal in modulus eigenvalue among all nodes all GcmMatrices of the mesh
		real minimalSpatialStep = 0.0; // minimal spatial step over all mesh

		virtual void initializeImpl(const Task &task) = 0;
		virtual void beforeStageImpl() = 0;
		virtual void afterStageImpl() = 0;
		virtual void beforeStepImpl() = 0;
		virtual void afterStepImpl() = 0;
		virtual void recalculateMinimalSpatialStep() = 0;
		virtual void recalculateMaximalLambda() = 0;
		virtual void applyInitialConditions(const Task &task) = 0;
		virtual void applyBorderConditions() = 0;
	};
}

#endif // LIBGCM_GRID_HPP
