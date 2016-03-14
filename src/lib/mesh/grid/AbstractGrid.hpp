#ifndef LIBGCM_ABSTRACTGRID_HPP
#define LIBGCM_ABSTRACTGRID_HPP

#include <mpi.h>

#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>

namespace gcm {
	/**
	 * Base class for all grids
	 */
	class AbstractGrid {
	public:
		/** @param task properties and initial conditions etc */
		AbstractGrid(const Task &task) {
			// TODO - should rank be here?
			rank = MPI::COMM_WORLD.Get_rank();
			numberOfWorkers = MPI::COMM_WORLD.Get_size();
			if (task.forceSequence) {
				rank = 0;
				numberOfWorkers = 1;
			}
		}
		virtual ~AbstractGrid() { }

		void beforeStatement(const Statement& statement) {
			beforeStatementImpl(statement);
		}

		/** @warning this is not always just MPI-rank (see forceSequence) */
		int getRank() const {
			assert_ge(rank, 0);
			assert_lt(rank, numberOfWorkers);
			return rank;
		}
		/** @warning this is not always just MPI-size (see forceSequence) */
		int getNumberOfWorkers() const {
			assert_gt(numberOfWorkers, 0);
			return numberOfWorkers;
		}

		/** @return maximal in modulus eigenvalue among all nodes all GcmMatrices of the mesh */
		real getMaximalLambda() const {
			assert_gt(maximalLambda, 0.0);
			return maximalLambda;
		}
		/** @return minimal spatial step for Courant condition */
		real getMinimalSpatialStep() const {
			assert_gt(minimalSpatialStep, 0.0);
			return minimalSpatialStep;
		}

		std::string getId() const { return id; }

		void beforeStep() {
			recalculateMinimalSpatialStep();
			recalculateMaximalLambda();
		}
		void afterStep() { }

	protected:
		std::string id = "mesh"; // name of the mesh
		int rank = -1; // index of core
		int numberOfWorkers = -1; // number of cores

		real maximalLambda = 0.0; // maximal in modulus eigenvalue among all nodes all GcmMatrices of the mesh
		real minimalSpatialStep = 0.0; // minimal spatial step over all mesh

		virtual void beforeStatementImpl(const Statement& statement) = 0;
		virtual void afterStatement() = 0;

		virtual void recalculateMinimalSpatialStep() = 0;
		virtual void recalculateMaximalLambda() = 0;
	};
}

#endif // LIBGCM_ABSTRACTGRID_HPP
