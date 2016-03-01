#ifndef LIBGCM_ENGINE_HPP
#define LIBGCM_ENGINE_HPP

#include <lib/util/snapshot/Snapshotter.hpp>
#include <lib/numeric/solvers/Solver.hpp>
#include <lib/util/Logging.hpp>


namespace gcm {
	/**
	 * Responsible for the whole process of calculation
	 */
	class Engine {
	public:
		void initialize(const Task& task);

		/**
		 * Perform calculation of the task
		 */
		void run();


		void addSnapshotter(Snapshotter* snapshotter) { snapshotters.push_back(snapshotter); };
		void setSolver(Solver* solver_) { solver = solver_; };

	protected:
		Solver* solver = nullptr;
		std::vector<Snapshotter*> snapshotters;
		real requiredTime = 0;

		USE_AND_INIT_LOGGER("gcm.Engine");
	};
}

#endif // LIBGCM_ENGINE_HPP
