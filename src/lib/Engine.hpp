#ifndef LIBGCM_ENGINE_HPP
#define LIBGCM_ENGINE_HPP

#include <lib/util/snapshot/VtkTextStructuredSnapshotter.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/util/Logging.hpp>


namespace gcm {
	/**
	 * Responsible for the whole process of calculation
	 */
	template<class TNode>
	class Engine {
	public:
		void initialize(const Task& task);
		~Engine();

		/**
		 * Perform calculation of the task
		 */
		void run();

	protected:
		DefaultSolver<TNode>* solver = nullptr;
		Snapshotter<TNode>* snapshotter = nullptr;
		real requiredTime = 0;

		USE_AND_INIT_LOGGER("gcm.Engine");
	};
}

#endif // LIBGCM_ENGINE_HPP
