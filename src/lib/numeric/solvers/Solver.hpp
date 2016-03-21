#ifndef LIBGCM_SOLVER_HPP
#define LIBGCM_SOLVER_HPP

#include <lib/numeric/gcm/GridCharacteristicMethod.hpp>
#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>

namespace gcm {
	class AbstractGrid;

	/**
	 * Class for handling complete time step
	 */
	class Solver {
	public:
		Solver(const Task&) { }
		virtual ~Solver() { }

		virtual void beforeStatement(const Statement& statement) = 0;
		virtual void afterStatement() = 0;

		virtual void nextTimeStep(const real timeStep) = 0;

		/** @return grid with actual values */
		virtual AbstractGrid* getActualGrid() const = 0;

		/** Calculate time step from Courant–Friedrichs–Lewy condition */
		virtual real calculateTimeStep() const = 0;
	};
}

#endif // LIBGCM_SOLVER_HPP
