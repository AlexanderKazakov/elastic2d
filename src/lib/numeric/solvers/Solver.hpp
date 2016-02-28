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
		virtual void initialize(const Task& task) {
			currentTime = 0.0; // solver can be used for several calculations
			initializeImpl(task);
		};
		virtual void initializeImpl(const Task &task) = 0;
		virtual void nextTimeStepImpl() = 0;
		virtual ~Solver() { };

		void nextTimeStep() {
			real tau = calculateTau();
			nextTimeStepImpl();
			currentTime += tau;
		};

		/** @return grid with actual values */
		virtual AbstractGrid* getGrid() const = 0;

		/** Calculate time step from Courant–Friedrichs–Lewy condition */
		virtual real calculateTau() const = 0;

		real getCurrentTime() const {
			return currentTime;
		};

	private:
		real currentTime = 0.0;
	};
}

#endif // LIBGCM_SOLVER_HPP
