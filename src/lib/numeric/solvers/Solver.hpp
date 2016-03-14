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
		void initialize(const Task& task) {
			initializeImpl(task);
		}
		void beforeStatement(const Statement& statement) {
			currentTime = 0;
			beforeStatementImpl(statement);
		}
		void nextTimeStep() {
			real tau = calculateTau();
			nextTimeStepImpl();
			currentTime += tau;
		}
		void afterStatement() {
			afterStatementImpl();
		}
		virtual ~Solver() { }
		
		/** @return grid with actual values */
		virtual AbstractGrid* getGrid() const = 0;
		/** Calculate time step from Courant–Friedrichs–Lewy condition */
		virtual real calculateTau() const = 0;
		real getCurrentTime() const { return currentTime; }
		
	protected:
		virtual void initializeImpl(const Task& task) = 0;
		virtual void beforeStatementImpl(const Statement& statement) = 0;
		virtual void nextTimeStepImpl() = 0;
		virtual void afterStatementImpl() = 0;
		
	private:
		real currentTime = 0;
	};
}

#endif // LIBGCM_SOLVER_HPP
