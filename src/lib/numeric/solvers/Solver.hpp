#ifndef LIBGCM_SOLVER_HPP
#define LIBGCM_SOLVER_HPP

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
	
	/**
	 * All necessary solver actions before a stage would be performed
	 */
	virtual void beforeStage() = 0;
	
	/**
	 * Calculate contact nodes in given stage
	 */
	virtual void contactStage(const int s, const real timeStep) = 0;
	
	/**
	 * Calculate inner and border nodes in given stage
	 */
	virtual void privateStage(const int s, const real timeStep) = 0;
	
	/**
	 * All necessary solver actions after stages performed
	 */
	virtual void afterStages(const real timeStep) = 0;
	
	
	/** @return grid with actual values */
	virtual AbstractGrid* getAbstractMesh() const = 0;
	
	/** Calculate time step from Courant–Friedrichs–Lewy condition */
	virtual real calculateTimeStep() const = 0;
	
};


}

#endif // LIBGCM_SOLVER_HPP
