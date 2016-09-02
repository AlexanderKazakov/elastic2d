#ifndef LIBGCM_ABSTRACTENGINE_HPP
#define LIBGCM_ABSTRACTENGINE_HPP

#include <libgcm/util/task/Task.hpp>
#include <libgcm/util/snapshot/Snapshotter.hpp>


/**
 * @mainpage
 *
 * GCM - grid-characteristic method
 * ================================
 * 
 * Library for numerical simulations of wave and associated processes
 * ------------------------------------------------------------------
 * 
 * ## About
 *
 * 1-2-3-dimensional simulation of wave and associated processes
 *
 */


namespace gcm {

/**
 * Main class. Responsible for the whole process of calculation
 */
class AbstractEngine {
public:
	typedef AbstractGrid::GridId GridId;
	
	
	AbstractEngine(const Task& task);
	virtual ~AbstractEngine() { }
	AbstractEngine(const AbstractEngine&) = delete;
	AbstractEngine& operator=(const AbstractEngine&) = delete;
	
	
	/**
	 * Perform all calculations
	 */
	void run();
	
	
protected:
	const real CourantNumber = 0;
	real requiredTime = 0;
	
	void afterConstruction(const Task& task);
	virtual void nextTimeStep() = 0;
	virtual real estimateTimeStep() = 0;
	virtual void writeSnapshots(const int step) = 0;
	
	
	USE_AND_INIT_LOGGER("gcm.AbstractEngine")
};


}

#endif // LIBGCM_ABSTRACTENGINE_HPP
