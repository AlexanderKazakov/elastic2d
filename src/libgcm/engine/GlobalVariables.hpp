#ifndef LIBGCM_GLOBALVARIABLES_HPP
#define LIBGCM_GLOBALVARIABLES_HPP

#include <mpi.h>
#include <libgcm/util/infrastructure/infrastructure.hpp>


namespace gcm {

class AbstractEngine;

/**
 * Current physical time and time step used by the program.
 * TODO - replace it to the AbstractEngine, don't use global vars
 */
struct Clock {
	/** Current physical time in the calculated statement */
	static real Time() {
		return time;
	}
	
	/** Time step (aka "tau") used by the whole program to calculate the time layer */
	static real TimeStep() {
		return timeStep;
	}
	
private:
	static real time;
	static real timeStep;
	
	static void setZero() {
		time = timeStep = 0;
	}
	
	static void tickTack() {
		time += timeStep;
	}
	
	/// DO NOT put another "friends" here!
	friend class AbstractEngine;
};


/**
 * MPI information
 */
struct Mpi {
	static int Rank() {
		return rank;
	}

	static int Size() {
		return size;
	}

	/**
	 * If true, independently from number of working processes
	 * calculation is performed in sequence. Useful for MPI testing.
	 */
	static bool ForceSequence() {
		return forceSequence;
	}
	
private:
	static int rank;
	static int size;
	static bool forceSequence;

	static void initialize(const bool forceSequence_) {
		rank = MPI::COMM_WORLD.Get_rank();
		size = MPI::COMM_WORLD.Get_size();
		forceSequence = forceSequence_;
	}
	
	/// DO NOT put another "friends" here!
	friend class AbstractEngine;
};


}

#endif // LIBGCM_GLOBALVARIABLES_HPP
