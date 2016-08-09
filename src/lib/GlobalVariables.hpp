#ifndef LIBGCM_GLOBALVARIABLES_HPP
#define LIBGCM_GLOBALVARIABLES_HPP

#include <mpi.h>
#include <lib/util/Types.hpp>


namespace gcm {


/**
 * Current physical time and time step used by the program.
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
	
	/** Physical time in the calculated statement at next time layer */
	static real TimeAtNextTimeLayer() {
		return Time() + TimeStep();
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
		friend class Engine;
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
	 * calculation of every concrete statement is performed in sequence.
	 * Useful for inverse problem calculation and MPI testing.
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
		friend class Engine;
};


}

#endif // LIBGCM_GLOBALVARIABLES_HPP