#ifndef LIBGCM_ENGINE_HPP
#define LIBGCM_ENGINE_HPP

#include <mpi.h>

#include <lib/util/task/Task.hpp>
#include <lib/mesh/grid/AbstractGlobalScene.hpp>
#include <lib/util/Logging.hpp>

/**
 * @mainpage
 *
 * GCM - grid-characteristic method
 * ================================
 * 
 * Engine for numerical simulations of wave and associated processes
 * -----------------------------------------------------------------
 * 
 * ## About
 *
 * 1-2-3-dimensional simulation of wave and associated processes
 *
 */


namespace gcm {
class Solver;
class Snapshotter;


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


/**
 * Main class. Responsible for the whole process of calculation
 */
class Engine {
public:
	Engine(const Task& task_);
	~Engine();
	Engine(const Engine&) = delete;
	Engine& operator=(const Engine&) = delete;
	
	
	/**
	 * Perform calculation of the whole task (it can be several statements)
	 */
	void run();
	
	
	/**
	 * Prepare to run statement
	 */
	void beforeStatement(const Statement& statement);
	
	
	/**
	 * Perform calculation of statement after preparation
	 */
	void runStatement();
	
	
	/// for tests
	const Solver* getSolver() const {
		assert_eq(1, bodies.size());
		return bodies.front().solver;
	}
	
private:
	
	struct Body {
		Solver* solver;
		std::vector<Snapshotter*> snapshotters;
	};
	
	
	AbstractGlobalScene* globalScene;
	std::vector<Body> bodies;
	Task task;
	real requiredTime = 0;
	
	
	void estimateTimeStep();
	
	
	USE_AND_INIT_LOGGER("gcm.Engine")
};


}

#endif // LIBGCM_ENGINE_HPP
