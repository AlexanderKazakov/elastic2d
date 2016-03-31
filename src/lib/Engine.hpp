#ifndef LIBGCM_ENGINE_HPP
#define LIBGCM_ENGINE_HPP

#include <mpi.h>

#include <lib/util/task/Task.hpp>
#include <lib/util/Logging.hpp>

namespace gcm {
class Solver;
class Snapshotter;

/**
 * Current physical time in the calculated statement
 */
struct Clock {
	static real Time() {
		return time;
	}

	private:
		static real time;
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

protected:
	Solver* solver = nullptr;
	std::vector<Snapshotter*> snapshotters;
	Task task;
	real requiredTime = 0;

	/**
	 * Prepare to run statement
	 */
	void beforeStatement(const Statement& statement);

	/**
	 * Perform calculation of statement
	 */
	void runStatement();

	real estimateTimeStep();

	USE_AND_INIT_LOGGER("gcm.Engine")
};


}

#endif // LIBGCM_ENGINE_HPP
