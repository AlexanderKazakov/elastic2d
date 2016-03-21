#ifndef LIBGCM_ENGINE_HPP
#define LIBGCM_ENGINE_HPP

#include <mpi.h>

#include <lib/util/task/Task.hpp>
#include <lib/util/Logging.hpp>

namespace gcm {
	class Solver;
	class Snapshotter;
	
	/**
	 * Main class. Responsible for the whole process of calculation
	 * from time = 0 to time = requiredTime.
	 * Clock of the program, singleton.
	 */
	class Engine {
	public:
		static Engine& Instance() {
			static Engine instance;
			return instance;
		}
	private:
		Engine() : MpiRank(MPI::COMM_WORLD.Get_rank()),
		           MpiSize(MPI::COMM_WORLD.Get_size()) { }
		~Engine() { }
	public:
		Engine(const Engine&) = delete;
		Engine& operator=(const Engine&) = delete;
		
		const int MpiRank;
		const int MpiSize;

		void initialize(const Task& task_);

		/**
		 * Perform calculation of the whole task (it can be several statements)
		 */
		void run();
		
		real getCurrentTime() const { return currentTime; }
		bool getForceSequence() const { return forceSequence; }

	protected:
		Solver* solver = nullptr;
		std::vector<Snapshotter*> snapshotters;

		/** 
		 * If true, independently from number of working processes
		 * calculation of a concrete statement is performed in sequence.
		 * Useful for inverse problem calculation and MPI testing.
		 */
		bool forceSequence = false;
		/**
		 * Required physical time to calculate the statement
		 */
		real requiredTime = 0;
		/**
		 * Current physical time in the calculated statement
		 */
		real currentTime = 0;

		Task task;

		/**
		 * Perform calculation of statement
         */
		void runStatement();
		
		void beforeStatement(const Statement& statement);

		real estimateTimeStep();
		
		USE_AND_INIT_LOGGER("gcm.Engine")
	};
}

#endif // LIBGCM_ENGINE_HPP
