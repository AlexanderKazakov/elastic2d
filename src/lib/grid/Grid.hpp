#ifndef LIBGCM_GRID_HPP
#define LIBGCM_GRID_HPP

#include <memory>
#include <mpi.h>

#include "lib/nodes/Node.hpp"
#include "lib/Task.hpp"

namespace gcm {
	class Grid {
	protected:
		int rank = 0; // index of core
		int numberOfWorkers = 0; // number of cores

		/* ------------------ Properties and conditions ------------------ */

		InitialConditions initialConditions = InitialConditions::Zero;
		std::map <std::string, BorderConditions> borderConditions;

		/* ------------------ Properties and conditions (end) ------------------ */

		virtual void initializeImpl(const Task &task) = 0;

	public:
		/** StructuredGrid factory
		 * @param task properties and initial conditions etc
		 */
		void initialize(const Task &task);

	private:
		virtual void applyBorderConditions() = 0;
		virtual void applyInitialConditions() = 0;
		virtual real getMinimalSpatialStep() const = 0;

	public:
		int getRank() const { return rank; };
		int getNumberOfWorkers() const { return numberOfWorkers; };

		/**
		 * Change rheology in some area
		 *
		 * @param rho2rho0 = (rho in the area) / (default rho)
		 * @param lambda2lambda0 = (lambda in the area) / (default lambda)
		 * @param mu2mu0 = (mu in the area) / (default mu)
		 */
		virtual void changeRheology(const real &rho2rho0, const real &lambda2lambda0, const real &mu2mu0) = 0;
	};
}

#endif // LIBGCM_GRID_HPP
