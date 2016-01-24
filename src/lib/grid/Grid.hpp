#ifndef LIBGCM_GRID_HPP
#define LIBGCM_GRID_HPP

#include <vector>
#include <mpi.h>

#include <lib/util/Logging.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/linal/Linal.hpp>

namespace gcm {
	template<class TModel> class MpiStructuredSolver;

	template <class TNode>
	class Grid {
	public:
		typedef TNode Node;
		typedef typename Node::Vector Vector;
		typedef typename Node::Matrix Matrix;

		/** @param task properties and initial conditions etc */
		void initialize(const Task &task);

		int getRank() const {
			assert_ge(rank, 0);
			assert_lt(rank, numberOfWorkers);
			return rank;
		};
		int getNumberOfWorkers() const {
			assert_gt(numberOfWorkers, 0);
			return numberOfWorkers;
		};
		real getMaximalLambda() const {
			assert_gt(maximalLambda, 0.0);
			return maximalLambda;
		};
		real getMinimalSpatialStep() const {
			real minH = this->getMinimalSpatialStepImpl();
			assert_gt(minH, 0.0);
			return minH;
		};

		/**
		 * Change rheology in some area
		 *
		 * @param rho2rho0 = (rho in the area) / (default rho)
		 * @param lambda2lambda0 = (lambda in the area) / (default lambda)
		 * @param mu2mu0 = (mu in the area) / (default mu)
		 */
		virtual void changeRheology(const real &rho2rho0, const real &lambda2lambda0, const real &mu2mu0) = 0;

	protected:
		int rank = -1; // index of core
		int numberOfWorkers = -1; // number of cores

		/* Node storage */
		std::vector<Node> nodes;

		real maximalLambda = 0.0; // maximal in modulus eigenvalue among all nodes all GcmMatrices of the mesh

		virtual void initializeImpl(const Task &task) = 0;
		virtual void applyInitialConditions(const Task &task) = 0;
		virtual real getMinimalSpatialStepImpl() const = 0;
		virtual void applyBorderConditions() = 0;
		
		friend class MpiStructuredSolver<TNode>;
	};
}

#endif // LIBGCM_GRID_HPP
