#ifndef LIBGCM_STRUCTUREDGRID_HPP
#define LIBGCM_STRUCTUREDGRID_HPP

#include <memory>
#include <mpi.h>

#include "lib/grid/Grid.hpp"
#include "lib/linal/Vector.hpp"
#include "lib/nodes/Node.hpp"
#include "lib/interpolation/Interpolator.hpp"
#include "lib/Task.hpp"

namespace gcm {
	template<class TModel> class MpiStructuredSolver;
	template<class TModel> class VtkTextStructuredSnapshotter;

	template <class TModel>
	class StructuredGrid : public Grid {
	public:
		typedef typename TModel::Node Node;
		typedef typename Node::Vector Vector;
		typedef typename Node::Matrix Matrix;
		/* ------------------ Properties and conditions ------------------ */

		int accuracyOrder = 0; // order of accuracy of spatial interpolation

		int X = 0; // number of nodes along x direction
		int Y = 0; // number of nodes along y direction of this grid (on this core)
		int Z = 0; // number of nodes along z direction
		int globalY = 0; // number of nodes along y direction of all meshes (all cores)

		real h[3] = {0.0,  /* x spatial step */
		             0.0,  /* y spatial step */
		             0.0,  /* z spatial step */ };

		/* ------------------ Properties and conditions (end) ------------------ */

		int startY = 0; // global y-index of first real node of the grid

		/* GcmMatrices that is common for majority of nodes */
		std::shared_ptr<typename TModel::GcmMatrices> defaultMatrix;
		real maximalLambda = 0.0; // maximal eigenvalue among all nodes all GcmMatrices of the mesh

		/* Data storage. Real nodes with assistant nodes on borders */
		Node *nodes = nullptr;

		/**
		 * Operator to iterate relatively to real nodes.
		 * Read / write access
		 * @param y y index < Y
		 * @param x x index < X
		 */
		inline Node &operator()(const int y, const int x) {
			return nodes[(2 * accuracyOrder + X) * (y + accuracyOrder) + (x + accuracyOrder)];
		};

		/**
		 * Operator to iterate relatively real nodes.
		 * Read only access
		 * @param y y index < Y
		 * @param x x index < X
		 */
		inline const Node &get(const int y, const int x) const {
			return nodes[(2 * accuracyOrder + X) * (y + accuracyOrder) + (x + accuracyOrder)];
		};

		/* Equal-distance spatial interpolator */
		Interpolator<Vector> interpolator;

	public:
		/**
		 * @param task properties and initial conditions etc
		 */
		virtual void initializeImpl(const Task &task) override;

		/**
		 * Interpolate nodal values in specified points.
		 * Interpolated value for k-th point in vector %dx are
		 * stored in k-th column of returned Matrix.
		 * @param stage direction
		 * @param y y-index of the reference node
		 * @param x x-index of the reference node
		 * @param dx Vector of distances from reference node on which
		 * values should be interpolated
		 * @return Matrix with interpolated nodal values in columns
		 */
		Matrix interpolateValuesAround
				(const int stage, const int y, const int x, const Vector& dx) const;

		/**
		 * Place in src nodal values which are necessary for
		 * interpolation in specified point. The number of placed
		 * in values is equal to src.size()
		 */
		void findSourcesForInterpolation(const int stage, const int y, const int x,
		                                 const real &dx, std::vector<Vector>& src) const;

		friend class VtkTextStructuredSnapshotter<TModel>;
		friend class MpiStructuredSolver<TModel>;
		// TODO - replace after implementing iterator
		friend class Grid;
		// TODO (end)

		/* ---------------- For testing purposes ---------------- */
	public:
		const real &getH0ForTest() const { return h[0]; };

		const real &getH1ForTest() const { return h[1]; };

		const int getYForTest() const { return Y; };

		const int getXForTest() const { return X; };

		const int getStartYForTest() const { return startY; };

		const Node &getNodeForTest(const int y, const int x) const { return get(y, x); };

		/* Change rheology in some area
		 *
		 * @param rho2rho0 = (rho in the area) / (default rho)
		 * @param lambda2lambda0 = (lambda in the area) / (default lambda)
		 * @param mu2mu0 = (mu in the area) / (default mu)
		 */
		virtual void changeRheology(const real &rho2rho0, const real &lambda2lambda0, const real &mu2mu0) override;

		/* ---------------- For testing purposes (end) ---------------- */

	private:
		virtual void applyBorderConditions() override;
		virtual void applyInitialConditions() override;
	public:
		virtual real getMinimalSpatialStep() const override;
		real getMaximalLambda() const { return maximalLambda; };
	};
}

#endif // LIBGCM_STRUCTUREDGRID_HPP
