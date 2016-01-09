#ifndef LIBGCM_STRUCTUREDGRID_HPP
#define LIBGCM_STRUCTUREDGRID_HPP

#include <memory>
#include <mpi.h>

#include "lib/util/Logging.hpp"
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
		int globalX = 0; // number of nodes along x direction of all meshes (from all cores)
		int startX = 0; // global x-index of first real node of the grid

		real h[3] = {0.0,  /* x spatial step */
		             0.0,  /* y spatial step */
		             0.0,  /* z spatial step */ };

		/* ------------------ Properties and conditions (end) ------------------ */


		/* GcmMatrices that is common for majority of nodes */
		std::shared_ptr<typename TModel::GcmMatrices> defaultMatrix;
		real maximalLambda = 0.0; // maximal eigenvalue among all nodes all GcmMatrices of the mesh

		/* Data storage. Real nodes with assistant nodes on borders */
		Node *nodes = nullptr;

		/**
		 * Operator to iterate relatively to real nodes.
		 * Read / write access
		 * @param x x index < X
		 * @param y y index < Y
		 * @param z z index < Z
		 */
		inline Node &operator()(const int x, const int y, const int z) {
			return nodes[(2 * accuracyOrder + Z) * (2 * accuracyOrder + Y) * (x + accuracyOrder) + (2 * accuracyOrder + Z) * (y + accuracyOrder) + (z + accuracyOrder)];
		};

		/**
		 * Operator to iterate relatively real nodes.
		 * Read only access
		 * @param x x index < X
		 * @param y y index < Y
		 * @param z z index < Z
		 */
		inline const Node &get(const int x, const int y, const int z) const {
			return nodes[(2 * accuracyOrder + Z) * (2 * accuracyOrder + Y) * (x + accuracyOrder) + (2 * accuracyOrder + Z) * (y + accuracyOrder) + (z + accuracyOrder)];
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
		 * @param x x-index of the reference node
		 * @param y y-index of the reference node
		 * @param z z-index of the reference node
		 * @param dx Vector of distances from reference node on which
		 * values should be interpolated
		 * @return Matrix with interpolated nodal values in columns
		 */
		Matrix interpolateValuesAround
				(const int stage, const int x, const int y, const int z, const Vector& dx) const;

		/**
		 * Place in src nodal values which are necessary for
		 * interpolation in specified point. The number of placed
		 * in values is equal to src.size()
		 */
		void findSourcesForInterpolation(const int stage, const int x, const int y, const int z,
		                                 const real &dx, std::vector<Vector>& src) const;

		friend class VtkTextStructuredSnapshotter<TModel>;
		friend class MpiStructuredSolver<TModel>;

		/* ---------------- For testing purposes ---------------- */
	public:
		const real &getH0ForTest() const { return h[0]; };

		const real &getH1ForTest() const { return h[1]; };

		const int getYForTest() const { return Y; };

		const int getXForTest() const { return X; };

		const int getStartXForTest() const { return startX; };

		const Node &getNodeForTest(const int x, const int y, const int z) const { return get(x, y, z); };

		/* Change rheology in some area
		 *
		 * @param rho2rho0 = (rho in the area) / (default rho)
		 * @param lambda2lambda0 = (lambda in the area) / (default lambda)
		 * @param mu2mu0 = (mu in the area) / (default mu)
		 */
		virtual void changeRheology(const real &rho2rho0, const real &lambda2lambda0, const real &mu2mu0) override;

		/* ---------------- For testing purposes (end) ---------------- */

	private:
		USE_AND_INIT_LOGGER("gcm.StructuredGrid");
		virtual void applyBorderConditions() override;
		virtual void applyInitialConditions() override;
	public:
		virtual real getMinimalSpatialStep() const override;
		real getMaximalLambda() const { return maximalLambda; };
	};
}

#endif // LIBGCM_STRUCTUREDGRID_HPP
