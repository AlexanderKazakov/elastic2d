#ifndef LIBGCM_STRUCTUREDGRID_HPP
#define LIBGCM_STRUCTUREDGRID_HPP

#include <mpi.h>

#include "lib/grid/Grid.hpp"
#include "lib/nodes/Node.hpp"
#include "lib/interpolation/Interpolator.hpp"

namespace gcm {
	template<class TModel> class MpiStructuredSolver;
	template<class TModel> class VtkTextStructuredSnapshotter;

	template <class TModel>
	class StructuredGrid : public Grid<TModel> {

	public:
		typedef typename Grid<TModel>::Node Node;
		typedef typename Grid<TModel>::Vector Vector;
		typedef typename Grid<TModel>::Matrix Matrix;

	private:
		int accuracyOrder = 0; // order of accuracy of spatial interpolation

		int X = 0; // number of nodes along x direction of this grid (on this core)
		int Y = 0; // number of nodes along y direction
		int Z = 0; // number of nodes along z direction

		// the grid is divided equally along x-axis among processes
		int globalX = 0; // number of nodes along x direction of all grids (from all cores)
		int globalStartXindex = 0; // global x-index of the first real node of the grid

		real startX = 0.0; // global x-coordinate of the first real node of the grid
		real startY = 0.0; // global y-coordinate of the first real node of the grid
		real startZ = 0.0; // global z-coordinate of the first real node of the grid

		real h[3] = { 0.0,  /* x spatial step */
		              0.0,  /* y spatial step */
		              0.0,  /* z spatial step */ };

		/**
		 * Operator to iterate relatively to real nodes.
		 * Read / write access
		 * @param x x index < X
		 * @param y y index < Y
		 * @param z z index < Z
		 */
		inline Node &operator()(const int x, const int y, const int z) {
			return this->nodes[(2 * accuracyOrder + Z) * (2 * accuracyOrder + Y) * (x + accuracyOrder) + (2 * accuracyOrder + Z) * (y + accuracyOrder) + (z + accuracyOrder)];
		};

		/**
		 * Operator to iterate relatively real nodes.
		 * Read only access
		 * @param x x index < X
		 * @param y y index < Y
		 * @param z z index < Z
		 */
		inline const Node &get(const int x, const int y, const int z) const {
			return this->nodes[(2 * accuracyOrder + Z) * (2 * accuracyOrder + Y) * (x + accuracyOrder) + (2 * accuracyOrder + Z) * (y + accuracyOrder) + (z + accuracyOrder)];
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

		/* ---------------- For testing purposes ---------------- */

		const real &getH0ForTest() const { return h[0]; };

		const real &getH1ForTest() const { return h[1]; };

		const int getYForTest() const { return Y; };

		const int getXForTest() const { return X; };

		const int getStartXForTest() const { return globalStartXindex; };

		const Node &getNodeForTest(const int x, const int y, const int z) const { return get(x, y, z); };

		/* Change rheology in some area
		 *
		 * @param rho2rho0 = (rho in the area) / (default rho)
		 * @param lambda2lambda0 = (lambda in the area) / (default lambda)
		 * @param mu2mu0 = (mu in the area) / (default mu)
		 */
		virtual void changeRheology(const real &rho2rho0, const real &lambda2lambda0, const real &mu2mu0) override;

		/* ---------------- For testing purposes (end) ---------------- */

		virtual real getMinimalSpatialStep() const override;

	private:
		USE_AND_INIT_LOGGER("gcm.StructuredGrid");

		std::map<CUBIC_BORDERS, BorderCondition::CONDITION> borderConditions;

		virtual void applyBorderConditions() override;
		virtual void applyInitialConditions(const Task& task) override;
		linal::Vector3 getCoordinates(const int x, const int y, const int z) const {
			// TODO - replace with Vector3r
			return {startX + x * h[0], startY + y * h[1], startZ + z * h[2]};
		};

		friend class VtkTextStructuredSnapshotter<TModel>;
		friend class MpiStructuredSolver<TModel>;
	};
}

#endif // LIBGCM_STRUCTUREDGRID_HPP
