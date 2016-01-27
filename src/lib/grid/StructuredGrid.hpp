#ifndef LIBGCM_STRUCTUREDGRID_HPP
#define LIBGCM_STRUCTUREDGRID_HPP

#include <mpi.h>

#include <lib/grid/Grid.hpp>
#include <lib/grid/nodes/Node.hpp>
#include <lib/linal/special/VectorInt.hpp>
#include <lib/numeric/interpolation/Interpolator.hpp>
#include <lib/numeric/border_conditions/StructuredGridBorderConditions.hpp>

namespace gcm {
	template<class TNode> class DefaultSolver;
	template<class TNode> class GridCharacteristicMethod;
	template<class TNode> class VtkTextStructuredSnapshotter;
	template<class TNode> class Binary2DSeismograph;

	template<class TNode>
	class StructuredGrid : public Grid<TNode> {
	public:
		typedef typename Grid<TNode>::Node Node;
		typedef typename Grid<TNode>::Vector Vector;
		typedef typename Grid<TNode>::Matrix Matrix;

	protected:
		int accuracyOrder = 0; // order of accuracy of spatial interpolation

		linal::VectorInt<3> sizes = {0, 0, 0}; // numbers of nodes along each direction (on this core)

		// the grid is divided equally along x-axis among processes
		int globalStartXindex = 0; // global x-index of the first real node of the grid

		linal::Vector<3> startR = {0, 0, 0}; // global coordinates of the first real node of the grid
		linal::Vector<3> h = {0, 0, 0}; // spatial steps along each direction

		StructuredGridBorderConditions<TNode> borderConditions;
		
		/**
		 * Operator to iterate relatively to real nodes.
		 * Read / write access
		 * @param x x index < sizes(0)
		 * @param y y index < sizes(1)
		 * @param z z index < sizes(2)
		 */
		inline Node &operator()(const int x, const int y, const int z) {
			return this->nodes[ (unsigned long)
			                     (2 * accuracyOrder + sizes(2)) * (2 * accuracyOrder + sizes(1)) * (x + accuracyOrder)
			                   + (2 * accuracyOrder + sizes(2)) * (y + accuracyOrder)
			                   + (z + accuracyOrder) ];
		};

		/**
		 * Operator to iterate relatively real nodes.
		 * Read only access
		 * @param x x index < sizes(0)
		 * @param y y index < sizes(1)
		 * @param z z index < sizes(2)
		 */
		inline const Node &get(const int x, const int y, const int z) const {
			return this->nodes[ (unsigned long)
			                     (2 * accuracyOrder + sizes(2)) * (2 * accuracyOrder + sizes(1)) * (x + accuracyOrder)
			                   + (2 * accuracyOrder + sizes(2)) * (y + accuracyOrder)
			                   + (z + accuracyOrder) ];
		};

		/* Equal-distance spatial interpolator */
		Interpolator<Vector> interpolator;

		/**
		 * @param task properties and initial conditions etc
		 */
		virtual void initializeImpl(const Task &task) override;
		virtual void applyInitialConditions(const Task& task) override;

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

		virtual real getMinimalSpatialStepImpl() const override;

		void beforeStage();
		void afterStage();
		void exchangeNodesWithNeighbors();
		virtual void applyBorderConditions() override;
		linal::Vector3 getCoordinates(const int x, const int y, const int z) const {
			return startR + linal::plainMultiply(linal::VectorInt<3>({x, y, z}), h);
		};


		USE_AND_INIT_LOGGER("gcm.StructuredGrid");
		friend class VtkTextStructuredSnapshotter<TNode>;
		friend class GridCharacteristicMethod<TNode>;
		friend class DefaultSolver<TNode>;
		friend class StructuredGridBorderConditions<TNode>;
		friend class Binary2DSeismograph<TNode>;
	};
}

#endif // LIBGCM_STRUCTUREDGRID_HPP
