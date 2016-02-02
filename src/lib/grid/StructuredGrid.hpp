#ifndef LIBGCM_STRUCTUREDGRID_HPP
#define LIBGCM_STRUCTUREDGRID_HPP

#include <mpi.h>

#include <lib/grid/Grid.hpp>
#include <lib/grid/nodes/Nodes.hpp>
#include <lib/linal/linal.hpp>
#include <lib/numeric/interpolation/Interpolator.hpp>
#include <lib/numeric/border_conditions/StructuredGridBorderConditions.hpp>

namespace gcm {
	template<class TGrid> class DefaultSolver;
	template<class TGrid> class GridCharacteristicMethod;
	template<class TGrid> class VtkTextStructuredSnapshotter;
	template<class TGrid> class Binary2DSeismograph;

	template<typename TModel, template<typename> class GcmMatricesStorage> class Node;

	/**
	 * Non-movable structured rectangular grid
	 * @tparam TNode building block of mesh
	 * @tparam TModel rheology model
	 */
	template<typename TModel,
			template<typename> class GcmMatricesStorage = DefaultGcmMatricesStorage>
	class StructuredGrid : public Grid {
	public:
		typedef TModel Model;
		typedef Node<TModel, GcmMatricesStorage> NODE;
		typedef typename NODE::Vector Vector;
		typedef typename NODE::GCM_MATRICES GCM_MATRICES;
		typedef typename GCM_MATRICES::Matrix Matrix;
		static const int DIMENSIONALITY = TModel::DIMENSIONALITY;

	protected:
		/* Node storage */
		std::vector<NODE> nodes;

		int accuracyOrder = 0; // order of accuracy of spatial interpolation

		linal::VectorInt<3> sizes = {0, 0, 0}; // numbers of nodes along each direction (on this core)

		// the grid is divided equally along x-axis among processes
		int globalStartXindex = 0; // global x-index of the first real node of the grid

		linal::Vector<3> startR = {0, 0, 0}; // global coordinates of the first real node of the grid
		linal::Vector<3> h = {0, 0, 0}; // spatial steps along each direction

		StructuredGridBorderConditions<StructuredGrid> borderConditions;
		
		/**
		 * Operator to iterate relatively to real nodes.
		 * Read / write access
		 * @param x x index < sizes(0)
		 * @param y y index < sizes(1)
		 * @param z z index < sizes(2)
		 */
		inline NODE &operator()(const int x, const int y, const int z) {
			return this->nodes[ (unsigned long)
			                     (2 * accuracyOrder + sizes(2)) * (2 * accuracyOrder + sizes(1)) * (x + accuracyOrder)
			                   + (2 * accuracyOrder + sizes(2)) * (y + accuracyOrder)
			                   + (z + accuracyOrder) ];
		};

	public:
		/**
		 * Operator to iterate relatively real nodes.
		 * Read only access
		 * @param x x index < sizes(0)
		 * @param y y index < sizes(1)
		 * @param z z index < sizes(2)
		 */
		inline const NODE &get(const int x, const int y, const int z) const {
			return this->nodes[ (unsigned long)
			                     (2 * accuracyOrder + sizes(2)) * (2 * accuracyOrder + sizes(1)) * (x + accuracyOrder)
			                   + (2 * accuracyOrder + sizes(2)) * (y + accuracyOrder)
			                   + (z + accuracyOrder) ];
		};

	protected:
		/* Equal-distance spatial interpolator */
		Interpolator<Vector> interpolator;

		/**
		 * Interpolate nodal values in specified points.
		 * Interpolated value for k-th point in vector dx are
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

		virtual void initializeImpl(const Task &task) override;
		virtual void applyInitialConditions(const Task& task) override;
		virtual void beforeStageImpl() override;
		virtual void afterStageImpl() override { };
		virtual void beforeStepImpl() override { };
		virtual void afterStepImpl() override { };
		virtual void recalculateMinimalSpatialStep() override { };
		virtual void recalculateMaximalLambda() override { /* TODO for non-linear materials */ };

		void exchangeNodesWithNeighbors();
		virtual void applyBorderConditions() override;
		linal::Vector3 getCoordinates(const int x, const int y, const int z) const {
			return startR + linal::plainMultiply(linal::VectorInt<3>({x, y, z}), h);
		};

		USE_AND_INIT_LOGGER("gcm.StructuredGrid");
		friend class GridCharacteristicMethod<StructuredGrid>;
		friend class DefaultSolver<StructuredGrid>;
		friend class StructuredGridBorderConditions<StructuredGrid>;
		// it's better than create million gets
		friend class VtkTextStructuredSnapshotter<StructuredGrid>;
		friend class Binary2DSeismograph<StructuredGrid>;
	};
}

#endif // LIBGCM_STRUCTUREDGRID_HPP
