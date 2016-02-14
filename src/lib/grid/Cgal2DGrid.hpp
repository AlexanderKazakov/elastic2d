#ifndef LIBGCM_CGAL2DGRID_HPP
#define LIBGCM_CGAL2DGRID_HPP

#include <mpi.h>

#include <lib/grid/Grid.hpp>
#include <lib/linal/linal.hpp>
#include <lib/util/storage/Cgal2DTriangulationStorage.hpp>
#include <lib/util/snapshot/VtkCgal2DSnapshotter.hpp>


namespace gcm {
	template<class TGrid> class DefaultSolver;
	template<class TGrid> class GridCharacteristicMethod;

	/**
	 * 2D movable unstructured triangle grid
	 * @tparam TModel rheology model
	 */
	template<typename TModel>
	class Cgal2DGrid : public Grid {
	public:
		typedef          TModel               Model;
		typedef typename Model::Vector        PdeVector;
		typedef typename Model::OdeVariables  OdeVariables;
		typedef typename Model::GCM_MATRICES  GCM_MATRICES;
		typedef typename GCM_MATRICES::Matrix Matrix;

		static const int DIMENSIONALITY = TModel::DIMENSIONALITY;

	protected:
		/**
		 * Data storage. Real values plus auxiliary values on borders.
		 * "...New" means on the next time layer.
		 */
		StdVectorStorage<PdeVector>     pdeVectors;
		StdVectorStorage<PdeVector>     pdeVectorsNew;
		StdVectorStorage<GCM_MATRICES*> gcmMatrices;
		StdVectorStorage<OdeVariables>  odeValues;
		/** Coordinates, cells, etc ... */
		Cgal2DTriangulationStorage      mesh;

		linal::VectorInt<3> sizes = {0, 0, 0}; // numbers of nodes along each direction (on this core)
		linal::Vector<3> startR = {0, 0, 0}; // global coordinates of the first real node of the grid
		linal::Vector<3> h = {0, 0, 0}; // spatial steps along each direction

		/** Read / write access to real actual PDE vectors */
		inline PdeVector &operator()(const size_t i) {
			return this->pdeVectors[i];
		};

	public:
		/** @return number of nodes */
		size_t size() const {
			return mesh.verticesSize();
		};
		/** Read-only access to real actual PDE vectors */
		inline const PdeVector& get(const size_t i) const {
			return this->pdeVectors[i];
		};
		/** Read-only access to real GCM matrices */
		inline GCM_MATRICES* getMatrix(const size_t i) const {
			return this->gcmMatrices[i];
		};

	protected:

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
		Matrix interpolateValuesAround(const int stage, const size_t i, const PdeVector& dx) const;

		/**
		 * Place in src nodal values which are necessary for
		 * interpolation in specified point. The number of placed
		 * in values is equal to src.size()
		 */
		void findSourcesForInterpolation(const int stage, const size_t i, const real &dx,
		                                 std::vector<PdeVector>& src) const;

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
		linal::Vector3 getCoordinates(const size_t i) const {
			return mesh.getCoordinates(i);
		};

		USE_AND_INIT_LOGGER("gcm.Cgal2DGrid");
		friend class GridCharacteristicMethod<Cgal2DGrid>;
		friend class DefaultSolver<Cgal2DGrid>;
		// it's better than create million gets
		friend class VtkCgal2DSnapshotter<Cgal2DGrid>;
	};
}

#endif // LIBGCM_CGAL2DGRID_HPP
