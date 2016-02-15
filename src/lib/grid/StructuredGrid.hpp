#ifndef LIBGCM_STRUCTUREDGRID_HPP
#define LIBGCM_STRUCTUREDGRID_HPP

#include <mpi.h>

#include <lib/grid/Grid.hpp>
#include <lib/linal/linal.hpp>
#include <lib/numeric/interpolation/Interpolator.hpp>
#include <lib/numeric/border_conditions/StructuredGridBorderConditions.hpp>

namespace gcm {
	template<class TGrid> class DefaultSolver;
	template<class TGrid> class GridCharacteristicMethod;
	template<class TGrid> class VtkStructuredSnapshotter;
	template<class TGrid> class Binary2DSeismograph;

	/**
	 * Non-movable structured rectangular grid
	 * @tparam TModel rheology model
	 */
	template<typename TModel>
	class StructuredGrid : public Grid {
	public:
		typedef          TModel               Model;
		typedef typename Model::Vector        PdeVector;
		typedef typename Model::OdeVariables  OdeVariables;
		typedef typename Model::GCM_MATRICES  GCM_MATRICES;
		typedef typename GCM_MATRICES::Matrix Matrix;

		static const int DIMENSIONALITY = TModel::DIMENSIONALITY;

		struct Iterator : public linal::VectorInt<3> {
			linal::VectorInt<3> bounds = {0, 0, 0};
			Iterator(const linal::VectorInt<3> indices_, const linal::VectorInt<3> bounds_) {
				(*this) = indices_;
				bounds = bounds_;
			};
			const Iterator& operator*() { return *this; };
			using linal::VectorInt<3>::Matrix;
			using linal::VectorInt<3>::operator=;
		};

		/** slow-X fast-Z memory access efficient iterator */
		struct ForwardIterator : public Iterator {
			using Iterator::Iterator;
			Iterator& operator++() {
				!increment(2) && !increment(1) && ((*this)(0)++);
				return (*this);
			};
		private:
			int increment(const int i) {
				(*this)(i) = ((*this)(i) + 1) % this->bounds(i);
				return (*this)(i);
			};
		};
		ForwardIterator begin() const { return ForwardIterator({0, 0, 0}, sizes); };
		ForwardIterator end() const { return ForwardIterator({sizes(0), 0, 0}, sizes); };

		/** slow-Z fast-X iterator */
		struct VtkIterator : public Iterator {
			using Iterator::Iterator;
			VtkIterator& operator++() {
				!increment(0) && !increment(1) && ((*this)(2)++);
				return (*this);
			};
		private:
			int increment(const int i) {
				(*this)(i) = ((*this)(i) + 1) % this->bounds(i);
				return (*this)(i);
			};
		};
		VtkIterator vtkBegin() const { return VtkIterator({0, 0, 0}, sizes); };
		VtkIterator vtkEnd() const { return VtkIterator({0, 0, sizes(2)}, sizes); };

		/** Read-only access to real actual PDE vectors */
		const PdeVector& pde(const int x, const int y, const int z) const {
			return this->pdeVectors[getIndex(x, y, z)];
		};
		const PdeVector& pde(const Iterator& it) const {
			return this->pdeVectors[getIndex(it)];
		};
		/** Read-only access to real actual ODE values */
		const OdeVariables& ode(const int x, const int y, const int z) const {
			return this->odeValues[getIndex(x, y, z)];
		};
		const OdeVariables& ode(const Iterator& it) const {
			return this->odeValues[getIndex(it)];
		};
		/** Read-only access to real GCM matrices */
		GCM_MATRICES* matrix(const int x, const int y, const int z) const {
			return this->gcmMatrices[getIndex(x, y, z)];
		};
		GCM_MATRICES* matrix(const Iterator& it) const {
			return this->gcmMatrices[getIndex(it)];
		};
		/** Read-only access to real coordinates */
		const linal::Vector3 coords(const int x, const int y, const int z) const {
			return startR + linal::plainMultiply(linal::VectorInt<3>({x, y, z}), h);
		};
		const linal::Vector3 coords(const Iterator& it) const {
			return startR + linal::plainMultiply(it, h);
		};

		size_t sizeOfRealNodes() const {
			return (size_t) linal::directProduct(sizes);
		};
		size_t sizeOfAllNodes() const {
			return (size_t)linal::directProduct(sizes + 2 * accuracyOrder * linal::VectorInt<3>({1, 1, 1}));
		};

	protected:
		/**
		 * Data storage. Real values plus auxiliary values on borders.
		 * "...New" means on the next time layer.
		 */
		StdVectorStorage<PdeVector>     pdeVectors;
		StdVectorStorage<PdeVector>     pdeVectorsNew;
		StdVectorStorage<GCM_MATRICES*> gcmMatrices;
		StdVectorStorage<OdeVariables>  odeValues;

		int accuracyOrder = 0; // order of accuracy of spatial interpolation

		linal::VectorInt<3> sizes = {0, 0, 0}; // numbers of nodes along each direction (on this core)

		linal::Vector<3> startR = {0, 0, 0}; // global coordinates of the first real node of the grid
		linal::Vector<3> h = {0, 0, 0}; // spatial steps along each direction

		StructuredGridBorderConditions<StructuredGrid> borderConditions;
		
		/** Read / write access to real actual PDE vectors */
		PdeVector& _pde(const int x, const int y, const int z) {
			return this->pdeVectors[getIndex(x, y, z)];
		};
		PdeVector& _pde(const Iterator& it) {
			return this->pdeVectors[getIndex(it)];
		};
		OdeVariables& _ode(const Iterator& it) {
			return this->odeValues[getIndex(it)];
		};
		PdeVector& _pdeNew(const int x, const int y, const int z) {
			return this->pdeVectorsNew[getIndex(x, y, z)];
		};
		PdeVector& _pdeNew(const Iterator& it) {
			return this->pdeVectorsNew[getIndex(it)];
		};
		GCM_MATRICES*& _matrix(const int x, const int y, const int z) {
			return this->gcmMatrices[getIndex(x, y, z)];
		};
		GCM_MATRICES*& _matrix(const Iterator& it) {
			return this->gcmMatrices[getIndex(it)];
		};

	public:
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
		Matrix interpolateValuesAround(const int stage, const Iterator& it, const PdeVector& dx) const;

		/**
		 * Place in src nodal values which are necessary for
		 * interpolation in specified point. The number of placed
		 * in values is equal to src.size()
		 */
		void findSourcesForInterpolation(const int stage, const Iterator& it,
		                                 const real &dx, std::vector<PdeVector>& src) const;

	protected:
		/**
		 * @param it begin() <= iterator < end()
		 * @return index in std::vector
		 */
		size_t getIndex(const Iterator& it) const {
			return getIndex(it(0), it(1), it(2));
		};
		/**
		 * @param x x index < sizes(0)
		 * @param y y index < sizes(1)
		 * @param z z index < sizes(2)
		 * @return index in std::vector
		 */
		size_t getIndex(const int x, const int y, const int z) const {
			return (size_t)
			         (2 * accuracyOrder + sizes(2)) * (2 * accuracyOrder + sizes(1)) * (x + accuracyOrder)
			       + (2 * accuracyOrder + sizes(2)) * (y + accuracyOrder)
			       + (z + accuracyOrder);
		};


		/* Equal-distance spatial interpolator */
		Interpolator<PdeVector> interpolator;

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

		USE_AND_INIT_LOGGER("gcm.StructuredGrid");
		friend class GridCharacteristicMethod<StructuredGrid>;
		friend class DefaultSolver<StructuredGrid>;
		friend class StructuredGridBorderConditions<StructuredGrid>;
		// it's better than create million gets
		friend class VtkStructuredSnapshotter<StructuredGrid>;
		friend class Binary2DSeismograph<StructuredGrid>;
	};
}

#endif // LIBGCM_STRUCTUREDGRID_HPP
