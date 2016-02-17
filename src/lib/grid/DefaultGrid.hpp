#ifndef LIBGCM_DEFAULTGRID_HPP
#define LIBGCM_DEFAULTGRID_HPP

#include <lib/util/storage/StdVectorStorage.hpp>
#include <lib/util/task/InitialCondition.hpp>

namespace gcm {
	template<typename TModel, typename TGrid> class GcmHandler;
	template<typename TGrid> class Binary2DSeismograph;
	template<typename TGrid> class VtkStructuredSnapshotter;
	template<typename TGrid> class DefaultSolver;
	template<typename TGrid> class GridCharacteristicMethod;

	/**
	 * Grid that implement the approach when all nodal data are stored 
	 * in separated vectors not in one node
	 * @tparam TModel  rheology model
	 * @tparam TGrid   geometric aspects
	 */
	template<typename TModel, typename TGrid>
	class DefaultGrid : public TGrid {
	public:
		typedef          TModel               Model;
		typedef typename Model::Vector        PdeVector;
		typedef typename Model::OdeVariables  OdeVariables;
		typedef typename Model::GCM_MATRICES  GCM_MATRICES;
		typedef typename GCM_MATRICES::Matrix Matrix;

		static const int DIMENSIONALITY = TModel::DIMENSIONALITY;
		
		typedef typename TGrid::Iterator      Iterator;

		typedef GcmHandler<TModel, TGrid>     GCM_HANDLER;
		
		/** Read-only access to real actual PDE vectors */
		const PdeVector& pde(const Iterator& it) const {
			return this->pdeVectors[this->getIndex(it)];
		};
		/** Read-only access to real actual ODE values */
		const OdeVariables& ode(const Iterator& it) const {
			return this->odeValues[this->getIndex(it)];
		};
		/** Read-only access to real GCM matrices */
		GCM_MATRICES* matrix(const Iterator& it) const {
			return this->gcmMatrices[this->getIndex(it)];
		};

	protected:
		/** Read / write access to real actual PDE vectors */
		PdeVector& _pde(const Iterator& it) {
			return this->pdeVectors[this->getIndex(it)];
		};
		/** Read / write access to real actual ODE vectors */
		OdeVariables& _ode(const Iterator& it) {
			return this->odeValues[this->getIndex(it)];
		};
		/** Read / write access to real PDE vectors in auxiliary "on next time layer" storage */
		PdeVector& _pdeNew(const Iterator& it) {
			return this->pdeVectorsNew[this->getIndex(it)];
		};
		/** Read / write access to real actual GCM matrices */
		GCM_MATRICES*& _matrix(const Iterator& it) {
			return this->gcmMatrices[this->getIndex(it)];
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

		virtual void initializeImplImpl(const Task& task) override {
			this->pdeVectors.zeroInitialize(this->sizeOfAllNodes());
			this->pdeVectorsNew.zeroInitialize(this->sizeOfAllNodes());
			this->gcmMatrices.zeroInitialize(this->sizeOfAllNodes());
			if (Model::InternalOde::NonTrivial)
				this->odeValues.zeroInitialize(this->sizeOfAllNodes());

			typename TModel::Material material;
			material.initialize(task);
			auto gcmMatricesPtr = new GCM_MATRICES(material);
			for (auto& gcmMatrix : this->gcmMatrices) {
				gcmMatrix = gcmMatricesPtr;
			}
			this->maximalLambda = gcmMatricesPtr->getMaximalEigenvalue();
		};
		void applyInitialConditions(const Task& task) {
			InitialCondition<TModel> initialCondition;
			initialCondition.initialize(task);
			for (auto it : *this) {
				initialCondition.apply(_pde(it), this->coords(it));
			}
		};

		virtual void recalculateMaximalLambda() override { /* TODO for non-linear materials */ };

		friend class GcmHandler<TModel, TGrid>;
		friend class Binary2DSeismograph<DefaultGrid>;
		friend class VtkStructuredSnapshotter<DefaultGrid>;
		friend class DefaultSolver<DefaultGrid>;
		friend class GridCharacteristicMethod<DefaultGrid>;
	};
}

#endif // LIBGCM_DEFAULTGRID_HPP
