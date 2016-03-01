#ifndef LIBGCM_DEFAULTMESH_HPP
#define LIBGCM_DEFAULTMESH_HPP

#include <lib/util/task/InitialCondition.hpp>

namespace gcm {
	template<typename TMesh> class DefaultSolver;
	template<typename TMesh> class GridCharacteristicMethod;
	template<typename TModel, typename TGrid> struct GcmHandler;
	template<typename TModel, typename TGrid> struct DataBus;
	template<typename TModel, typename TGrid> struct MeshMover;
	template<typename TModel, typename TGrid> struct BorderConditions;

	/**
	 * Mesh that implement the approach when all nodal data are stored
	 * in separated vectors not in one node
	 * @tparam TModel  rheology model
	 * @tparam TGrid   geometric aspects
	 */
	template<typename TModel, typename TGrid>
	class DefaultMesh : public TGrid {
	public:
		typedef          TModel               Model;
		typedef typename Model::PdeVector     PdeVector;
		typedef typename Model::OdeVariables  OdeVariables;
		typedef typename Model::GCM_MATRICES  GCM_MATRICES;
		typedef typename GCM_MATRICES::Matrix Matrix;

		static const int DIMENSIONALITY = TModel::DIMENSIONALITY;

		typedef          TGrid                Grid;
		typedef typename Grid::Iterator       Iterator;

		typedef GcmHandler<TModel, TGrid>     GCM_HANDLER;
		
		/** Read-only access to real actual PDE vectors */
		const PdeVector& pde(const Iterator& it) const {
			return this->pdeVectors[this->getIndex(it)];
		};
		/** Read-only access to real actual ODE values */
		const OdeVariables& ode(const Iterator& it) const {
			return this->odeValues[this->getIndex(it)];
		};
		/** Read-only access to real PDE vectors on next time layer */
		const PdeVector& pdeNew(const Iterator& it) const {
			return this->pdeVectorsNew[this->getIndex(it)];
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
		std::vector<PdeVector>        pdeVectors;
		std::vector<PdeVector>        pdeVectorsNew;
		std::vector<GCM_MATRICES*>    gcmMatrices;
		std::vector<OdeVariables>     odeValues;
		
		virtual void beforeStatementImpl(const Statement& statement) override {
			// TODO - for movable meshes it should reconstruct the grid
			allocate();

			typename TModel::Material material;
			material.initialize(statement);
			auto gcmMatricesPtr = new GCM_MATRICES(material);
			for (auto& gcmMatrix : this->gcmMatrices) {
				gcmMatrix = gcmMatricesPtr;
			}
			this->maximalLambda = gcmMatricesPtr->getMaximalEigenvalue();
		};
		void applyInitialConditions(const Statement& statement) {
			InitialCondition<TModel> initialCondition;
			initialCondition.initialize(statement);
			for (auto it : *this) {
				initialCondition.apply(_pde(it), this->coords(it));
			}
		};
		virtual void afterStatement() override { };


		virtual void recalculateMaximalLambda() override { /* TODO for non-linear materials */ };

		void allocate() {
			zeroInitialize(pdeVectors, this->sizeOfAllNodes());
			zeroInitialize(pdeVectorsNew, this->sizeOfAllNodes());
			zeroInitialize(gcmMatrices, this->sizeOfAllNodes());
			if (Model::InternalOde::NonTrivial)
				zeroInitialize(odeValues, this->sizeOfAllNodes());
		};
		template<typename T>
		static void zeroInitialize(std::vector<T>& stdVec, size_t initSize) {
			stdVec.resize(initSize);
			memset(&(stdVec[0]), 0, stdVec.size() * sizeof(T));
		};
		

		friend class DefaultSolver<DefaultMesh>;
		friend class GridCharacteristicMethod<DefaultMesh>;
		friend class GcmHandler<TModel, TGrid>;
		friend class DataBus<TModel, TGrid>;
		friend class MeshMover<TModel, TGrid>;
		friend class BorderConditions<TModel, TGrid>;
	};
}

#endif // LIBGCM_DEFAULTMESH_HPP
