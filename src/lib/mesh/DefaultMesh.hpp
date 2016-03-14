#ifndef LIBGCM_DEFAULTMESH_HPP
#define LIBGCM_DEFAULTMESH_HPP

#include <lib/util/task/InitialCondition.hpp>

namespace gcm {
	template<typename TMesh> class DefaultSolver;
	template<typename TMesh> class GridCharacteristicMethod;
	template<typename TModel, typename TGrid, typename TMaterial> struct GcmHandler;
	template<typename TModel, typename TGrid, typename TMaterial> struct DataBus;
	template<typename TModel, typename TGrid, typename TMaterial> struct MeshMover;
	template<typename TModel, typename TGrid, typename TMaterial> struct BorderConditions;

	/**
	 * Mesh that implement the approach when all nodal data are stored
	 * in separated vectors not in one node.
	 * All nodes have the same type of rheology model, variables and material.
	 * @tparam TModel     rheology model
	 * @tparam TGrid      geometric aspects
	 * @tparam TMaterial  type of material
	 */
	template<typename TModel, typename TGrid, typename TMaterial>
	class DefaultMesh : public TGrid {
	public:
		typedef          TModel                        Model;
		typedef typename Model::PdeVector              PdeVector;
		typedef typename Model::OdeVariables           OdeVariables;
		typedef typename Model::GCM_MATRICES           GCM_MATRICES;
		typedef typename Model::GcmMatricesPtr         GcmMatricesPtr;
		typedef typename Model::ConstGcmMatricesPtr    ConstGcmMatricesPtr;
		typedef typename GCM_MATRICES::Matrix          Matrix;

		typedef          TGrid                         Grid;
		typedef typename Grid::Iterator                Iterator;

		typedef          TMaterial                     Material;
		typedef std::shared_ptr<Material>              MaterialPtr;
		typedef std::shared_ptr<const Material>        ConstMaterialPtr;

		typedef GcmHandler<Model, Grid, Material>      GCM_HANDLER; // todo - unify with others

		// Dimensionality of rheology model, the grid can have different
		static const int DIMENSIONALITY = TModel::DIMENSIONALITY;


		DefaultMesh(const Task& task) : Grid(task) { }
		virtual ~DefaultMesh() { }

		/** Read-only access to actual PDE vectors */
		const PdeVector& pde(const Iterator& it) const {
			return this->pdeVectors[this->getIndex(it)];
		}
		/** Read-only access to actual ODE values */
		const OdeVariables& ode(const Iterator& it) const {
			return this->odeValues[this->getIndex(it)];
		}
		/** Read-only access to PDE vectors on next time layer */
		const PdeVector& pdeNew(const Iterator& it) const {
			return this->pdeVectorsNew[this->getIndex(it)];
		}
		/** Read-only access to GCM matrices */
		ConstGcmMatricesPtr matrix(const Iterator& it) const {
			return this->gcmMatrices[this->getIndex(it)];
		}
		/** Read-only access to material */
		ConstMaterialPtr material(const Iterator& it) const {
			return this->materials[this->getIndex(it)];
		}

	protected:
		/** Read / write access to actual PDE vectors */
		PdeVector& _pde(const Iterator& it) {
			return this->pdeVectors[this->getIndex(it)];
		}
		/** Read / write access to actual ODE vectors */
		OdeVariables& _ode(const Iterator& it) {
			return this->odeValues[this->getIndex(it)];
		}
		/** Read / write access to PDE vectors in auxiliary "on next time layer" storage */
		PdeVector& _pdeNew(const Iterator& it) {
			return this->pdeVectorsNew[this->getIndex(it)];
		}
		/** Read / write access to actual GCM matrices */
		GcmMatricesPtr& _matrix(const Iterator& it) {
			return this->gcmMatrices[this->getIndex(it)];
		}
		/** Read / write access to actual GCM matrices */
		MaterialPtr& _material(const Iterator& it) {
			return this->materials[this->getIndex(it)];
		}

		/**
		 * Data storage. Real values plus auxiliary values on borders.
		 * "...New" means on the next time layer.
		 */
		std::vector<PdeVector>        pdeVectors;
		std::vector<PdeVector>        pdeVectorsNew;
		std::vector<GcmMatricesPtr>   gcmMatrices;
		std::vector<MaterialPtr>      materials;
		std::vector<OdeVariables>     odeValues;
		
		virtual void beforeStatementImpl(const Statement& statement) override {
			// TODO - for movable meshes it should reconstruct the grid
			allocate();
			applyRheologyConditions(statement);
			applyInitialConditions(statement);
		}
		void allocate() {
			pdeVectors.resize(this->sizeOfAllNodes(), PdeVector::zeros());
			pdeVectorsNew.resize(this->sizeOfAllNodes(), PdeVector::zeros());
			gcmMatrices.resize(this->sizeOfAllNodes(), GcmMatricesPtr());
			materials.resize(this->sizeOfAllNodes(), MaterialPtr());
			if (Model::InternalOde::NonTrivial) {
				odeValues.resize(this->sizeOfAllNodes());
			}
		}
		void applyRheologyConditions(const Statement& statement) {
			Material* material_ = new Material();
			material_->initialize(statement);
			std::shared_ptr<Material> mp(material_);
			for (auto& materialPtr : this->materials) {
				materialPtr = mp;
			}

			auto gcmMatricesPtr = std::make_shared<GCM_MATRICES>();
			Model::constructGcmMatrices(gcmMatricesPtr, PdeVector::zeros(), *material_);
			for (auto& gcmMatrix : this->gcmMatrices) {
				gcmMatrix = gcmMatricesPtr;
			}
			this->maximalLambda = gcmMatricesPtr->getMaximalEigenvalue();
		}
		void applyInitialConditions(const Statement& statement) {
			InitialCondition<Model, Material> initialCondition;
			initialCondition.initialize(statement);
			for (auto it : *this) {
				initialCondition.apply(_pde(it), this->coords(it));
			}
		}

		virtual void afterStatement() override { }

		virtual void recalculateMaximalLambda() override { /* TODO for non-linear materials */ }
		

		friend class DefaultSolver<DefaultMesh>;
		friend class GridCharacteristicMethod<DefaultMesh>;
		friend class GcmHandler<Model, Grid, Material>;
		friend class DataBus<Model, Grid, Material>;
		friend class MeshMover<Model, Grid, Material>;
		friend class BorderConditions<Model, Grid, Material>;
	};
}

#endif // LIBGCM_DEFAULTMESH_HPP
