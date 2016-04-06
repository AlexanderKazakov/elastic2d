#ifndef LIBGCM_DEFAULTMESH_HPP
#define LIBGCM_DEFAULTMESH_HPP

#include <lib/util/task/BorderCondition.hpp>
#include <lib/util/task/InitialCondition.hpp>
#include <lib/util/task/MaterialsCondition.hpp>

namespace gcm {
template<typename TMesh> class DefaultSolver;
template<typename TMesh> class GridCharacteristicMethod;
template<typename TModel, typename TGrid, typename TMaterial> struct GcmHandler;
template<typename TModel, typename TGrid, typename TMaterial> struct DataBus;
template<typename TModel, typename TGrid, typename TMaterial> struct MeshMover;
template<typename TModel, typename TGrid, typename TMaterial> struct OldBorderConditions;

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
	typedef TModel                              Model;
	typedef BorderCondition<Model>              BORDER_CONDITION;
	typedef typename Model::PdeVector           PdeVector;
	typedef typename Model::OdeVariables        OdeVariables;
	typedef typename Model::GCM_MATRICES        GCM_MATRICES;
	typedef typename Model::GcmMatricesPtr      GcmMatricesPtr;
	typedef typename Model::ConstGcmMatricesPtr ConstGcmMatricesPtr;
	typedef typename GCM_MATRICES::Matrix       Matrix;

	typedef TGrid                               Grid;
	typedef typename Grid::Iterator             Iterator;

	typedef TMaterial                           Material;
	typedef std::shared_ptr<Material>           MaterialPtr;
	typedef std::shared_ptr<const Material>     ConstMaterialPtr;

	typedef GcmHandler<Model, Grid, Material>   GCM_HANDLER; 

	// Dimensionality of rheology model, the grid can have different
	static const int DIMENSIONALITY = TModel::DIMENSIONALITY;

	struct Node { // TODO - node is some sort of handle, it can be 
		// realized in iterator
		Node(const Iterator& iterator, DefaultMesh* const mesh_) :
			it(iterator), mesh(mesh_) { }

		const PdeVector& pde() const { return mesh->pde(it); }
		PdeVector& _pde() { return mesh->_pde(it); }

		const OdeVariables& ode() const { return mesh->ode(it); }
		OdeVariables& _ode() { return mesh->_ode(it); }

		ConstGcmMatricesPtr matrices() const { return mesh->matrices(it); }
		GcmMatricesPtr& _matrices() { return mesh->_matrices(it); }

		ConstMaterialPtr material() const { return mesh->material(it); }
		MaterialPtr& _material() { return mesh->_material(it); }

		Real3 coords() const { return mesh->coords(it); }

		/** @warning coords aren't copied */
		template<typename TNodePtr>
		void copyFrom(TNodePtr origin) {
			_pde() = origin->pde();
			_ode() = origin->ode();         // TODO!!! if ode is not present??
			_matrices() = origin->_matrices();
			_material() = origin->_material();
		}

	private:
		const Iterator it;
		DefaultMesh* const mesh;
	};

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
	ConstGcmMatricesPtr matrices(const Iterator& it) const {
		return this->gcmMatrices[this->getIndex(it)];
	}

	/** Read-only access to material */
	ConstMaterialPtr material(const Iterator& it) const {
		return this->materials[this->getIndex(it)];
	}

	real getMaximalEigenvalue() const {
		assert_gt(maximalEigenvalue, 0);
		return maximalEigenvalue;
	}

protected:
	/** Read / write "node" wrapper */
	std::shared_ptr<Node> node(const Iterator& it) {
		return std::make_shared<Node>(it, this);
	}

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
	GcmMatricesPtr& _matrices(const Iterator& it) {
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
	///@{
	std::vector<PdeVector> pdeVectors;
	std::vector<PdeVector> pdeVectorsNew;
	std::vector<GcmMatricesPtr> gcmMatrices;
	std::vector<MaterialPtr> materials;
	std::vector<OdeVariables> odeValues;
	///@}

	real maximalEigenvalue = 0; ///< maximal in modulus eigenvalue of all gcm matrices
	
	/// list of border conditions
	std::vector<std::pair<std::shared_ptr<Area>,
	                      BORDER_CONDITION> > borderConditions;
	
	/**
	 * @return border condition for specified point. 
	 * If several border conditions are overlapped here,
	 * the last in borderConditions is returned.
	 * If no border conditions in this point, nullptr returned.
	 */
	const BORDER_CONDITION* getBorderCondition(const Iterator& it) {
		const BORDER_CONDITION* ans = nullptr;
		for (const auto& bc : borderConditions) {
			if (bc.first->contains(this->coords(it))) {
				ans = &(bc.second);
			}
		}
		return ans;
	}
	                      
	void beforeStatement(const Statement& statement) {
		// TODO - for movable meshes it should reconstruct the grid
		allocate();
		applyMaterialConditions(statement);
		applyInitialConditions(statement);
		setBorderConditions(statement);
	}

private:
	void allocate();
	void applyMaterialConditions(const Statement& statement);
	void applyInitialConditions(const Statement& statement);
	void setBorderConditions(const Statement& statement);

	void recalculateMaximalLambda() { /* TODO for non-linear materials */ }
	void afterStatement() { }

	friend class DefaultSolver<DefaultMesh>;
	friend class GridCharacteristicMethod<DefaultMesh>;
	friend class GcmHandler<Model, Grid, Material>;
	friend class DataBus<Model, Grid, Material>;
	friend class MeshMover<Model, Grid, Material>;
	friend class OldBorderConditions<Model, Grid, Material>;
};


template<typename TModel, typename TGrid, typename TMaterial>
void DefaultMesh<TModel, TGrid, TMaterial>::
allocate() {
	pdeVectors.resize(this->sizeOfAllNodes(), PdeVector::zeros());
	pdeVectorsNew.resize(this->sizeOfAllNodes(), PdeVector::zeros());
	gcmMatrices.resize(this->sizeOfAllNodes(), GcmMatricesPtr());
	materials.resize(this->sizeOfAllNodes(), MaterialPtr());
	if (Model::InternalOde::NonTrivial) {
		odeValues.resize(this->sizeOfAllNodes());
	}
}


template<typename TModel, typename TGrid, typename TMaterial>
void DefaultMesh<TModel, TGrid, TMaterial>::
applyMaterialConditions(const Statement& statement) {
	MaterialsCondition<Model, Material> condition(statement);
	for (auto it :* this) {
		condition.apply(node(it));
	}
	maximalEigenvalue = condition.getMaximalEigenvalue();
}


template<typename TModel, typename TGrid, typename TMaterial>
void DefaultMesh<TModel, TGrid, TMaterial>::
applyInitialConditions(const Statement& statement) {
	InitialCondition<Model, Material> initialCondition;
	initialCondition.initialize(statement);
	for (auto it :* this) {
		initialCondition.apply(_pde(it), this->coords(it));
	}
}


template<typename TModel, typename TGrid, typename TMaterial>
void DefaultMesh<TModel, TGrid, TMaterial>::
setBorderConditions(const Statement& statement) {
	borderConditions.clear();
	for (const auto& bc : statement.borderConditions) {
		borderConditions.push_back({bc.area, BORDER_CONDITION(bc)});
	}
}


}

#endif // LIBGCM_DEFAULTMESH_HPP
