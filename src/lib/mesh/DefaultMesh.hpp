#ifndef LIBGCM_DEFAULTMESH_HPP
#define LIBGCM_DEFAULTMESH_HPP

#include <lib/util/task/InitialCondition.hpp>
#include <lib/util/task/MaterialsCondition.hpp>

namespace gcm {

template<typename> class DefaultSolver;
template<typename, typename, typename> class DataBus;
template<typename, typename, typename> class GridCharacteristicMethod;

/**
 * Mesh that implements the approach when all nodal data are stored
 * in separated vectors not in one node.
 * (Mesh with any other approach (CalcNodes struct, ..) can be done,
 * but it must have the same interfaces)
 * All nodes have the same type of rheology model, variables and material.
 * @tparam TModel     rheology model
 * @tparam TGrid      geometric aspects
 * @tparam TMaterial  type of material
 */
template<typename TModel, typename TGrid, typename TMaterial>
class DefaultMesh : public TGrid {
public:
	typedef TModel                              Model;
	typedef typename Model::PdeVariables        PdeVariables;
	typedef typename Model::PdeVector           PdeVector;
//	typedef typename Model::OdeVariables        OdeVariables;
	typedef typename Model::GCM_MATRICES        GCM_MATRICES;
	typedef typename Model::GcmMatricesPtr      GcmMatricesPtr;
	typedef typename Model::ConstGcmMatricesPtr ConstGcmMatricesPtr;
	typedef typename GCM_MATRICES::Matrix       Matrix;
	
	typedef TGrid                               Grid;
	typedef typename Grid::GridId               GridId;
	typedef typename Grid::GlobalScene          GlobalScene;
	typedef typename Grid::Iterator             Iterator;
	typedef typename Grid::MatrixDD             MatrixDD;
	
	typedef TMaterial                           Material;
	typedef std::shared_ptr<Material>           MaterialPtr;
	typedef std::shared_ptr<const Material>     ConstMaterialPtr;
	
	/// Dimensionality of rheology model, the grid can have different FIXME - make the same only
	static const int DIMENSIONALITY = Model::DIMENSIONALITY;
	/// Dimensionality of the grid, rheology model can have different
	static const int GRID_DIMENSIONALITY = Grid::DIMENSIONALITY;
	
//	struct Node {
//		Node(const Iterator& iterator, DefaultMesh* const mesh_) :
//			it(iterator), mesh(mesh_) { }
	
//		const PdeVector& pde() const { return mesh->pde(it); }
//		PdeVector& _pde() { return mesh->_pde(it); }
	
//		const OdeVariables& ode() const { return mesh->ode(it); }
//		OdeVariables& _ode() { return mesh->_ode(it); }
	
//		ConstGcmMatricesPtr matrices() const { return mesh->matrices(it); }
//		GcmMatricesPtr& _matrices() { return mesh->_matrices(it); }
	
//		ConstMaterialPtr material() const { return mesh->material(it); }
//		MaterialPtr& _material() { return mesh->_material(it); }
	
//		Real3 coords() const { return mesh->coords(it); }
	
//		/** @warning coords aren't copied */
//		template<typename TNodePtr>
//		void copyFrom(TNodePtr origin) {
//			_pde() = origin->pde();
//			_ode() = origin->ode();         // TODO !!! if ode is not present??
//			_matrices() = origin->_matrices();
//			_material() = origin->_material();
//		}
	
//	private:
//		const Iterator it;
//		DefaultMesh* const mesh;
//	};
	
	DefaultMesh(const Task& task, GlobalScene* gS, const GridId gridId_) :
			Grid(task, gS, gridId_) {
		allocate();
		MaterialsCondition<Model, Grid, Material, DefaultMesh>::apply(task, this);
		InitialCondition<Model, Grid, Material, DefaultMesh>::apply(task, this);
	}
	virtual ~DefaultMesh() { }
	
	
	/** Read-only access to actual PDE variables */
	const PdeVariables& pdeVars(const Iterator& it) const {
		return this->pdeVariables[this->getIndex(it)];
	}
	
	/** 
	 * Read-only access to actual PDE vectors.
	 * Yes, it has to be a different function from pdeVars.
	 */
	const PdeVector& pde(const Iterator& it) const {
		return this->pdeVariables[this->getIndex(it)];
	}
	
//	/** Read-only access to actual ODE values */
//	const OdeVariables& ode(const Iterator& it) const {
//		assert_eq(odeVariables.size(), this->sizeOfAllNodes()); // TODO - less expensive
//		return this->odeVariables[this->getIndex(it)];
//	}
	
	/** Read-only access to PDE vectors on next time layer */
	const PdeVector& pdeNew(const Iterator& it) const {
		return this->pdeVariablesNew[this->getIndex(it)];
	}
	
//	/** Read-only access to PDE vectors on current time layer */
//	const PdeVector& pdePrev(const Iterator& it) const {
//		return this->pdeVariablesPrev[this->getIndex(it)];
//	}
	
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
	
//	/** Read / write "node" wrapper */
//	std::shared_ptr<Node> node(const Iterator& it) {
//		return std::make_shared<Node>(it, this);
//	}
	
	/** Read-only access to actual PDE variables */
	PdeVariables& _pdeVars(const Iterator& it) {
		return this->pdeVariables[this->getIndex(it)];
	}
	
	/** Read / write access to actual PDE vectors */
	PdeVector& _pde(const Iterator& it) {
		return this->pdeVariables[this->getIndex(it)];
	}
	
//	/** Read / write access to actual ODE vectors */
//	OdeVariables& _ode(const Iterator& it) {
//		assert_eq(odeVariables.size(), this->sizeOfAllNodes()); // TODO - less expensive
//		return this->odeVariables[this->getIndex(it)];
//	}
	
	/** Read / write access to PDE vectors in auxiliary "on next time layer" storage */
	PdeVector& _pdeNew(const Iterator& it) {
		return this->pdeVariablesNew[this->getIndex(it)];
	}
	
	/** Read / write access to actual GCM matrices */
	GcmMatricesPtr& _matrices(const Iterator& it) {
		return this->gcmMatrices[this->getIndex(it)];
	}
	
	/** Read / write access to actual GCM matrices */
	MaterialPtr& _material(const Iterator& it) {
		return this->materials[this->getIndex(it)];
	}
	
	
	
protected:
	/**
	 * Data storage
	 */
	///@{
	std::vector<PdeVariables> pdeVariables;
	std::vector<PdeVariables> pdeVariablesNew;
//	std::vector<PdeVariables> pdeVariablesPrev;
	std::vector<GcmMatricesPtr> gcmMatrices;
	std::vector<MaterialPtr> materials;
//	std::vector<OdeVariables> odeVariables;
	///@}
	
	MatrixDD calculationBasis = MatrixDD::Identity();
	real maximalEigenvalue = 0; ///< maximal in modulus eigenvalue of all gcm matrices
	
	virtual void changeCalculationBasis(const MatrixDD& basis) override;
	
private:
	
	void allocate();
	
	friend class DefaultSolver<DefaultMesh>;
	friend class DataBus<Model, Grid, Material>;
	friend class GridCharacteristicMethod<Model, Grid, Material>;
	friend class MaterialsCondition<Model, Grid, Material, DefaultMesh>;
	
};


template<typename TModel, typename TGrid, typename TMaterial>
void DefaultMesh<TModel, TGrid, TMaterial>::
allocate() {
	pdeVariables.resize(this->sizeOfAllNodes(), PdeVariables::Zeros());
	pdeVariablesNew.resize(this->sizeOfAllNodes(), PdeVariables::Zeros());
//	pdeVariablesPrev.resize(this->sizeOfAllNodes(), PdeVariables::Zeros());
	gcmMatrices.resize(this->sizeOfAllNodes(), GcmMatricesPtr());
	materials.resize(this->sizeOfAllNodes(), MaterialPtr());
//	if (Model::InternalOde::NonTrivial) {
//		odeVariables.resize(this->sizeOfAllNodes());
//	}
}


template<typename TModel, typename TGrid, typename TMaterial>
void DefaultMesh<TModel, TGrid, TMaterial>::
changeCalculationBasis(const MatrixDD& basis) {
	
	calculationBasis = basis;
	
	// FIXME - here we suppose mesh homogenity!
	Iterator someNode = *(this->begin());
	
	Model::constructGcmMatrices(this->_matrices(someNode),
			this->material(someNode), calculationBasis);
}


}

#endif // LIBGCM_DEFAULTMESH_HPP
