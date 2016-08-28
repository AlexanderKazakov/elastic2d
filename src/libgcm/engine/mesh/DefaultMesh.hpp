#ifndef LIBGCM_DEFAULTMESH_HPP
#define LIBGCM_DEFAULTMESH_HPP

#include <libgcm/engine/mesh/AbstractMesh.hpp>
#include <libgcm/util/task/InitialCondition.hpp>
#include <libgcm/util/task/MaterialsCondition.hpp>


namespace gcm {

namespace cubic {
template<typename> class DataBus;
}


/**
 * Mesh implements the approach when data are stored in separate vectors.
 * All nodes have the same type of rheology model and material.
 * @tparam TModel     rheology model
 * @tparam TGrid      geometric aspects
 * @tparam TMaterial  type of material
 */
template<typename TModel, typename TGrid, typename TMaterial>
class DefaultMesh : public AbstractMesh<TGrid> {
public:
	typedef TModel                              Model;
	typedef typename Model::PdeVariables        PdeVariables;
	typedef typename Model::PdeVector           PdeVector;
//	typedef typename Model::OdeVariables        OdeVariables;
	typedef typename Model::GCM_MATRICES        GCM_MATRICES;
	typedef typename Model::GcmMatricesPtr      GcmMatricesPtr;
	typedef typename Model::ConstGcmMatricesPtr ConstGcmMatricesPtr;
	typedef typename GCM_MATRICES::Matrix       Matrix;
	
	typedef AbstractMesh<TGrid>                 Base;
	typedef typename Base::Grid                 Grid;
	typedef typename Base::GridId               GridId;
	typedef typename Base::ConstructionPack     ConstructionPack;
	typedef typename Grid::Iterator             Iterator;
	typedef typename Grid::MatrixDD             MatrixDD;
	
	typedef TMaterial                           Material;
	typedef std::shared_ptr<Material>           MaterialPtr;
	typedef std::shared_ptr<const Material>     ConstMaterialPtr;
	
	/// Dimensionality of rheology model and grid
	static const int DIMENSIONALITY = Model::DIMENSIONALITY;
	
	
	DefaultMesh(const Task& task, const GridId gridId_,
			const ConstructionPack& constructionPack) :
					Base(gridId_, constructionPack) {
		static_assert(Grid::DIMENSIONALITY == Model::DIMENSIONALITY, "");
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
	
	virtual real getMaximalEigenvalue() const override {
		assert_gt(maximalEigenvalue, 0);
		return maximalEigenvalue;
	}
	
	virtual void changeCalculationBasis(const MatrixDD& basis) override;
	
	virtual void swapPdeTimeLayers() override {
		std::swap(pdeVariables, pdeVariablesNew);
	}
	
	MatrixDD getCalculationBasis() const {
		return calculationBasis;
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
	
	
private:
	
	void allocate();
	
	friend class cubic::DataBus<DefaultMesh<Model, Grid, Material>>;
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
