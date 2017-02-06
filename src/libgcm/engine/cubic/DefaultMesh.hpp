#ifndef LIBGCM_CUBIC_DEFAULTMESH_HPP
#define LIBGCM_CUBIC_DEFAULTMESH_HPP

#include <libgcm/engine/cubic/AbstractMesh.hpp>
#include <libgcm/util/task/MaterialsCondition.hpp>
#include <libgcm/util/task/InitialCondition.hpp>


namespace gcm {
namespace cubic {

/**
 * Mesh implements the approach when data are stored in separate std::vectors.
 * All nodes have the same type of rheology model and material.
 * @tparam TModel     rheology model
 * @tparam TGrid      instantiation of gcm::CubicGrid
 * @tparam TMaterial  type of material
 * @see comments in AbstractMesh
 */
template<typename TModel, typename TGrid, typename TMaterial>
class DefaultMesh : public AbstractMesh<TGrid> {
public:
	typedef TModel                              Model;
	typedef typename Model::PdeVariables        PdeVariables;
	typedef typename Model::PdeVector           PdeVector;
	typedef typename Model::GCM_MATRICES        GCM_MATRICES;
	typedef typename Model::GcmMatricesPtr      GcmMatricesPtr;
	typedef typename Model::ConstGcmMatricesPtr ConstGcmMatricesPtr;
	typedef typename GCM_MATRICES::Matrix       Matrix;
	static const Models::T ModelType = Model::Type;
	
	typedef AbstractMesh<TGrid>                 Base;
	typedef typename Base::Grid                 Grid;
	typedef typename Base::GridId               GridId;
	typedef typename Base::ConstructionPack     ConstructionPack;
	typedef typename Grid::Iterator             Iterator;
	typedef typename Grid::MatrixDD             MatrixDD;
	
	typedef TMaterial                           Material;
	typedef std::shared_ptr<Material>           MaterialPtr;
	typedef std::shared_ptr<const Material>     ConstMaterialPtr;
	static const Materials::T MaterialType = Material::Type;
	
	/// Dimensionality of rheology model and grid
	static const int DIMENSIONALITY = Model::DIMENSIONALITY;
	
	virtual Models::T getModelType() const override { return ModelType; }
	virtual Materials::T getMaterialType() const override { return MaterialType; }
	
	
	DefaultMesh(const Task&, const GridId gridId_,
			const ConstructionPack& constructionPack,
			const size_t numberOfNextPdeTimeLayers_) :
					Base(gridId_, constructionPack),
					numberOfNextPdeTimeLayers(numberOfNextPdeTimeLayers_),
					pdeIsSetUp(false) {
		static_assert(Grid::DIMENSIONALITY == Model::DIMENSIONALITY, "");
	}
	virtual ~DefaultMesh() { }
	
	virtual void setUpPde(const Task& task) override {
		assert_false(pdeIsSetUp);
		pdeIsSetUp = true;
		allocate();
		MaterialsCondition<Model, Grid, Material, DefaultMesh>::apply(task, this);
		InitialCondition<Model, Grid, Material, DefaultMesh>::apply(task, this);
	}
	
	
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
	
	/**
	 * Read-only access to PDE vectors on next time layer.
	 * @param s -- stage -- we have individual next time layer for each stage
	 */
	const PdeVector& pdeNew(const int s, const Iterator& it) const {
		return this->pdeVariablesNew[(size_t)s][this->getIndex(it)];
	}
	
	/** Read-only access to GCM matrices */
	ConstGcmMatricesPtr matrices(const Iterator& it) const {
		return this->gcmMatrices[this->getIndex(it)];
	}
	
	/** Read-only access to material */
	ConstMaterialPtr material(const Iterator& it) const {
		return this->materials[this->getIndex(it)];
	}
	
	/** Read / write access to actual PDE variables */
	PdeVariables& _pdeVars(const Iterator& it) {
		return this->pdeVariables[this->getIndex(it)];
	}
	
	/** Read / write access to actual PDE vectors */
	PdeVector& _pde(const Iterator& it) {
		return this->pdeVariables[this->getIndex(it)];
	}
	
	/**
	 * Read / write access to PDE vectors on next time layer.
	 * @param s -- stage -- we have individual next time layer for each stage
	 */
	PdeVector& _pdeNew(const int s, const Iterator& it) {
		return this->pdeVariablesNew[(size_t)s][this->getIndex(it)];
	}
	
	/** Read / write access to actual GCM matrices */
	GcmMatricesPtr& _matrices(const Iterator& it) {
		return this->gcmMatrices[this->getIndex(it)];
	}
	
	/** Read / write access to actual GCM matrices */
	MaterialPtr& _material(const Iterator& it) {
		return this->materials[this->getIndex(it)];
	}
	
	
	virtual real getMaximalEigenvalue() const override {
		assert_gt(maximalEigenvalue, 0);
		return maximalEigenvalue;
	}
	
	virtual void swapCurrAndNextPdeTimeLayer(const int indexOfNextPde) override {
		assert_lt(indexOfNextPde, (int)numberOfNextPdeTimeLayers);
		std::swap(pdeVariables, pdeVariablesNew[(size_t)indexOfNextPde]);
	}
	
	
protected:
	/// Data storage @{
	std::vector<PdeVariables> pdeVariables;
	std::vector<std::vector<PdeVariables>> pdeVariablesNew;
	std::vector<GcmMatricesPtr> gcmMatrices;
	std::vector<MaterialPtr> materials;
	/// @}
	
	/// there is only one "current" PDE time layer, but several "next"(new) layers
	size_t numberOfNextPdeTimeLayers = 0;
	/// maximal in modulus eigenvalue of all gcm matrices
	real maximalEigenvalue = 0;
	/// a way to delay with PDE data allocation
	bool pdeIsSetUp = false;
	
	
private:
	friend class MaterialsCondition<Model, Grid, Material, DefaultMesh>;
	
	void allocate() {
		pdeVariables.resize(this->sizeOfAllNodes(), PdeVariables::Zeros());
		pdeVariablesNew.resize(numberOfNextPdeTimeLayers);
		for (auto& pdeNew : pdeVariablesNew) {
			pdeNew.resize(this->sizeOfAllNodes(), PdeVariables::Zeros());
		}
		gcmMatrices.resize(this->sizeOfAllNodes(), GcmMatricesPtr());
		materials.resize(this->sizeOfAllNodes(), MaterialPtr());
	}
};


} // namespace cubic
} // namespace gcm

#endif // LIBGCM_CUBIC_DEFAULTMESH_HPP
