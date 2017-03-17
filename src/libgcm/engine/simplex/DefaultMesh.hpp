#ifndef LIBGCM_SIMPLEX_DEFAULTMESH_HPP
#define LIBGCM_SIMPLEX_DEFAULTMESH_HPP

#include <libgcm/engine/simplex/AbstractMesh.hpp>
#include <libgcm/util/task/InitialCondition.hpp>


namespace gcm {
namespace simplex {

/**
 * Mesh implements the approach when data are stored in separate std::vectors.
 * All nodes have the same type of rheology model and material.
 * @tparam TModel     rheology model
 * @tparam TGrid      instantiation of gcm::SimplexGrid
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
	/// Type of auxiliary info sent from inner gcm to border/contact correctors
	typedef typename Model::WaveIndices         WaveIndices;
	static const Models::T ModelType = Model::Type;
	
	typedef AbstractMesh<TGrid>                 Base;
	typedef typename Base::Grid                 Grid;
	typedef typename Base::GridId               GridId;
	typedef typename Base::ConstructionPack     ConstructionPack;
	typedef typename Base::RealD                RealD;
	typedef typename Base::MatrixDD             MatrixDD;
	typedef typename Grid::Iterator             Iterator;
	
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
	
	
	virtual void setUpPde(const Task& task, const MatrixDD& innerBasis,
			const BorderCalcMode borderCalcMode) override {
		assert_false(pdeIsSetUp);
		pdeIsSetUp = true;
		allocate();
		applyMaterialsCondition(task, innerBasis, borderCalcMode);
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
	
	/** Read-only access to WaveIndices */
	WaveIndices waveIndices(const Iterator& it) const {
		return this->waveIndicesData[this->getIndex(it)];
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
	
	/** Read / write access to WaveIndices */
	WaveIndices& _waveIndices(const Iterator& it) {
		return this->waveIndicesData[this->getIndex(it)];
	}
	
	
	virtual real getMaximalEigenvalue() const override {
		assert_gt(maximalEigenvalue, 0);
		return maximalEigenvalue;
	}
	
	virtual void averageNewPdeLayersToCurrent() override {
		for (Iterator it : *this) {
			_pde(it) = PdeVector::Zeros();
		}
		for (int s = 0; s < (int)numberOfNextPdeTimeLayers; s++) {
			for (Iterator it : *this) {
				_pde(it) += pdeNew(s, it) / numberOfNextPdeTimeLayers;
			}
		}
	}
	
	virtual void swapCurrAndNextPdeTimeLayer(const int indexOfNextPde) override {
		assert_lt(indexOfNextPde, (int)numberOfNextPdeTimeLayers);
		std::swap(pdeVariables, pdeVariablesNew[(size_t)indexOfNextPde]);
	}
	
	virtual void setInnerCalculationBasis(const MatrixDD& basis) override {
		Iterator someInnerNode = *(this->innerBegin());
		Model::constructGcmMatrices(_matrices(someInnerNode),
				material(someInnerNode), basis);
	}

	MatrixDD getInnerCalculationBasis() const {
		Iterator someInnerNode = *(this->innerBegin());
		return matrices(someInnerNode)->basis;
	}
	
	
protected:
	/// Data storage @{
	std::vector<PdeVariables> pdeVariables;
	std::vector<std::vector<PdeVariables>> pdeVariablesNew;
	std::vector<GcmMatricesPtr> gcmMatrices;
	std::vector<MaterialPtr> materials;
	std::vector<WaveIndices> waveIndicesData;
	/// @}
	
	/// there is only one "current" PDE time layer, but several "next"(new) layers
	size_t numberOfNextPdeTimeLayers = 0;
	/// maximal in modulus eigenvalue of all gcm matrices
	real maximalEigenvalue = 0;
	/// a way to delay with PDE data allocation
	bool pdeIsSetUp = false;
	
	
private:
	void allocate() {
		pdeVariables.resize(this->sizeOfAllNodes(), PdeVariables::Zeros());
		pdeVariablesNew.resize(numberOfNextPdeTimeLayers);
		for (auto& pdeNew : pdeVariablesNew) {
			pdeNew.resize(this->sizeOfAllNodes(), PdeVariables::Zeros());
		}
		gcmMatrices.resize(this->sizeOfAllNodes(), GcmMatricesPtr());
		materials.resize(this->sizeOfAllNodes(), MaterialPtr());
		// yes, it's not used in inner nodes at all.
		// but saving this not a big amount of memory requires much pain
		waveIndicesData.resize(this->sizeOfAllNodes());
	}
	
	void applyMaterialsCondition(const Task& task, const MatrixDD& innerBasis,
			const BorderCalcMode borderCalcMode) {
		assert_true(task.materialConditions.type ==
				Task::MaterialCondition::Type::BY_BODIES); // TODO BY_AREAS
		std::shared_ptr<AbstractMaterial> abstractMaterial =
				task.materialConditions.byBodies.bodyMaterialMap.at(this->id);
		std::shared_ptr<Material> concreteMaterial =
				std::dynamic_pointer_cast<Material>(abstractMaterial);
		assert_true(concreteMaterial);
		std::shared_ptr<GCM_MATRICES> innerGcmMatrices = std::make_shared<GCM_MATRICES>();
		Model::constructGcmMatrices(innerGcmMatrices, concreteMaterial, innerBasis);
		maximalEigenvalue = innerGcmMatrices->getMaximalEigenvalue();
		
		for (Iterator it : *this) {
			this->_material(it) = concreteMaterial;
			this->_matrices(it) = innerGcmMatrices;
		}
		
		if (borderCalcMode == BorderCalcMode::GLOBAL_BASIS) { return; }
		
		for (auto it = this->borderBegin(); it != this->borderEnd(); ++it) {
			const RealD borderNormal = this->borderNormal(*it);
			_matrices(*it) = std::make_shared<GCM_MATRICES>();
			Model::constructGcmMatrices(_matrices(*it), concreteMaterial,
					linal::createLocalBasisWithX(borderNormal));
			if (matrices(*it)->getMaximalEigenvalue() > maximalEigenvalue) {
				maximalEigenvalue = matrices(*it)->getMaximalEigenvalue();
			}
		}
		
		for (auto it = this->contactBegin(); it != this->contactEnd(); ++it) {
			const RealD contactNormal = this->contactNormal(*it);
			_matrices(*it) = std::make_shared<GCM_MATRICES>();
			Model::constructGcmMatrices(_matrices(*it), concreteMaterial,
					linal::createLocalBasisWithX(contactNormal));
			if (matrices(*it)->getMaximalEigenvalue() > maximalEigenvalue) {
				maximalEigenvalue = matrices(*it)->getMaximalEigenvalue();
			}
		}
	}
};

} // namespace simplex 
} // namespace gcm

#endif // LIBGCM_SIMPLEX_DEFAULTMESH_HPP
