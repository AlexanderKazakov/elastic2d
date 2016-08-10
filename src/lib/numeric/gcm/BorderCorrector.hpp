#ifndef LIBGCM_BORDERCORRECTOR_HPP
#define LIBGCM_BORDERCORRECTOR_HPP

#include <list>
#include <lib/mesh/grid/AbstractGrid.hpp>

#include <lib/mesh/DefaultMesh.hpp>
#include <lib/util/task/BorderCondition.hpp>
#include <lib/rheology/models/ElasticModel.hpp>
#include <lib/rheology/models/AcousticModel.hpp>
#include <lib/rheology/materials/IsotropicMaterial.hpp>


namespace gcm {

/**
 * @defgroup Border correctors
 * Classes for applying "outer-waves"-correction on border nodes
 * in order to satisfy some border condition.
 */

template<typename TGrid>
class AbstractBorderCorrector {
public:
	
	typedef typename TGrid::Iterator    Iterator;
	typedef typename TGrid::RealD       RealD;
	
	struct NodeBorder {
		/// iterator of the node in grid
		Iterator iterator;
		/// border normal (the direction is outside the grid)
		RealD normal;
	};
	
	
	/**
	 * Apply contact corrector for all nodes from the list
	 */
	virtual void apply(AbstractGrid* grid,
			std::list<NodeBorder> borderNodes) = 0;
	
	
protected:
	/// maximal found condition number of corrector matrices
	real maxConditionNumber = 0;
	
	
	/**
	 * General expression of linear border condition is:
	 *     B * u = b,
	 * where u is pde-vector and B is border matrix.
	 * Given with inner-calculated pde vector, we correct them
	 * with outer waves combination (Omega) in order to satisfy border condition.
	 * @see BorderCondition
	 */
	template<typename PdeVector,
			typename MatrixOmega, typename MatrixB, typename VectorB>
	void
	correctBorder(PdeVector& u,
			const MatrixOmega& Omega, const MatrixB& B, const VectorB& b) {
		
		const auto M = B * Omega;
		checkConditionNumber(M, maxConditionNumber);
		const auto alpha = linal::solveLinearSystem(M, b - B * u);
		u += Omega * alpha;
	}
	
	
private:
	
	template<typename MatrixT>
	void
	checkConditionNumber(const MatrixT& m, real& currentMax) {
		const real currentValue = linal::conditionNumber(m);
		if (currentValue > currentMax) {
			currentMax = currentValue;
			LOG_INFO("New maximal condition number in border matrix: " << currentMax);
		}
	}
	
	
	USE_AND_INIT_LOGGER("gcm.BorderCorrector")
};



template<typename Model, typename Material, typename TGrid>
class FixedForceBorderCorrector : public AbstractBorderCorrector<TGrid> {
public:
	
	typedef DefaultMesh<Model, TGrid, Material> Mesh;
	
	typedef AbstractBorderCorrector<TGrid>      Base;
	typedef typename Base::NodeBorder           NodeBorder;
	
	FixedForceBorderCorrector(const Statement::BorderCondition& bc) :
			borderCondition(bc) { }
	
	virtual void apply(AbstractGrid* grid, 
			std::list<NodeBorder> borderNodes) override {
		
		Mesh* mesh = dynamic_cast<Mesh*>(grid);
		assert_true(mesh);
		const auto b = borderCondition.b();
		
		for (const NodeBorder& nodeBorder: borderNodes) {
			
			const auto Omega = Model::constructOuterEigenvectors(
					mesh->material(nodeBorder.iterator),
					linal::createLocalBasis(nodeBorder.normal));
			const auto B = Model::borderMatrixFixedForce(
					nodeBorder.normal);
			auto& u = mesh->_pdeNew(nodeBorder.iterator);
			
			this->correctBorder(u, Omega, B, b);
		}
		
	}
	
private:
	const BorderCondition<Model> borderCondition;
};



template<typename Model, typename Material, typename TGrid>
class FixedVelocityBorderCorrector : public AbstractBorderCorrector<TGrid> {
public:
	
	typedef DefaultMesh<Model, TGrid, Material> Mesh;
	
	typedef AbstractBorderCorrector<TGrid>      Base;
	typedef typename Base::NodeBorder           NodeBorder;
	
	FixedVelocityBorderCorrector(const Statement::BorderCondition& bc) :
			borderCondition(bc) { }
	
	virtual void apply(AbstractGrid* grid, 
			std::list<NodeBorder> borderNodes) override {
		
		Mesh* mesh = dynamic_cast<Mesh*>(grid);
		assert_true(mesh);
		const auto b = borderCondition.b();
		
		for (const NodeBorder& nodeBorder: borderNodes) {
			
			const auto Omega = Model::constructOuterEigenvectors(
					mesh->material(nodeBorder.iterator),
					linal::createLocalBasis(nodeBorder.normal));
			const auto B = Model::borderMatrixFixedVelocity(
					nodeBorder.normal);
			auto& u = mesh->_pdeNew(nodeBorder.iterator);
			
			this->correctBorder(u, Omega, B, b);
		}
		
	}
	
private:
	const BorderCondition<Model> borderCondition;
};



template<typename TGrid>
class BorderCorrectorFactory {
public:
	
	static const int DIMENSIONALITY = TGrid::DIMENSIONALITY;
	
	typedef ElasticModel<DIMENSIONALITY>     ElasticModelD;
	typedef AcousticModel<DIMENSIONALITY>    AcousticModelD;
	
	
	static std::shared_ptr<AbstractBorderCorrector<TGrid>> create(
			const Statement::BorderCondition& condition,
			const Models::T model, const Materials::T material) {
		
		if (material != Materials::T::ISOTROPIC) {
			THROW_UNSUPPORTED("Unsupported material");
		}
		
		switch (condition.type) {
			case BorderConditions::T::FIXED_FORCE:
				
				switch (model) {
					case Models::T::ELASTIC:
						return std::make_shared<FixedForceBorderCorrector<
								ElasticModelD, IsotropicMaterial, TGrid>>(condition);
					case Models::T::ACOUSTIC:
						return std::make_shared<FixedForceBorderCorrector<
								AcousticModelD, IsotropicMaterial, TGrid>>(condition);
					default:
						THROW_INVALID_ARG("Unknown type of model");
				}
				
			case BorderConditions::T::FIXED_VELOCITY:
				
				switch (model) {
					case Models::T::ELASTIC:
						return std::make_shared<FixedVelocityBorderCorrector<
								ElasticModelD, IsotropicMaterial, TGrid>>(condition);
					case Models::T::ACOUSTIC:
						return std::make_shared<FixedVelocityBorderCorrector<
								AcousticModelD, IsotropicMaterial, TGrid>>(condition);
					default:
						THROW_INVALID_ARG("Unknown type of model");
				}
				
			default:
				THROW_INVALID_ARG("Unknown type of border condition");
		}
	}
	
};


}


#endif // LIBGCM_BORDERCORRECTOR_HPP
