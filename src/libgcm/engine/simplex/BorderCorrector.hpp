#ifndef LIBGCM_SIMPLEX_BORDERCORRECTOR_HPP
#define LIBGCM_SIMPLEX_BORDERCORRECTOR_HPP

#include <list>

#include <libgcm/engine/mesh/AbstractMesh.hpp>
#include <libgcm/util/task/BorderCondition.hpp>
#include <libgcm/rheology/models/ElasticModel.hpp>
#include <libgcm/rheology/models/AcousticModel.hpp>
#include <libgcm/rheology/materials/IsotropicMaterial.hpp>


namespace gcm {
namespace simplex {


/**
 * @defgroup Border correctors
 * Classes for applying "outer-waves"-correction on border nodes
 * in order to satisfy some border condition.
 * @see BorderCondition
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
	 * Apply border corrector for all nodes from the list 
	 * along given direction
	 */
	virtual void apply(std::shared_ptr<AbstractMesh<TGrid>> grid,
			std::list<NodeBorder> borderNodes, const RealD& direction) = 0;
	
	
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
	 * @return correction outer wave
	 */
	template<typename PdeVector,
			typename MatrixOmega, typename MatrixB, typename VectorB>
	PdeVector
	calculateOuterWaveCorrection(const PdeVector& u,
			const MatrixOmega& Omega, const MatrixB& B, const VectorB& b) {
		
		const auto M = B * Omega;
//		if (linal::determinant(M) == 0) { return u.Zeros(); }
		checkConditionNumber(M, maxConditionNumber);
		const auto alpha = linal::solveLinearSystem(M, b - B * u);
		return Omega * alpha;
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
	
	
	USE_AND_INIT_LOGGER("gcm.simplex.BorderCorrector")
};



template<typename Model, typename Material, typename TGrid,
         typename BorderMatrixCreator>
class ConcreteBorderCorrector : public AbstractBorderCorrector<TGrid> {
public:
	
	typedef DefaultMesh<Model, TGrid, Material> Mesh;
	
	typedef AbstractBorderCorrector<TGrid>      Base;
	typedef typename Base::NodeBorder           NodeBorder;
	typedef typename Base::RealD                RealD;
	
	ConcreteBorderCorrector(const Task::BorderCondition& bc) :
			borderCondition(bc) { }
	
	virtual void apply(std::shared_ptr<AbstractMesh<TGrid>> grid,
			std::list<NodeBorder> borderNodes, const RealD& direction) override {
		
		std::shared_ptr<Mesh> mesh = std::dynamic_pointer_cast<Mesh>(grid);
		assert_true(mesh);
		const auto b = borderCondition.b();
		
		for (const NodeBorder& nodeBorder: borderNodes) {
			
			const RealD reflectionDirection = direction;
//					linal::reflectionDirection(nodeBorder.normal, direction);
			const real projection = linal::dotProduct(reflectionDirection, nodeBorder.normal);
			if (projection == 0) { continue; }
			const RealD outerDirection = reflectionDirection * Utils::sign(projection);
			const auto Omega = Model::constructOuterEigenvectors(
					mesh->material(nodeBorder.iterator),
					linal::createLocalBasis(outerDirection));
			const auto B = BorderMatrixCreator::create(
//					outerDirection);
					nodeBorder.normal);
			auto& u = mesh->_pdeNew(nodeBorder.iterator);
			
			u += this->calculateOuterWaveCorrection(u, Omega, B, b);
		}
		
	}
	
	
private:
	const BorderCondition<Model> borderCondition;
	
};



template<typename TModel>
struct FixedForceBorderMatrixCreator {
	typedef typename TModel::RealD        RealD;
	typedef typename TModel::BorderMatrix BorderMatrix;
	static BorderMatrix create(const RealD& normal) {
		return TModel::borderMatrixFixedForce(normal);
	}
};
template<typename TModel>
struct FixedVelocityBorderMatrixCreator {
	typedef typename TModel::RealD        RealD;
	typedef typename TModel::BorderMatrix BorderMatrix;
	static BorderMatrix create(const RealD& normal) {
		return TModel::borderMatrixFixedVelocity(normal);
	}
};



template<typename TGrid>
class BorderCorrectorFactory {
public:
	
	static const int DIMENSIONALITY = TGrid::DIMENSIONALITY;
	
	typedef ElasticModel<DIMENSIONALITY>     ElasticModelD;
	typedef AcousticModel<DIMENSIONALITY>    AcousticModelD;
	
	
	static std::shared_ptr<AbstractBorderCorrector<TGrid>> create(
			const Task::BorderCondition& condition,
			const Models::T model, const Materials::T material) {
		
		if (material != Materials::T::ISOTROPIC) {
			THROW_UNSUPPORTED("Unsupported material");
		}
		
		switch (condition.type) {
			case BorderConditions::T::FIXED_FORCE:
				
				switch (model) {
					case Models::T::ELASTIC:
						return std::make_shared<ConcreteBorderCorrector<
								ElasticModelD, IsotropicMaterial, TGrid,
								FixedForceBorderMatrixCreator<ElasticModelD>>>(
										condition);
					case Models::T::ACOUSTIC:
						return std::make_shared<ConcreteBorderCorrector<
								AcousticModelD, IsotropicMaterial, TGrid,
								FixedForceBorderMatrixCreator<AcousticModelD>>>(
										condition);
					default:
						THROW_INVALID_ARG("Unknown type of model");
				}
				
			case BorderConditions::T::FIXED_VELOCITY:
				
				switch (model) {
					case Models::T::ELASTIC:
						return std::make_shared<ConcreteBorderCorrector<
								ElasticModelD, IsotropicMaterial, TGrid,
								FixedVelocityBorderMatrixCreator<ElasticModelD>>>(
										condition);
					case Models::T::ACOUSTIC:
						return std::make_shared<ConcreteBorderCorrector<
								AcousticModelD, IsotropicMaterial, TGrid,
								FixedVelocityBorderMatrixCreator<AcousticModelD>>>(
										condition);
					default:
						THROW_INVALID_ARG("Unknown type of model");
				}
				
			default:
				THROW_INVALID_ARG("Unknown type of border condition");
		}
	}
	
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_BORDERCORRECTOR_HPP
