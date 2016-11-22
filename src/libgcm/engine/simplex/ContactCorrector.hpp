#ifndef LIBGCM_SIMPLEX_CONTACTCORRECTOR_HPP
#define LIBGCM_SIMPLEX_CONTACTCORRECTOR_HPP

#include <list>

#include <libgcm/engine/mesh/DefaultMesh.hpp>
#include <libgcm/rheology/models/ElasticModel.hpp>
#include <libgcm/rheology/models/AcousticModel.hpp>
#include <libgcm/rheology/materials/IsotropicMaterial.hpp>


namespace gcm {
namespace simplex {


/**
 * @defgroup Contact correctors
 * Classes for applying "outer-waves"-correction on contact nodes
 * in order to satisfy some contact condition.
 */

template<typename TGrid>
class AbstractContactCorrector {
public:
	
	typedef typename TGrid::Iterator    Iterator;
	typedef typename TGrid::RealD       RealD;
	
	struct NodesContact {
		/// pair of nodes in contact (their iterators in grids)
		Iterator first, second;
		/// contact normal (the direction is from first to second)
		RealD normal;
	};
	
	
	/**
	 * Apply contact corrector for all nodes from the list
	 * along given direction
	 */
	virtual void apply(const int s,
			std::shared_ptr<AbstractMesh<TGrid>> a,
			std::shared_ptr<AbstractMesh<TGrid>> b,
			std::list<NodesContact> nodesInContact,
			const RealD& direction) = 0;
	
	
protected:
	/// maximal found condition numbers of correctors matrices
	real maxConditionR = 0;
	real maxConditionA = 0;
	
	
	/**
	 * General expression of linear contact condition is:
	 *     B1_A * u_A = B1_B * u_B,
	 *     B2_A * u_A = B2_B * u_B,
	 * where u is pde-vector and B is border matrices.
	 * Given with inner-calculated pde vectors, we correct them
	 * with outer waves combination (Omega) in order to satisfy contact condition.
	 * @see BorderCorrector
	 */
	template<typename PdeVector, typename MatrixOmega, typename MatrixB>
	void
	correctNodesContact(
			PdeVector& uA,
			const MatrixOmega& OmegaA, const MatrixB& B1A, const MatrixB& B2A,
			PdeVector& uB,
			const MatrixOmega& OmegaB, const MatrixB& B1B, const MatrixB& B2B) {
		
		const auto R = linal::invert(B1A * OmegaA);
//		checkConditionNumber(R, "R", maxConditionR);
		const auto p = R * (B1B * uB - B1A * uA);
		const auto Q = R * (B1B * OmegaB);
		
		const auto A = (B2B * OmegaB) - ((B2A * OmegaA) * Q);
//		checkConditionNumber(A, "A", maxConditionA);
		const auto f = ((B2A * OmegaA) * p) + (B2A * uA) - (B2B * uB);
		
		const auto alphaB = linal::solveLinearSystem(A, f);
		const auto alphaA = p + Q * alphaB;
		
		uA += OmegaA * alphaA;
		uB += OmegaB * alphaB;
	}
	
	
private:
	
	template<typename MatrixT>
	void
	checkConditionNumber(const MatrixT& m, const std::string name, real& currentMax) {
		const real currentValue = linal::conditionNumber(m);
		if (currentValue > currentMax) {
			currentMax = currentValue;
			LOG_INFO("New maximal condition number in matrix "
					<< name << ": " << currentMax);
		}
	}
	
	
	USE_AND_INIT_LOGGER("gcm.simplex.ContactCorrector")
};



template<typename ModelA, typename MaterialA,
         typename ModelB, typename MaterialB,
         typename TGrid,
         typename ContactMatrixCreator>
class ConcreteContactCorrector : public AbstractContactCorrector<TGrid> {
public:
	
	typedef DefaultMesh<ModelA, TGrid, MaterialA> MeshA;
	typedef DefaultMesh<ModelA, TGrid, MaterialB> MeshB;
	
	typedef AbstractContactCorrector<TGrid> Base;
	typedef typename Base::NodesContact     NodesContact;
	typedef typename Base::RealD            RealD;
	
	virtual void apply(const int s,
			std::shared_ptr<AbstractMesh<TGrid>> a,
			std::shared_ptr<AbstractMesh<TGrid>> b,
			std::list<NodesContact> nodesInContact,
			const RealD& direction) override {
		
		std::shared_ptr<MeshA> meshA = std::dynamic_pointer_cast<MeshA>(a);
		assert_true(meshA);
		std::shared_ptr<MeshB> meshB = std::dynamic_pointer_cast<MeshB>(b);
		assert_true(meshB);
		
		for (const NodesContact& nodesContact : nodesInContact) {
			
			const real projection = linal::dotProduct(direction, nodesContact.normal);
			if (std::fabs(projection) < EQUALITY_TOLERANCE) { continue; }
			
			const RealD directionFromAToB = direction * Utils::sign(projection);
			const auto OmegaA = ModelA::constructOuterEigenvectors(
					meshA->material(nodesContact.first),
					linal::createLocalBasis(  directionFromAToB));
			const auto OmegaB = ModelB::constructOuterEigenvectors(
					meshB->material(nodesContact.second),
					linal::createLocalBasis( -directionFromAToB));
			
			const auto B1A = ContactMatrixCreator::createB1A(nodesContact.normal);
			const auto B1B = ContactMatrixCreator::createB1B(nodesContact.normal);
			const auto B2A = ContactMatrixCreator::createB2A(nodesContact.normal);
			const auto B2B = ContactMatrixCreator::createB2B(nodesContact.normal);
			
			auto& uA = meshA->_pdeNew(s, nodesContact.first);
			auto& uB = meshB->_pdeNew(s, nodesContact.second);
			
			this->correctNodesContact(uA, OmegaA, B1A, B2A,
			                          uB, OmegaB, B1B, B2B);
		}
	}
	
};



template<typename ModelA, typename ModelB>
struct AdhesionContactMatrixCreator {
	typedef typename ModelA::RealD        RealD;
	typedef typename ModelA::BorderMatrix BorderMatrix;
	
	static BorderMatrix createB1A(const RealD& normal) {
		return ModelA::borderMatrixFixedVelocityGlobalBasis(normal);
	}
	static BorderMatrix createB1B(const RealD& normal) {
		return ModelB::borderMatrixFixedVelocityGlobalBasis(normal);
	}
	static BorderMatrix createB2A(const RealD& normal) {
		return ModelA::borderMatrixFixedForceGlobalBasis(normal);
	}
	static BorderMatrix createB2B(const RealD& normal) {
		return ModelB::borderMatrixFixedForceGlobalBasis(normal);
	}
};
template<typename ModelA, typename ModelB>
struct SlideContactMatrixCreator {
	typedef typename ModelA::RealD        RealD;
	typedef typename ModelA::BorderMatrix BorderMatrix;
	
	// FIXME - this is valid for acoustic model only
	static BorderMatrix createB1A(const RealD& normal) {
		return ModelA::borderMatrixFixedVelocity(normal);
	}
	static BorderMatrix createB1B(const RealD& normal) {
		return ModelB::borderMatrixFixedVelocity(normal);
	}
	static BorderMatrix createB2A(const RealD& normal) {
		return ModelA::borderMatrixFixedForce(normal);
	}
	static BorderMatrix createB2B(const RealD& normal) {
		return ModelB::borderMatrixFixedForce(normal);
	}
};



template<typename TGrid>
class ContactCorrectorFactory {
public:
	
	static const int DIMENSIONALITY = TGrid::DIMENSIONALITY;
	
	typedef ElasticModel<DIMENSIONALITY>     ElasticModelD;
	typedef AcousticModel<DIMENSIONALITY>    AcousticModelD;
	
	
	static std::shared_ptr<AbstractContactCorrector<TGrid>> create(
			const ContactConditions::T condition,
			const Models::T model1, const Materials::T material1,
			const Models::T model2, const Materials::T material2) {
		
		switch (condition) {
			case ContactConditions::T::ADHESION:
				if (model1 == Models::T::ELASTIC &&
				    model2 == Models::T::ELASTIC &&
				    material1 == Materials::T::ISOTROPIC &&
				    material2 == Materials::T::ISOTROPIC) {
					
					return std::make_shared<ConcreteContactCorrector<
							ElasticModelD, IsotropicMaterial,
							ElasticModelD, IsotropicMaterial, TGrid,
							AdhesionContactMatrixCreator<ElasticModelD, ElasticModelD>>>();
					
				} else {
					THROW_UNSUPPORTED("Incompatible or unsupported contact conditions, \
							models and materials combination");
				}
				
			case ContactConditions::T::SLIDE:
				if (model1 == Models::T::ACOUSTIC &&
				    model2 == Models::T::ACOUSTIC &&
				    material1 == Materials::T::ISOTROPIC &&
				    material2 == Materials::T::ISOTROPIC) {
					
					return std::make_shared<ConcreteContactCorrector<
							AcousticModelD, IsotropicMaterial,
							AcousticModelD, IsotropicMaterial, TGrid,
							SlideContactMatrixCreator<AcousticModelD, AcousticModelD>>>();
					
				} else {
					THROW_UNSUPPORTED("Incompatible or unsupported contact conditions, \
							models and materials combination");
				}
				
			default:
				THROW_INVALID_ARG("Unknown type of contact condition");
		}
	}
	
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_CONTACTCORRECTOR_HPP
