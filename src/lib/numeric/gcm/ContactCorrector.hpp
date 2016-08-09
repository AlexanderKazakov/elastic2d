#ifndef LIBGCM_CONTACTCORRECTOR_HPP
#define LIBGCM_CONTACTCORRECTOR_HPP

#include <list>
#include <lib/mesh/grid/AbstractGrid.hpp>

#include <lib/mesh/DefaultMesh.hpp>
#include <lib/rheology/models/ElasticModel.hpp>
#include <lib/rheology/models/AcousticModel.hpp>
#include <lib/rheology/materials/IsotropicMaterial.hpp>


namespace gcm {

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
	 */
	virtual void apply(AbstractGrid* a, AbstractGrid* b,
			std::list<NodesContact> nodesInContact) = 0;
	
	
protected:
	/// maximal found condition numbers of contact corrector matrices
	real maxConditionR = 0;
	real maxConditionA = 0;
	
	
	/**
	 * General expression of linear contact condition is:
	 *     B1_A * u_A = B1_B * u_B,
	 *     B2_A * u_A = B2_B * u_B,
	 * where u is pde-vector and B is border matrices.
	 * Given with inner-calculated pde vectors, we correct them
	 * with outer waves combination in order to satisfy contact conditions.
	 * @see BorderCondition
	 */
	template<typename PdeVector, typename MatrixOmega, typename MatrixB>
	void
	correctNodesContact(
			PdeVector& uA,
			const MatrixOmega& omegaA, const MatrixB& B1A, const MatrixB& B2A,
			PdeVector& uB,
			const MatrixOmega& omegaB, const MatrixB& B1B, const MatrixB& B2B) {
		
		const auto R = linal::invert(B1A * omegaA);
		const real currentConditionR = linal::conditionNumber(R);
		if (currentConditionR > maxConditionR) {
			maxConditionR = currentConditionR;
			LOG_INFO("New maximal condition number in matrix R: " << maxConditionR);
		}
		
		const auto p = R * (B1B * uB - B1A * uA);
		const auto Q = R * (B1B * omegaB);
		
		const auto A = (B2B * omegaB) - ((B2A * omegaA) * Q);
		const real currentConditionA = linal::conditionNumber(A);
		if (currentConditionA > maxConditionA) {
			maxConditionA = currentConditionA;
			LOG_INFO("New maximal condition number in matrix A: " << maxConditionA);
		}
		
		const auto f = ((B2A * omegaA) * p) + (B2A * uA) - (B2B * uB);
		
		const auto alphaB = linal::solveLinearSystem(A, f);
		const auto alphaA = p + Q * alphaB;
		
		uA += omegaA * alphaA;
		uB += omegaB * alphaB;
	}
	
	
	USE_AND_INIT_LOGGER("gcm.ContactCorrector")
};



template<typename ModelA, typename MaterialA,
         typename ModelB, typename MaterialB,
         typename TGrid>
class AdhesionContactCorrector : public AbstractContactCorrector<TGrid> {
public:
	
	typedef DefaultMesh<ModelA, TGrid, MaterialA> MeshA;
	typedef DefaultMesh<ModelA, TGrid, MaterialB> MeshB;
	
	typedef AbstractContactCorrector<TGrid> Base;
	typedef typename Base::NodesContact     NodesContact;
	
	virtual void apply(AbstractGrid* a, AbstractGrid* b, 
			std::list<NodesContact> nodesInContact) override {
		
		MeshA* meshA = dynamic_cast<MeshA*>(a);
		assert_true(meshA);
		MeshB* meshB = dynamic_cast<MeshB*>(b);
		assert_true(meshB);
		
		for (const NodesContact& nodesContact : nodesInContact) {
			
			const auto omegaA = ModelA::constructOuterEigenvectors(
					meshA->material(nodesContact.first),
					linal::createLocalBasis(  nodesContact.normal));
			const auto omegaB = ModelB::constructOuterEigenvectors(
					meshB->material(nodesContact.second),
					linal::createLocalBasis( -nodesContact.normal));
			
			const auto B1A = ModelA::borderMatrixFixedVelocityGlobalBasis(
					nodesContact.normal);
			const auto B1B = ModelB::borderMatrixFixedVelocityGlobalBasis(
					nodesContact.normal);
			
			const auto B2A = ModelA::borderMatrixFixedForceGlobalBasis(
					nodesContact.normal);
			const auto B2B = ModelB::borderMatrixFixedForceGlobalBasis(
					nodesContact.normal);
			
			auto& uA = meshA->_pdeNew(nodesContact.first);
			auto& uB = meshB->_pdeNew(nodesContact.second);
			
			this->correctNodesContact(uA, omegaA, B1A, B2A,
			                          uB, omegaB, B1B, B2B);
		}
	}
	
};



template<typename ModelA, typename MaterialA,
         typename ModelB, typename MaterialB,
         typename TGrid>
class SlidingContactCorrector : public AbstractContactCorrector<TGrid> {
public:
	
	typedef DefaultMesh<ModelA, TGrid, MaterialA> MeshA;
	typedef DefaultMesh<ModelA, TGrid, MaterialB> MeshB;
	
	typedef AbstractContactCorrector<TGrid> Base;
	typedef typename Base::NodesContact     NodesContact;
	
	virtual void apply(AbstractGrid* a, AbstractGrid* b, 
			std::list<NodesContact> nodesInContact) override {
		
		MeshA* meshA = dynamic_cast<MeshA*>(a);
		assert_true(meshA);
		MeshB* meshB = dynamic_cast<MeshB*>(b);
		assert_true(meshB);
		
		for (const NodesContact& nodesContact : nodesInContact) {
			
			const auto omegaA = ModelA::constructOuterEigenvectors(
					meshA->material(nodesContact.first),
					linal::createLocalBasis(  nodesContact.normal));
			const auto omegaB = ModelB::constructOuterEigenvectors(
					meshB->material(nodesContact.second),
					linal::createLocalBasis( -nodesContact.normal));
			
			const auto B1A = ModelA::borderMatrixFixedVelocity(
					nodesContact.normal);
			const auto B1B = ModelB::borderMatrixFixedVelocity(
					nodesContact.normal);
			
			const auto B2A = ModelA::borderMatrixFixedForce(
					nodesContact.normal);
			const auto B2B = ModelB::borderMatrixFixedForce(
					nodesContact.normal);
			
			auto& uA = meshA->_pdeNew(nodesContact.first);
			auto& uB = meshB->_pdeNew(nodesContact.second);
			
			this->correctNodesContact(uA, omegaA, B1A, B2A,
			                          uB, omegaB, B1B, B2B);
		}
	}
	
};



template<typename TGrid>
class ContactCorrectorFactory {
public:
	
	static const int DIMENSIONALITY = TGrid::DIMENSIONALITY;
	
	typedef ElasticModel<DIMENSIONALITY>     ElasticModelD;
	typedef AcousticModel<DIMENSIONALITY>    AcousticModelD;
	
	static std::shared_ptr<AbstractContactCorrector<TGrid>> create(
			const ContactConditions::T type,
			const Models::T model1, const Materials::T material1,
			const Models::T model2, const Materials::T material2) {
		
		switch (type) {
			case ContactConditions::T::ADHESION:
				if (model1 == Models::T::ELASTIC &&
				    model2 == Models::T::ELASTIC &&
				    material1 == Materials::T::ISOTROPIC &&
				    material2 == Materials::T::ISOTROPIC) {
					
					return std::make_shared<AdhesionContactCorrector<
							ElasticModelD, IsotropicMaterial,
							ElasticModelD, IsotropicMaterial, TGrid>>();
					
				} else {
					THROW_UNSUPPORTED("Incompatible or unsupported contact conditions, \
							models and materials combination");
				}
				
			case ContactConditions::T::SLIDE:
				if (model1 == Models::T::ACOUSTIC &&
				    model2 == Models::T::ACOUSTIC &&
				    material1 == Materials::T::ISOTROPIC &&
				    material2 == Materials::T::ISOTROPIC) {
					
					return std::make_shared<SlidingContactCorrector<
							AcousticModelD, IsotropicMaterial,
							AcousticModelD, IsotropicMaterial, TGrid>>();
					
				} else {
					THROW_UNSUPPORTED("Incompatible or unsupported contact conditions, \
							models and materials combination");
				}
				
			default:
				THROW_INVALID_ARG("Unknown type of contact condition");
		}
	}
	
	
};


}


#endif // LIBGCM_CONTACTCORRECTOR_HPP
