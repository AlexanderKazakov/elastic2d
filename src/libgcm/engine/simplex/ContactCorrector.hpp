#ifndef LIBGCM_SIMPLEX_CONTACTCORRECTOR_HPP
#define LIBGCM_SIMPLEX_CONTACTCORRECTOR_HPP

#include <list>

#include <libgcm/engine/simplex/DefaultMesh.hpp>
#include <libgcm/engine/simplex/common.hpp>
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
	typedef typename TGrid::MatrixDD    MatrixDD;
	
	struct NodesContact {
		/// pair of nodes in contact (their iterators in grids)
		Iterator first, second;
		/// contact normal (the direction is from first to second)
		RealD normal;
	};
	
	/**
	 * Apply contact correction for all node pairs from the list
	 * along the contact normal direction.
	 * It's supposed that gcm-matrices in contact nodes are written
	 * in local basis and the first direction calculation (stage 0)
	 * performed along contact normal direction.
	 * Thus, this correction must be called after first stage only,
	 * because other directions are degenerate a priori.
	 */
	virtual void applyInLocalBasis(
			std::shared_ptr<AbstractMesh<TGrid>> a,
			std::shared_ptr<AbstractMesh<TGrid>> b,
			std::list<NodesContact> nodesInContact) const = 0;
	
	/**
	 * Apply contact correction for all node pairs from the list
	 * along the direction of the given stage.
	 * It's supposed that gcm-matrices in contact nodes are written
	 * in global basis as in inner nodes.
	 * Thus, this correction must be called after all stages.
	 * @param nextPdeLayerIndex is equal to stage if 
	 * splitting by directions is performed by summ, but if
	 * splitting by directions is done by product, it equals to 0 on all stages
	 */
	virtual void applyInGlobalBasis(
			const int nextPdeLayerIndex,
			const int stage,
			std::shared_ptr<AbstractMesh<TGrid>> a,
			std::shared_ptr<AbstractMesh<TGrid>> b,
			std::list<NodesContact> nodesInContact) const = 0;
	
	
protected:
	template<typename PdeVector>
	struct CorrectionResult {
		bool isSuccessful;
		PdeVector valueA, valueB;
	};
	
	/**
	 * General expression of linear contact condition is:
	 *     B1_A * u_A = B1_B * u_B,
	 *     B2_A * u_A = B2_B * u_B,
	 * where u is pde-vector and B is border matrices.
	 * Given with inner-calculated pde vectors, we correct them
	 * with outer waves combination (Omega) in order to satisfy contact condition.
	 * But this is not always possible due to degeneracies
	 * @see BorderCorrector
	 */
	template<typename PdeVector, typename MatrixOmega, typename MatrixB>
	CorrectionResult<PdeVector>
	calculateOuterWaveCorrection(
			const PdeVector& uA,
			const MatrixOmega& OmegaA, const MatrixB& B1A, const MatrixB& B2A,
			const PdeVector& uB,
			const MatrixOmega& OmegaB, const MatrixB& B1B, const MatrixB& B2B) const {
		const auto R1 = B1A * OmegaA;
//		if (linal::determinant(R1) == 0) {
		if (std::fabs(linal::determinant(R1)) < 0.1) {
			return { false, PdeVector::Zeros(), PdeVector::Zeros() };
		}
		const auto R = linal::invert(R1);
		const auto p = R * (B1B * uB - B1A * uA);
		const auto Q = R * (B1B * OmegaB);
		
		const auto A = (B2B * OmegaB) - ((B2A * OmegaA) * Q);
		const auto f = ((B2A * OmegaA) * p) + (B2A * uA) - (B2B * uB);
//		if (linal::determinant(A) == 0) {
		if (std::fabs(linal::determinant(A)) < 0.1) {
			return { false, PdeVector::Zeros(), PdeVector::Zeros() };
		}
		
		const auto alphaB = linal::solveLinearSystem(A, f);
		const auto alphaA = p + Q * alphaB;
		return { true, OmegaA * alphaA, OmegaB * alphaB };
	}
	
	USE_AND_INIT_LOGGER("gcm.simplex.ContactCorrector")
};



template<typename ModelA, typename MaterialA,
         typename ModelB, typename MaterialB,
         typename TGrid,
         typename ContactMatrixCreator>
class ContactCorrectorInPdeVectors : public AbstractContactCorrector<TGrid> {
public:
	typedef DefaultMesh<ModelA, TGrid, MaterialA> MeshA;
	typedef DefaultMesh<ModelB, TGrid, MaterialB> MeshB;
	typedef typename MeshA::PdeVector            PdeVector;
	typedef typename MeshA::PdeVariables         PdeVariables;
	typedef typename MeshA::WaveIndices          WaveIndices;
	static const int DIMENSIONALITY = MeshA::DIMENSIONALITY;
	static const int OUTER_NUMBER = ModelA::OUTER_NUMBER;
	
	typedef AbstractContactCorrector<TGrid> Base;
	typedef typename Base::NodesContact     NodesContact;
	typedef typename Base::RealD            RealD;
	typedef typename Base::MatrixDD         MatrixDD;
	
	ContactCorrectorInPdeVectors(const ContactConditions::T condition) :
			_condition(condition) { }
	
	virtual void applyInLocalBasis(
			std::shared_ptr<AbstractMesh<TGrid>> a,
			std::shared_ptr<AbstractMesh<TGrid>> b,
			std::list<NodesContact> nodesInContact) const override {
		const int stage = 0; ///< the only valid stage number
		std::shared_ptr<MeshA> meshA = std::dynamic_pointer_cast<MeshA>(a);
		assert_true(meshA);
		std::shared_ptr<MeshB> meshB = std::dynamic_pointer_cast<MeshB>(b);
		assert_true(meshB);
		
		for (const NodesContact& nodesContact : nodesInContact) {
			const auto OmegaA = getColumnsFromGcmMatrices<ModelA>(
					stage, ModelA::RIGHT_INVARIANTS, meshA->matrices(nodesContact.first));
			const auto OmegaB = getColumnsFromGcmMatrices<ModelB>(
					stage, ModelB::RIGHT_INVARIANTS, meshB->matrices(nodesContact.second));
			const auto B1A = ContactMatrixCreator::createB1A(nodesContact.normal);
			const auto B1B = ContactMatrixCreator::createB1B(nodesContact.normal);
			const auto B2A = ContactMatrixCreator::createB2A(nodesContact.normal);
			const auto B2B = ContactMatrixCreator::createB2B(nodesContact.normal);
			
			PdeVector& uA = meshA->_pdeNew(stage, nodesContact.first);
			PdeVector& uB = meshB->_pdeNew(stage, nodesContact.second);
			
			const auto correction = this->calculateOuterWaveCorrection(
					uA, OmegaA, B1A, B2A,
					uB, OmegaB, B1B, B2B);
			assert_true(correction.isSuccessful);
			uA += correction.valueA;
			uB += correction.valueB;
		}
	}
	
	virtual void applyInGlobalBasis(
			const int nextPdeLayerIndex,
			const int stage,
			std::shared_ptr<AbstractMesh<TGrid>> a,
			std::shared_ptr<AbstractMesh<TGrid>> b,
			std::list<NodesContact> nodesInContact) const override {
		std::shared_ptr<MeshA> meshA = std::dynamic_pointer_cast<MeshA>(a);
		assert_true(meshA);
		std::shared_ptr<MeshB> meshB = std::dynamic_pointer_cast<MeshB>(b);
		assert_true(meshB);
		
		for (const NodesContact& nodesContact : nodesInContact) {
			const WaveIndices outersA = meshA->waveIndices(nodesContact.first);
			const WaveIndices outersB = meshB->waveIndices(nodesContact.second);
			const auto B1A = ContactMatrixCreator::createB1A(nodesContact.normal);
			const auto B1B = ContactMatrixCreator::createB1B(nodesContact.normal);
			const auto B2A = ContactMatrixCreator::createB2A(nodesContact.normal);
			const auto B2B = ContactMatrixCreator::createB2B(nodesContact.normal);
			PdeVariables& uA = meshA->_pdeVarsNew(nextPdeLayerIndex, nodesContact.first);
			PdeVariables& uB = meshB->_pdeVarsNew(nextPdeLayerIndex, nodesContact.second);
			
			if ((outersA == ModelA::RIGHT_INVARIANTS && outersB == ModelB::LEFT_INVARIANTS) ||
				(outersB == ModelB::RIGHT_INVARIANTS && outersA == ModelA::LEFT_INVARIANTS)) {
			/// Normal case for contact corrector
				const auto OmegaA = getColumnsFromGcmMatrices<ModelA>(
						stage, outersA, meshA->matrices(nodesContact.first));
				const auto OmegaB = getColumnsFromGcmMatrices<ModelB>(
						stage, outersB, meshB->matrices(nodesContact.second));
				const auto correction = this->calculateOuterWaveCorrection(
						uA, OmegaA, B1A, B2A,
						uB, OmegaB, B1B, B2B);
				if (correction.isSuccessful) {
					uA += correction.valueA;
					uB += correction.valueB;
				} else {
					ModelA::applyPlainContactCorrection(
							uA, uB, _condition, nodesContact.normal);
				}
				
			} else if ((outersA.size() == 2 * OUTER_NUMBER && outersB.empty()) ||
					(outersB.size() == 2 * OUTER_NUMBER && outersA.empty())) {
			/// Double-outer case: apply correction as average from both sides
				// TODO - here should be double-outer border case,
				// but it requires to solve 6x6 SLE
				ModelA::applyPlainContactCorrection(
						uA, uB, _condition, nodesContact.normal);
				
			} else {
//				THROW_UNSUPPORTED("TODO");
				ModelA::applyPlainContactCorrection(
						uA, uB, _condition, nodesContact.normal);
				
			}
		}
	}
	
	
private:
	const ContactConditions::T _condition;
};



template<typename ModelA, typename MaterialA,
         typename ModelB, typename MaterialB,
         typename TGrid,
         typename ContactMatrixCreator>
class ContactCorrectorInRiemannInvariants : public AbstractContactCorrector<TGrid> {
public:
	typedef DefaultMesh<ModelA, TGrid, MaterialA> MeshA;
	typedef DefaultMesh<ModelB, TGrid, MaterialB> MeshB;
	typedef typename MeshA::GCM_MATRICES          GcmMatrices;
	typedef typename MeshA::PdeVector             PdeVector;
	typedef typename MeshA::PdeVariables          PdeVariables;
	typedef typename MeshA::WaveIndices           WaveIndices;
	static const int DIMENSIONALITY = MeshA::DIMENSIONALITY;
	static const int OUTER_NUMBER = ModelA::OUTER_NUMBER;
	
	typedef AbstractContactCorrector<TGrid> Base;
	typedef typename Base::NodesContact     NodesContact;
	typedef typename Base::RealD            RealD;
	typedef typename Base::MatrixDD         MatrixDD;
	
	ContactCorrectorInRiemannInvariants(const ContactConditions::T condition) :
			pdeCorrector(condition) { }
	
	virtual void applyInLocalBasis(
			std::shared_ptr<AbstractMesh<TGrid>> a,
			std::shared_ptr<AbstractMesh<TGrid>> b,
			std::list<NodesContact> nodesInContact) const override {
		const int stage = 0; ///< the only valid stage number
		std::shared_ptr<MeshA> meshA = std::dynamic_pointer_cast<MeshA>(a);
		assert_true(meshA);
		std::shared_ptr<MeshB> meshB = std::dynamic_pointer_cast<MeshB>(b);
		assert_true(meshB);
		
		convertToPdeVariables(stage, stage, meshA, meshB, nodesInContact);
		pdeCorrector.applyInLocalBasis(a, b, nodesInContact);
		convertToRiemannInvariants(stage, stage, meshA, meshB, nodesInContact);
	}
	
	virtual void applyInGlobalBasis(
			const int nextPdeLayerIndex,
			const int stage,
			std::shared_ptr<AbstractMesh<TGrid>> a,
			std::shared_ptr<AbstractMesh<TGrid>> b,
			std::list<NodesContact> nodesInContact) const override {
		std::shared_ptr<MeshA> meshA = std::dynamic_pointer_cast<MeshA>(a);
		assert_true(meshA);
		std::shared_ptr<MeshB> meshB = std::dynamic_pointer_cast<MeshB>(b);
		assert_true(meshB);
		
		convertToPdeVariables(nextPdeLayerIndex, stage, meshA, meshB, nodesInContact);
		pdeCorrector.applyInGlobalBasis(nextPdeLayerIndex, stage, a, b, nodesInContact);
		convertToRiemannInvariants(nextPdeLayerIndex, stage, meshA, meshB, nodesInContact);
	}
	
	
private:
	const ContactCorrectorInPdeVectors<ModelA, MaterialA, ModelB, MaterialB,
			TGrid, ContactMatrixCreator> pdeCorrector;
	
	
	void convertToPdeVariables(const int nextPdeLayerIndex, const int stage,
			std::shared_ptr<MeshA>& meshA, std::shared_ptr<MeshB>& meshB,
			std::list<NodesContact>& nodesInContact) const {
		for (const NodesContact& nodesContact : nodesInContact) {
			PdeVector& uA = meshA->_pdeNew(nextPdeLayerIndex, nodesContact.first);
			const GcmMatrices& gcmMatricesA = *(meshA->matrices(nodesContact.first));
			uA = gcmMatricesA(stage).U1 * uA;
			PdeVector& uB = meshB->_pdeNew(nextPdeLayerIndex, nodesContact.second);
			const GcmMatrices& gcmMatricesB = *(meshB->matrices(nodesContact.second));
			uB = gcmMatricesB(stage).U1 * uB;
		}
	}
	
	void convertToRiemannInvariants(const int nextPdeLayerIndex, const int stage,
			std::shared_ptr<MeshA>& meshA, std::shared_ptr<MeshB>& meshB,
			std::list<NodesContact>& nodesInContact) const {
		for (const NodesContact& nodesContact : nodesInContact) {
			PdeVector& uA = meshA->_pdeNew(nextPdeLayerIndex, nodesContact.first);
			const GcmMatrices& gcmMatricesA = *(meshA->matrices(nodesContact.first));
			uA = gcmMatricesA(stage).U * uA;
			PdeVector& uB = meshB->_pdeNew(nextPdeLayerIndex, nodesContact.second);
			const GcmMatrices& gcmMatricesB = *(meshB->matrices(nodesContact.second));
			uB = gcmMatricesB(stage).U * uB;
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
			const GcmType gcmType,
			const ContactConditions::T condition,
			const Models::T model1, const Materials::T material1,
			const Models::T model2, const Materials::T material2) {
		
		// TODO - this is a horrible shit:
		switch (gcmType) {
		
		case GcmType::ADVECT_RIEMANN_INVARIANTS:
		switch (condition) {
			case ContactConditions::T::ADHESION:
				if (model1 == Models::T::ELASTIC &&
				    model2 == Models::T::ELASTIC &&
				    material1 == Materials::T::ISOTROPIC &&
				    material2 == Materials::T::ISOTROPIC) {
					
					return std::make_shared<ContactCorrectorInRiemannInvariants<
							ElasticModelD, IsotropicMaterial,
							ElasticModelD, IsotropicMaterial, TGrid,
							AdhesionContactMatrixCreator<ElasticModelD, ElasticModelD>>>(condition);
					
				} else {
					THROW_UNSUPPORTED("Incompatible or unsupported contact conditions, \
							models and materials combination");
				}
				
			case ContactConditions::T::SLIDE:
				if (model1 == Models::T::ACOUSTIC &&
				    model2 == Models::T::ACOUSTIC &&
				    material1 == Materials::T::ISOTROPIC &&
				    material2 == Materials::T::ISOTROPIC) {
					
					return std::make_shared<ContactCorrectorInRiemannInvariants<
							AcousticModelD, IsotropicMaterial,
							AcousticModelD, IsotropicMaterial, TGrid,
							SlideContactMatrixCreator<AcousticModelD, AcousticModelD>>>(condition);
					
				} else {
					THROW_UNSUPPORTED("Incompatible or unsupported contact conditions, \
							models and materials combination");
				}
				
			default:
				THROW_INVALID_ARG("Unknown type of contact condition");
		}
		
		case GcmType::ADVECT_PDE_VECTORS:
		switch (condition) {
			case ContactConditions::T::ADHESION:
				if (model1 == Models::T::ELASTIC &&
					model2 == Models::T::ELASTIC &&
					material1 == Materials::T::ISOTROPIC &&
					material2 == Materials::T::ISOTROPIC) {
					
					return std::make_shared<ContactCorrectorInPdeVectors<
							ElasticModelD, IsotropicMaterial,
							ElasticModelD, IsotropicMaterial, TGrid,
							AdhesionContactMatrixCreator<ElasticModelD, ElasticModelD>>>(condition);
					
				} else {
					THROW_UNSUPPORTED("Incompatible or unsupported contact conditions, \
							models and materials combination");
				}
				
			case ContactConditions::T::SLIDE:
				if (model1 == Models::T::ACOUSTIC &&
					model2 == Models::T::ACOUSTIC &&
					material1 == Materials::T::ISOTROPIC &&
					material2 == Materials::T::ISOTROPIC) {
					
					return std::make_shared<ContactCorrectorInPdeVectors<
							AcousticModelD, IsotropicMaterial,
							AcousticModelD, IsotropicMaterial, TGrid,
							SlideContactMatrixCreator<AcousticModelD, AcousticModelD>>>(condition);
					
				} else {
					THROW_UNSUPPORTED("Incompatible or unsupported contact conditions, \
							models and materials combination");
				}
				
			default:
				THROW_INVALID_ARG("Unknown type of contact condition");
		}
		
		default:
		THROW_UNSUPPORTED("Unknown type of gcm-method");
		
		}
	}
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_CONTACTCORRECTOR_HPP
