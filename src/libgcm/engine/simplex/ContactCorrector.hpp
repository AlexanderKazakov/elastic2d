#ifndef LIBGCM_SIMPLEX_CONTACTCORRECTOR_HPP
#define LIBGCM_SIMPLEX_CONTACTCORRECTOR_HPP

#include <list>

#include <libgcm/engine/mesh/DefaultMesh.hpp>
#include <libgcm/rheology/models/ElasticModel.hpp>
#include <libgcm/rheology/models/AcousticModel.hpp>
#include <libgcm/rheology/materials/IsotropicMaterial.hpp>


namespace gcm {
namespace simplex {

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
	 * For all node pairs from the list, along given calculation direction,
	 * finalize calculation of contact nodes
	 * (which was begun in simplex::GridCharacteristicMethod)
	 */
	virtual void apply(
			const int stage,
			std::shared_ptr<AbstractMesh<TGrid>> a,
			std::shared_ptr<AbstractMesh<TGrid>> b,
			std::list<NodesContact> nodesInContact,
			const RealD& direction) = 0;
	
	USE_AND_INIT_LOGGER("gcm.simplex.ContactCorrector")
};


template<typename ModelA, typename MaterialA,
         typename ModelB, typename MaterialB,
         typename TGrid,
         typename RefractionCalculator>
class ConcreteContactCorrector : public AbstractContactCorrector<TGrid> {
public:
	typedef AbstractContactCorrector<TGrid>       Base;
	typedef typename Base::NodesContact           NodesContact;
	typedef typename Base::RealD                  RealD;
	typedef typename Base::Iterator               Iterator;
	
	typedef DefaultMesh<ModelA, TGrid, MaterialA> MeshA;
	typedef typename MeshA::WaveIndices           WaveIndicesA;
	typedef typename MeshA::PdeVector             PdeVectorA;
	typedef typename MeshA::Model::GcmMatrix      GcmMatrixA;
	
	typedef DefaultMesh<ModelB, TGrid, MaterialB> MeshB;
	typedef typename MeshB::WaveIndices           WaveIndicesB;
	typedef typename MeshB::PdeVector             PdeVectorB;
	typedef typename MeshB::Model::GcmMatrix      GcmMatrixB;
	
	virtual void apply(
			const int stage,
			std::shared_ptr<AbstractMesh<TGrid>> a,
			std::shared_ptr<AbstractMesh<TGrid>> b,
			std::list<NodesContact> nodesInContact,
			const RealD& calcDirection) override {
		std::shared_ptr<MeshA> meshA = std::dynamic_pointer_cast<MeshA>(a);
		assert_true(meshA);
		std::shared_ptr<MeshB> meshB = std::dynamic_pointer_cast<MeshB>(b);
		assert_true(meshB);
		
		for (const NodesContact& nodesContact : nodesInContact) {
			const GcmMatrixA& gcmMatrixA =
					(*(meshA->matrices(nodesContact.first)))(stage);
			const GcmMatrixB& gcmMatrixB =
					(*(meshB->matrices(nodesContact.second)))(stage);
			const std::pair<PdeVectorA, PdeVectorB> fromA =
				calcOneSide(nodesContact.first, nodesContact.second,
					*meshA, *meshB, gcmMatrixA, nodesContact.normal, calcDirection);
			const std::pair<PdeVectorB, PdeVectorA> fromB =
				calcOneSide(nodesContact.second, nodesContact.first,
					*meshB, *meshA, gcmMatrixB, -nodesContact.normal, calcDirection);
			
			meshA->_pdeNew(nodesContact.first) =
				gcmMatrixA.U1 * meshA->pdeNew(nodesContact.first) /* initial */ +
				fromA.first /* reflection */ +
				fromB.second /* refraction */;
			meshB->_pdeNew(nodesContact.second) =
				gcmMatrixB.U1 * meshB->pdeNew(nodesContact.second) /* initial */ +
				fromB.first /* reflection */ +
				fromA.second /* refraction */;
		}
	}
	
private:
	/// Calculate reflection wave in mesh1 (first) and
	/// refraction wave in mesh2 (second), caused by initial wave come from mesh1.
	/// Mesh{i} can be either MeshA or MeshB
	template<typename Mesh1, typename Mesh2>
	std::pair<typename Mesh1::PdeVector, typename Mesh2::PdeVector>
	calcOneSide(const Iterator one, const Iterator two,
			const Mesh1& mesh1, const Mesh2& mesh2, 
			const typename Mesh1::Model::GcmMatrix& gcmMatrix1,
			const RealD& normalOuterForMesh1, const RealD& calcDirection) {
		typedef typename Mesh1::PdeVector Pde1;
		typedef typename Mesh2::PdeVector Pde2;
		
		const Pde1 initialWavesInRiemannVariables = mesh1.pdeNew(one);
		const typename Mesh1::WaveIndices
				initialWavesIndicesInGcmMatrix1 = mesh1.waveIndices(one);
		
		Pde1 reflectionInPdeVariables = Pde1::Zeros();
		Pde2 refractionInPdeVariables = Pde2::Zeros();
		for (const int waveIndex : initialWavesIndicesInGcmMatrix1) {
			std::pair<Pde1, Pde2> reflectionAndRefraction =
				RefractionCalculator::calculate(
					initialWavesInRiemannVariables(waveIndex),
					waveIndex, normalOuterForMesh1, calcDirection, gcmMatrix1,
					mesh1.material(one), mesh2.material(two));
			reflectionInPdeVariables += reflectionAndRefraction.first;
			refractionInPdeVariables += reflectionAndRefraction.second;
		}
		
		return {reflectionInPdeVariables, refractionInPdeVariables};
	}
};


template<int Dimensionality>
struct AcousticToAcousticRefractionCalculator {
	typedef AcousticModel<Dimensionality> Model;
	typedef typename Model::RealD        RealD;
	typedef typename Model::GcmMatrix    GcmMatrix;
	typedef typename Model::Matrix       Matrix;
	typedef typename Model::PdeVector    PdeVector;
	
	/// Solve reflection-refraction problem
	/// the first item in answer is reflection wave,
	/// the second item -- refraction wave
	static std::pair<PdeVector, PdeVector> calculate(
			const real waveAmplitude, const int waveIndex,
			const RealD& contactNormal, const RealD& calcDirection,
			const GcmMatrix& matrixAlongCalcDirection,
			std::shared_ptr<const IsotropicMaterial> material1,
			std::shared_ptr<const IsotropicMaterial> material2) {
		const real rho1 = material1->rho;
		const real rho2 = material2->rho;
		const real c1 = Model::acousticVelocity(material1);
		const real c2 = Model::acousticVelocity(material2);
		
		const int invariantSign = Utils::sign(matrixAlongCalcDirection.L(waveIndex));
		const RealD initDirection = invariantSign * calcDirection;
		// cosinus of incident angle
		const real cos1 = linal::dotProduct(contactNormal, initDirection);
		assert_ge(cos1, 0);
		
		const RealD reflDirection =
				linal::reflectionDirection(contactNormal, initDirection);
		Matrix u1; //< helping temporary
		Model::constructEigenvectors(u1, material1,
				linal::createLocalBasis(invariantSign * reflDirection));
		const PdeVector reflectionWave = u1.getColumn(waveIndex);
		
		if (c1*c1 / (c2*c2) + cos1*cos1 - 1 >= 0) {
		/// normal refraction
			const RealD refrDirection = linal::refractionDirection(
					contactNormal, initDirection, c1, c2);
			const real cos2 = linal::dotProduct(contactNormal, refrDirection);
			Model::constructEigenvectors(u1, material2,
					linal::createLocalBasis(invariantSign * refrDirection));
			const PdeVector refractionWave = u1.getColumn(waveIndex);
			// reflection coefficient
			const real V = (c2 * rho2 * cos1 - c1 * rho1 * cos2) /
					(c2 * rho2 * cos1 + c1 * rho1 * cos2);
			// refraction coefficient
			const real W = 2 * rho1 * c2 * cos1 /
					(c2 * rho2 * cos1 + c1 * rho1 * cos2);
			return {V * waveAmplitude * reflectionWave,
					W * waveAmplitude * refractionWave};
			
		} else {
		/// total internal reflection
		/// FIXME (-waveAmplitude) is not correct
			return {-waveAmplitude * reflectionWave, PdeVector::Zeros()};
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
			const ContactConditions::T condition,
			const Models::T model1, const Materials::T material1,
			const Models::T model2, const Materials::T material2) {
		
		if (material1 != Materials::T::ISOTROPIC ||
				material2 != Materials::T::ISOTROPIC) {
			THROW_UNSUPPORTED("Unsupported material");
		}
		
		if (model1 != Models::T::ACOUSTIC ||
				model2 != Models::T::ACOUSTIC) {
			THROW_UNSUPPORTED("Unsupported model");
		}
		
		if (condition != ContactConditions::T::SLIDE) {
			THROW_UNSUPPORTED("Unsupported condition");
		}
		
		return std::make_shared<ConcreteContactCorrector<
				AcousticModelD, IsotropicMaterial,
				AcousticModelD, IsotropicMaterial, TGrid,
				AcousticToAcousticRefractionCalculator<DIMENSIONALITY>>>();
	}
};


} // namespace simplex
} // namespace gcm


#endif // LIBGCM_SIMPLEX_CONTACTCORRECTOR_HPP
