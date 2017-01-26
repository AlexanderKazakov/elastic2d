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
	 * For all nodes from the list, along given calculation direction,
	 * finalize calculation of border nodes
	 * (which was begun in simplex::GridCharacteristicMethod)
	 */
	virtual void apply(const int stage, std::shared_ptr<AbstractMesh<TGrid>> grid,
			std::list<NodeBorder> borderNodes, const RealD& calcDirection) = 0;
	
	USE_AND_INIT_LOGGER("gcm.simplex.BorderCorrector")
};


template<typename Model, typename Material, typename TGrid,
         typename BorderReflectionCalculator>
class ConcreteBorderCorrector : public AbstractBorderCorrector<TGrid> {
public:
	typedef DefaultMesh<Model, TGrid, Material> Mesh;
	typedef typename Mesh::WaveIndices          WaveIndices;
	typedef typename Mesh::PdeVector            PdeVector;
	typedef typename Mesh::Model::GcmMatrix     GcmMatrix;
	
	typedef AbstractBorderCorrector<TGrid>      Base;
	typedef typename Base::NodeBorder           NodeBorder;
	typedef typename Base::RealD                RealD;
	
	ConcreteBorderCorrector(const Task::BorderCondition& bc) :
			borderCondition(bc) { }
	
	virtual void apply(const int stage, std::shared_ptr<AbstractMesh<TGrid>> grid,
			std::list<NodeBorder> borderNodes, const RealD& calcDirection) override {
		std::shared_ptr<Mesh> mesh = std::dynamic_pointer_cast<Mesh>(grid);
		assert_true(mesh);
		// TODO - non-zero condition
//		const auto b = borderCondition.b();
		
		for (const NodeBorder& nodeBorder: borderNodes) {
			const PdeVector initialWavesInRiemannVariables =
					mesh->pdeNew(nodeBorder.iterator);
			const GcmMatrix& gcmMatrix =
					(*(mesh->matrices(nodeBorder.iterator)))(stage);
			const WaveIndices initialWavesIndicesInGcmMatrix =
					mesh->waveIndices(nodeBorder.iterator);
			
			PdeVector reflectionWavesInPdeVaiables = PdeVector::Zeros();
			for (const int waveIndex : initialWavesIndicesInGcmMatrix) {
				reflectionWavesInPdeVaiables +=
					BorderReflectionCalculator::calculate(
						initialWavesInRiemannVariables(waveIndex),
						waveIndex, nodeBorder.normal, calcDirection,
						gcmMatrix, mesh->material(nodeBorder.iterator));
			}
			
			const PdeVector initialWavesInPdeVariables =
					gcmMatrix.U1 * initialWavesInRiemannVariables;
			mesh->_pdeNew(nodeBorder.iterator) =
					initialWavesInPdeVariables + reflectionWavesInPdeVaiables;
		}
	}
	
	
private:
	const BorderCondition<Model> borderCondition;
};


// TODO replace this shit to templating by static member functions
template<typename TModel>
struct FixedForceBorderReflectionCalculator {
	typedef typename TModel::RealD        RealD;
	typedef typename TModel::GcmMatrix    GcmMatrix;
	typedef typename TModel::PdeVector    PdeVector;
	static PdeVector calculate(
			const real waveAmplitude, const int waveIndex,
			const RealD& borderNormal, const RealD& calcDirection,
			const GcmMatrix& matrixAlongCalcDirection,
			std::shared_ptr<const IsotropicMaterial> material) {
		return TModel::calculateReflectionZeroBorderForce(
				waveAmplitude, waveIndex, borderNormal, calcDirection,
				matrixAlongCalcDirection, material);
	}
};
template<typename TModel>
struct FixedVelocityBorderReflectionCalculator {
	typedef typename TModel::RealD        RealD;
	typedef typename TModel::GcmMatrix    GcmMatrix;
	typedef typename TModel::PdeVector    PdeVector;
	static PdeVector calculate(
			const real waveAmplitude, const int waveIndex,
			const RealD& borderNormal, const RealD& calcDirection,
			const GcmMatrix& matrixAlongCalcDirection,
			std::shared_ptr<const IsotropicMaterial> material) {
		return TModel::calculateReflectionZeroBorderVelocity(
				waveAmplitude, waveIndex, borderNormal, calcDirection,
				matrixAlongCalcDirection, material);
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
								FixedForceBorderReflectionCalculator<ElasticModelD>>>(
										condition);
					case Models::T::ACOUSTIC:
						return std::make_shared<ConcreteBorderCorrector<
								AcousticModelD, IsotropicMaterial, TGrid,
								FixedForceBorderReflectionCalculator<AcousticModelD>>>(
										condition);
					default:
						THROW_INVALID_ARG("Unknown type of model");
				}
				
			case BorderConditions::T::FIXED_VELOCITY:
				
				switch (model) {
					case Models::T::ELASTIC:
						return std::make_shared<ConcreteBorderCorrector<
								ElasticModelD, IsotropicMaterial, TGrid,
								FixedVelocityBorderReflectionCalculator<ElasticModelD>>>(
										condition);
					case Models::T::ACOUSTIC:
						return std::make_shared<ConcreteBorderCorrector<
								AcousticModelD, IsotropicMaterial, TGrid,
								FixedVelocityBorderReflectionCalculator<AcousticModelD>>>(
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
