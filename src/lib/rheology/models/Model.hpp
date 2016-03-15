#ifndef LIBGCM_MODEL_HPP
#define LIBGCM_MODEL_HPP

#include <lib/rheology/materials/materials.hpp>
#include <lib/rheology/variables/variables.hpp>
#include <lib/rheology/ode/ode.hpp>
#include <lib/rheology/correctors/correctors.hpp>
#include <lib/numeric/gcm/GcmMatrices.hpp>


namespace gcm {

	/**
	 * Map between wave types and corresponding columns in matrix U1 in GcmMatrix.
	 */
	typedef std::map<Waves::T, int> WavesEigenvectorsMap;
	/**
	 * Map between material types and corresponding WavesEigenvectorsMap for concrete Model.
	 * Note that in concrete Model for concrete Material the order of eigenvalues
	 * in GcmMatrices has to be the same for all spatial directions.
	 */
	typedef std::map<Materials::T, WavesEigenvectorsMap> MaterialsWavesMap;

	/**
	 * The list of rheology models
	 */

	class Elastic1DModel {
	public:
		static const int DIMENSIONALITY = 1;

		typedef VelocitySigmaVariables<DIMENSIONALITY>   PdeVariables;
		typedef DummyOde                                 InternalOde;
		typedef DummyCorrector                           Corrector;

		static const int PDE_SIZE = PdeVariables::SIZE;

		typedef linal::Vector<PDE_SIZE, PdeVariables>    PdeVector;
		typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY>    GCM_MATRICES;
		typedef typename InternalOde::Variables          OdeVariables;
		typedef std::shared_ptr<GCM_MATRICES>            GcmMatricesPtr;
		typedef std::shared_ptr<const GCM_MATRICES>      ConstGcmMatricesPtr;

		static const MaterialsWavesMap MATERIALS_WAVES_MAP;

		// template<Node>
		static void constructGcmMatrices
		(GcmMatricesPtr m, const PdeVector& pde, std::shared_ptr<const IsotropicMaterial> material);
	};

	class Elastic2DModel {
	public:
		static const int DIMENSIONALITY = 2;

		typedef VelocitySigmaVariables<DIMENSIONALITY>   PdeVariables;
		typedef DummyOde                                 InternalOde;
		typedef DummyCorrector                           Corrector;

		static const int PDE_SIZE = PdeVariables::SIZE;

		typedef linal::Vector<PDE_SIZE, PdeVariables>    PdeVector;
		typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY>    GCM_MATRICES;
		typedef typename InternalOde::Variables          OdeVariables;
		typedef std::shared_ptr<GCM_MATRICES>            GcmMatricesPtr;
		typedef std::shared_ptr<const GCM_MATRICES>      ConstGcmMatricesPtr;

		static const MaterialsWavesMap MATERIALS_WAVES_MAP;

		static void constructGcmMatrices
		(GcmMatricesPtr m, const PdeVector& pde, std::shared_ptr<const IsotropicMaterial> material);
	};

	class Elastic3DModel {
	public:
		static const int DIMENSIONALITY = 3;

		typedef VelocitySigmaVariables<DIMENSIONALITY>   PdeVariables;
		typedef DummyOde                                 InternalOde;
		typedef DummyCorrector                           Corrector;

		static const int PDE_SIZE = PdeVariables::SIZE;

		typedef linal::Vector<PDE_SIZE, PdeVariables>    PdeVector;
		typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY>    GCM_MATRICES;
		typedef typename InternalOde::Variables          OdeVariables;
		typedef std::shared_ptr<GCM_MATRICES>            GcmMatricesPtr;
		typedef std::shared_ptr<const GCM_MATRICES>      ConstGcmMatricesPtr;

		static const MaterialsWavesMap MATERIALS_WAVES_MAP;

		static void constructGcmMatrices
		(GcmMatricesPtr m, const PdeVector& pde, std::shared_ptr<const IsotropicMaterial> material);
		static void constructGcmMatrices
		(GcmMatricesPtr m, const PdeVector& pde, std::shared_ptr<const OrthotropicMaterial> material);
	};

	class SuperDuperModel {
	public:
		static const int DIMENSIONALITY = 3;

		typedef VelocitySigmaVariables<DIMENSIONALITY>   PdeVariables;
		typedef ContinualDamageOde                       InternalOde;
		typedef IdealPlasticFlowCorrector                Corrector;

		static const int PDE_SIZE = PdeVariables::SIZE;

		typedef linal::Vector<PDE_SIZE, PdeVariables>    PdeVector;
		typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY>    GCM_MATRICES;
		typedef typename InternalOde::Variables          OdeVariables;
		typedef std::shared_ptr<GCM_MATRICES>            GcmMatricesPtr;
		typedef std::shared_ptr<const GCM_MATRICES>      ConstGcmMatricesPtr;

		static const MaterialsWavesMap MATERIALS_WAVES_MAP;

		template<typename TMaterialPtr>
		static void constructGcmMatrices
		(GcmMatricesPtr m, const PdeVector& pde, const TMaterialPtr material) {
			Elastic3DModel::constructGcmMatrices(m, pde, material);
		}

	};

}

#endif //LIBGCM_MODEL_HPP
