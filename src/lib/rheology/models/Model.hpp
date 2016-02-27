#ifndef LIBGCM_MODEL_HPP
#define LIBGCM_MODEL_HPP

#include <lib/rheology/materials/materials.hpp>
#include <lib/rheology/variables/variables.hpp>
#include <lib/rheology/ode/ode.hpp>
#include <lib/rheology/correctors/correctors.hpp>

namespace gcm {

	template<
			typename TVariables, 
			typename TMaterial, 
			typename TInternalOde = DummyOde,
			typename TCorrector = DummyCorrector
			>
	struct Model {
		typedef TVariables           PdeVariables;
		typedef TMaterial            Material;
		typedef TInternalOde         InternalOde;
		typedef TCorrector           Corrector;

		static const int DIMENSIONALITY = PdeVariables::DIMENSIONALITY;
		static const int PDE_SIZE = PdeVariables::SIZE;

		typedef linal::Vector<PDE_SIZE, PdeVariables>         PdeVector;
		typedef GcmMatrices<PdeVariables, Material>           GCM_MATRICES;
		typedef typename InternalOde::Variables               OdeVariables;
	};


	typedef Model<VelocitySigmaVariables<1>, IsotropicMaterial> Elastic1DModel;

	typedef Model<VelocitySigmaVariables<2>, IsotropicMaterial> Elastic2DModel;
	typedef Model<VelocitySigmaVariables<2>, IsotropicMaterial, ContinualDamageOde> ContinualDamageElastic2DModel;
	typedef Model<VelocitySigmaVariables<2>, IsotropicMaterial, DummyOde, IdealPlasticFlowCorrector> IdealPlastic2DModel;

	typedef Model<VelocitySigmaVariables<3>, IsotropicMaterial> Elastic3DModel;
	typedef Model<VelocitySigmaVariables<3>, OrthotropicMaterial> OrthotropicElastic3DModel;

	typedef Model<VelocitySigmaVariables<3>, OrthotropicMaterial,
			ContinualDamageOde, IdealPlasticFlowCorrector> SuperDuperModel;
}

#endif //LIBGCM_MODEL_HPP
