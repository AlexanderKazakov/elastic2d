#ifndef LIBGCM_MODEL_HPP
#define LIBGCM_MODEL_HPP

#include <lib/rheology/materials/Materials.hpp>
#include <lib/rheology/variables/Variables.hpp>

namespace gcm {

	/**
	 * Covers rheology aspects
	 */
	template<typename TVariables, typename TMaterial>
	struct Model {
		typedef TVariables Variables;
		typedef TMaterial  Material;

		static const int DIMENSIONALITY = Variables::DIMENSIONALITY;
		static const int PDE_SIZE = Variables::SIZE;

		typedef linal::Vector<PDE_SIZE, Variables> Vector;
		typedef GcmMatrices<Variables, Material> GCM_MATRICES;
	};

	// type definitions

	typedef Model<VelocitySigmaVariables<1>, IsotropicMaterial> Elastic1DModel;

	typedef Model<VelocitySigmaVariables<2>, IsotropicMaterial> Elastic2DModel;

	typedef Model<VelocitySigmaVariables<3>, IsotropicMaterial> Elastic3DModel;
	typedef Model<VelocitySigmaVariables<3>, OrthotropicMaterial> OrthotropicElastic3DModel;
}

#endif //LIBGCM_MODEL_HPP
