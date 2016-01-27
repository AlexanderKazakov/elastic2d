#ifndef LIBGCM_MODEL_HPP
#define LIBGCM_MODEL_HPP

#include <memory>
#include <mpi.h>

#include <lib/rheology/variables/VelocitySigmaVariables.hpp>
#include <lib/rheology/materials/IsotropicMaterial.hpp>

namespace gcm {

	/**
	 * Generator for rheology models
	 * @tparam Components building blocks to construct some model
	 */
	template<typename... Components>
	struct Model : Components... {

	};


	/**
	 * Using these envelopes as template arguments of Model we can construct any necessary Model.
	 */

	template<typename TVector>
	struct VectorEnvelope {
		typedef TVector Vector;
		typedef typename Vector::ContainerType Variables;
		static const int M = Vector::M;

		TVector u;
	};

	template<typename TGcmMatrices>
	struct GcmMatricesPtrEnvelope {
		typedef TGcmMatrices GcmMatrices;
		typedef typename GcmMatrices::Matrix Matrix;
		static const int DIMENSIONALITY = GcmMatrices::DIMENSIONALITY;
		std::shared_ptr<GcmMatrices> matrix;
	};

	template<typename TCoords>
	struct CoordsEnvelope {
		TCoords coords;
	};


	typedef Model<VectorEnvelope<linal::Vector<2, VelocitySigmaVariables<1>>>,
			GcmMatricesPtrEnvelope<GcmMatrices<VelocitySigmaVariables<1>, IsotropicMaterial>>> IdealElastic1DModel;
	typedef Model<VectorEnvelope<linal::Vector<5, VelocitySigmaVariables<2>>>,
			GcmMatricesPtrEnvelope<GcmMatrices<VelocitySigmaVariables<2>, IsotropicMaterial>>> IdealElastic2DModel;
	typedef Model<VectorEnvelope<linal::Vector<9, VelocitySigmaVariables<3>>>,
			GcmMatricesPtrEnvelope<GcmMatrices<VelocitySigmaVariables<3>, IsotropicMaterial>>> IdealElastic3DModel;


}

#endif //LIBGCM_MODEL_HPP
