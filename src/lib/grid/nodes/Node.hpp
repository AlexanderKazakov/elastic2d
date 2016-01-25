#ifndef LIBGCM_NODE_HPP
#define LIBGCM_NODE_HPP

#include <memory>
#include <mpi.h>

#include <lib/rheology/variables/VelocitySigmaVariables.hpp>
#include <lib/rheology/materials/IsotropicMaterial.hpp>

namespace gcm {

	/**
	 * Generator for nodes
	 * @tparam Components building blocks to construct some node
	 */
	template<typename... Components>
	struct Node : Components... {
		static MPI::Datatype MPI_NODE_TYPE; // Special type for node for MPI connection
	};


	/**
	 * Using these envelopes as template arguments of Node we can construct any necessary Node.
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


	typedef Node<VectorEnvelope<linal::Vector<2, VelocitySigmaVariables<1>>>,
			GcmMatricesPtrEnvelope<GcmMatrices<VelocitySigmaVariables<1>, IsotropicMaterial>>> IdealElastic1DNode;
	typedef Node<VectorEnvelope<linal::Vector<5, VelocitySigmaVariables<2>>>,
			GcmMatricesPtrEnvelope<GcmMatrices<VelocitySigmaVariables<2>, IsotropicMaterial>>> IdealElastic2DNode;
	typedef Node<VectorEnvelope<linal::Vector<9, VelocitySigmaVariables<3>>>,
			GcmMatricesPtrEnvelope<GcmMatrices<VelocitySigmaVariables<3>, IsotropicMaterial>>> IdealElastic3DNode;


}

#endif //LIBGCM_NODE_HPP
