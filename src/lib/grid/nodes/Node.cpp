#include <lib/grid/nodes/Node.hpp>
#include <lib/rheology/materials/IsotropicMaterial.hpp>

using namespace gcm;

// instantiations

template struct Node<VelocitySigmaVariables<1>>;
static_assert(std::is_pod<Node<VelocitySigmaVariables<1>>>::value, "This type is designed to be the plain old data");
template struct Node<VelocitySigmaVariables<2>>;
static_assert(std::is_pod<Node<VelocitySigmaVariables<2>>>::value, "This type is designed to be the plain old data");
template struct Node<VelocitySigmaVariables<3>>;
static_assert(std::is_pod<Node<VelocitySigmaVariables<3>>>::value, "This type is designed to be the plain old data");

template struct Node<VectorEnvelope<linal::Vector<2, VelocitySigmaVariables<1>>>,
		GcmMatricesPtrEnvelope<GcmMatrices<2, 1, IsotropicMaterial>>>;
template struct Node<VectorEnvelope<linal::Vector<5, VelocitySigmaVariables<2>>>,
		GcmMatricesPtrEnvelope<GcmMatrices<5, 2, IsotropicMaterial>>>;
template struct Node<VectorEnvelope<linal::Vector<9, VelocitySigmaVariables<3>>>,
        GcmMatricesPtrEnvelope<GcmMatrices<9, 3, IsotropicMaterial>>>;



// static initialization to prevent undefined reference
template<typename... Components> MPI::Datatype Node<Components ...>::MPI_NODE_TYPE;