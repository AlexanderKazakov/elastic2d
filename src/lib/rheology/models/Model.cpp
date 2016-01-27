#include <lib/rheology/models/Model.hpp>

using namespace gcm;

// instantiations

template struct Model<VelocitySigmaVariables<1>>;
static_assert(std::is_pod<Model<VelocitySigmaVariables<1>>>::value, "This type is designed to be the plain old data");
template struct Model<VelocitySigmaVariables<2>>;
static_assert(std::is_pod<Model<VelocitySigmaVariables<2>>>::value, "This type is designed to be the plain old data");
template struct Model<VelocitySigmaVariables<3>>;
static_assert(std::is_pod<Model<VelocitySigmaVariables<3>>>::value, "This type is designed to be the plain old data");

template struct Model<VectorEnvelope<linal::Vector<2, VelocitySigmaVariables<1>>>,
		GcmMatricesPtrEnvelope<GcmMatrices<VelocitySigmaVariables<1>, IsotropicMaterial>>>;
template struct Model<VectorEnvelope<linal::Vector<5, VelocitySigmaVariables<2>>>,
		GcmMatricesPtrEnvelope<GcmMatrices<VelocitySigmaVariables<2>, IsotropicMaterial>>>;
template struct Model<VectorEnvelope<linal::Vector<9, VelocitySigmaVariables<3>>>,
        GcmMatricesPtrEnvelope<GcmMatrices<VelocitySigmaVariables<3>, IsotropicMaterial>>>;

