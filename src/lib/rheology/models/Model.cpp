#include <lib/rheology/models/Model.hpp>

using namespace gcm;

// instantiations

template class Model<VelocitySigmaVariables<1>, IsotropicMaterial>;
template class Model<VelocitySigmaVariables<2>, IsotropicMaterial>;
template class Model<VelocitySigmaVariables<3>, IsotropicMaterial>;
template class Model<VelocitySigmaVariables<3>, OrthotropicMaterial>;

