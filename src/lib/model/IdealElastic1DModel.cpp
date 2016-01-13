#include "IdealElastic1DModel.hpp"

using namespace gcm;


const std::map<const std::string /* name */, const std::pair<const int /* index */, const int /* size */>>
		IdealElastic1DModel::VECTORS = {{"Velocity", {0, IdealElastic1DModel::Node::V_SIZE}}};

const std::map<const std::string /* name */, const int /* index */>
		IdealElastic1DModel::SCALARS = {{"Sxx", IdealElastic1DModel::Node::V_SIZE}};