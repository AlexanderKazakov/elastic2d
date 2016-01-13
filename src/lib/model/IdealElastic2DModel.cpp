#include "IdealElastic2DModel.hpp"

using namespace gcm;


const std::map<const std::string /* name */, const std::pair<const int /* index */, const int /* size */>>
		IdealElastic2DModel::VECTORS = {{"Velocity", {0, IdealElastic2DModel::Node::Vector::V_SIZE}}};

const std::map<const std::string /* name */, const int /* index */>
		IdealElastic2DModel::SCALARS = {{"Sxx", IdealElastic2DModel::Node::Vector::V_SIZE},
		                                {"Sxy", IdealElastic2DModel::Node::Vector::V_SIZE + 1},
		                                {"Syy", IdealElastic2DModel::Node::Vector::V_SIZE + 2}};