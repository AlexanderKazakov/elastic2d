#include "IdealElastic3DModel.hpp"

using namespace gcm;


const std::map<const std::string /* name */, const std::pair<const int /* index */, const int /* size */>>
		IdealElastic3DModel::VECTORS = {{"Velocity", {0, IdealElastic3DModel::Node::Vector::V_SIZE}}};

const std::map<const std::string /* name */, const int /* index */>
		IdealElastic3DModel::SCALARS = {{"Sxx", IdealElastic3DModel::Node::Vector::V_SIZE},
		                                {"Sxy", IdealElastic3DModel::Node::Vector::V_SIZE + 1},
		                                {"Sxz", IdealElastic3DModel::Node::Vector::V_SIZE + 2},
		                                {"Syy", IdealElastic3DModel::Node::Vector::V_SIZE + 3},
		                                {"Syz", IdealElastic3DModel::Node::Vector::V_SIZE + 4},
		                                {"Szz", IdealElastic3DModel::Node::Vector::V_SIZE + 5}};