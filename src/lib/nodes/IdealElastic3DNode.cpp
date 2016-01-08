#include "lib/nodes/IdealElastic3DNode.hpp"

using namespace gcm;

// initialization of static members

const int IdealElastic3DContainer::SIZE;
const int IdealElastic3DContainer::V_SIZE;
const int IdealElastic3DContainer::S_SIZE;

MPI::Datatype IdealElastic3DNode::MPI_NODE_TYPE;

const std::map<const std::string /* name */, const std::pair<const int /* index */, const int /* size */>> IdealElastic3DNode::VECTORS = {{"Velocity", {0, V_SIZE}}};

const std::map<const std::string /* name */, const int /* index */> IdealElastic3DNode::SCALARS = {{"Sxx", V_SIZE}, {"Sxy", V_SIZE + 1}, {"Sxz", V_SIZE + 2},
                                                                                                   {"Syy", V_SIZE + 3}, {"Syz", V_SIZE + 4}, {"Szz", V_SIZE + 5}};