#include "lib/nodes/IdealElastic2DNode.hpp"

using namespace gcm;

// initialization of static members

const int IdealElastic2DContainer::SIZE;
const int IdealElastic2DContainer::V_SIZE;
const int IdealElastic2DContainer::S_SIZE;

MPI::Datatype IdealElastic2DNode::MPI_NODE_TYPE;

const std::map<const std::string /* name */, const std::pair<const int /* index */, const int /* size */>> IdealElastic2DNode::VECTORS = {{"Velocity", {0, V_SIZE}}};

const std::map<const std::string /* name */, const int /* index */> IdealElastic2DNode::SCALARS = {{"Sxx", V_SIZE}, {"Sxy", V_SIZE + 1}, {"Syy", V_SIZE + 2}};