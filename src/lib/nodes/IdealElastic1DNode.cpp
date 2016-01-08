#include "lib/nodes/IdealElastic1DNode.hpp"

using namespace gcm;

// initialization of static members

const int IdealElastic1DContainer::SIZE;
const int IdealElastic1DContainer::V_SIZE;
const int IdealElastic1DContainer::S_SIZE;

MPI::Datatype IdealElastic1DNode::MPI_NODE_TYPE;

const std::map<const std::string /* name */, const std::pair<const int /* index */, const int /* size */>> IdealElastic1DNode::VECTORS = {{"Velocity", {0, V_SIZE}}};

const std::map<const std::string /* name */, const int /* index */> IdealElastic1DNode::SCALARS = {{"Sxx", V_SIZE}};