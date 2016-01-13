#include "lib/nodes/IdealElastic1DNode.hpp"

using namespace gcm;

// initialization of static members

const int IdealElastic1DContainer::SIZE;
const int IdealElastic1DContainer::V_SIZE;
const int IdealElastic1DContainer::S_SIZE;

MPI::Datatype IdealElastic1DNode::MPI_NODE_TYPE;
