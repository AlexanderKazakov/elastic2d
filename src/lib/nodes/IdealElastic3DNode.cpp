#include "lib/nodes/IdealElastic3DNode.hpp"

using namespace gcm;

// initialization of static members

const int IdealElastic3DContainer::SIZE;
const int IdealElastic3DContainer::V_SIZE;
const int IdealElastic3DContainer::S_SIZE;

MPI::Datatype IdealElastic3DNode::MPI_NODE_TYPE;
