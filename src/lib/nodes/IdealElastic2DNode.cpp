#include "lib/nodes/IdealElastic2DNode.hpp"

using namespace gcm;

// initialization of static members

const int IdealElastic2DContainer::SIZE;
const int IdealElastic2DContainer::V_SIZE;
const int IdealElastic2DContainer::S_SIZE;

MPI::Datatype IdealElastic2DNode::MPI_NODE_TYPE;
