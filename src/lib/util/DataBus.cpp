#include <lib/nodes/IdealElastic2DNode.hpp>
#include "lib/DataBus.hpp"

using namespace gcm;


void DataBus::createStaticTypes() {
	createStaticType<IdealElastic2DNode>();
}

template<class TNode>
void createStaticType() {
	TNode* nodes = new TNode[2];

	// Node
	MPI::Datatype node_types[] = {
			MPI::LB,
#if LIBGCM_DOUBLE_PRECISION
			MPI::DOUBLE,
#else
			MPI::FLOAT
#endif
			MPI::UB
	};

	int node_lengths[] = {1, TNode::M, 1};

	MPI::Aint node_disp[] = {
			MPI::Get_address(&nodes[0]),
			MPI::Get_address(&(nodes[0](0))),
			MPI::Get_address(&nodes[1])
	};
	for(int i = 3; i >= 0; i--)
		node_disp[i] -= MPI::Get_address(&nodes[0]);

	TNode::MPI_NODE_TYPE = MPI::Datatype::Create_struct(3, node_lengths, node_disp, node_types);

	TNode::MPI_NODE_TYPE.Commit();
}


