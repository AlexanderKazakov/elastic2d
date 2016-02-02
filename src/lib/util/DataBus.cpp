#include <lib/util/DataBus.hpp>
#include <lib/grid/nodes/Nodes.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;


template<class TNode>
void createStaticType() {
	TNode* nodes = new TNode[2];

	MPI::Datatype node_types[] = {
			MPI::LB,
#if LIBGCM_DOUBLE_PRECISION
			MPI::DOUBLE,
#else
			MPI::FLOAT,
#endif
			MPI::UB
	};

	int node_lengths[] = {1, TNode::Vector::M, 1};

	MPI::Aint node_disp[] = {
			MPI::Get_address(&nodes[0]),
			MPI::Get_address(&(nodes[0].u(0))),
			MPI::Get_address(&nodes[1])
	};
	for(int i = 3; i >= 0; i--)
		node_disp[i] -= MPI::Get_address(&nodes[0]);

	TNode::MPI_NODE_TYPE = MPI::Datatype::Create_struct(3, node_lengths, node_disp, node_types);

	TNode::MPI_NODE_TYPE.Commit();

	delete [] nodes;
}

void DataBus::createStaticTypes() {
	createStaticType<Node<Elastic1DModel>>();
	createStaticType<Node<Elastic2DModel>>();
	createStaticType<Node<Elastic3DModel>>();
	createStaticType<Node<OrthotropicElastic3DModel>>();
	createStaticType<Node<ContinualDamageElastic2DModel>>();
	createStaticType<Node<SuperDuperModel>>();
}
