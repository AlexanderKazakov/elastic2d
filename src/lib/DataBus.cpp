#include "lib/DataBus.hpp"
#include "lib/Node.hpp"

using namespace gcm;

MPI::Datatype DataBus::MPI_NODE; // zero initialization of static member


void DataBus::createStaticTypes() {
	Node* nodes = new Node[2];

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

	int node_lengths[] = {1, 5, 1};

	MPI::Aint node_disp[] = {
			MPI::Get_address(&nodes[0]),
			MPI::Get_address(&(nodes[0].u(0))),
			MPI::Get_address(&nodes[1])
	};
	for(int i = 3; i >= 0; i--)
		node_disp[i] -= MPI::Get_address(&nodes[0]);

	MPI_NODE = MPI::Datatype::Create_struct(3, node_lengths, node_disp, node_types);

	MPI_NODE.Commit();
}



