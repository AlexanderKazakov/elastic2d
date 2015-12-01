#include <gtest/gtest.h>

#include "lib/config.hpp"
#include "lib/DataBus.hpp"
#include "lib/Node.hpp"


TEST(MPI, CustomDatatype)
{
	MPI_Init(nullptr, nullptr);
	DataBus::createStaticTypes();

	int rank = MPI::COMM_WORLD.Get_rank();
	int numberOfWorkers = MPI::COMM_WORLD.Get_size();

	const int testNumberOfNodes = 1000;
	
	Node leftNodes[testNumberOfNodes];
	Node rightNodes[testNumberOfNodes];

	for (uint k = 0; k < testNumberOfNodes; k++) {
		for (uint i = 0; i < N; i++) {
			leftNodes[k].u(i) = rightNodes[k].u(i) = rank;
		}
	}

	if (numberOfWorkers > 1) {
		if (rank == 0) {
			MPI_Sendrecv_replace(rightNodes, testNumberOfNodes, DataBus::MPI_NODE,
			                     rank + 1, 1, rank + 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		} else if (rank == numberOfWorkers - 1) {
			MPI_Sendrecv_replace(leftNodes, testNumberOfNodes, DataBus::MPI_NODE,
			                     rank - 1, 1, rank - 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		} else {
			MPI_Sendrecv_replace(rightNodes, testNumberOfNodes, DataBus::MPI_NODE,
			                     rank + 1, 1, rank + 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv_replace(leftNodes, testNumberOfNodes, DataBus::MPI_NODE,
			                     rank - 1, 1, rank - 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		}

		for (uint k = 0; k < testNumberOfNodes; k++) {
			for (uint i = 0; i < N; i++) {
				if (rank == 0) {
					ASSERT_EQ(rightNodes[k].u(i), rank + 1);
				} else if (rank == numberOfWorkers - 1) {
					ASSERT_EQ(leftNodes[k].u(i), rank - 1);
				} else {
					ASSERT_EQ(leftNodes[k].u(i), rank - 1);
					ASSERT_EQ(rightNodes[k].u(i), rank + 1);
				}
			}
		}
	}

	MPI_Finalize();
}
