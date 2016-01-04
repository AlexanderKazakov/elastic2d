#include <gtest/gtest.h>

#include "lib/config.hpp"
#include "lib/DataBus.hpp"
#include "lib/Task.hpp"
#include "lib/mesh/Mesh.hpp"
#include "lib/MPISolver.hpp"

using namespace gcm;

TEST(MPI, CustomDatatype)
{
	int rank = MPI::COMM_WORLD.Get_rank();
	int numberOfWorkers = MPI::COMM_WORLD.Get_size();

	const int testNumberOfNodes = 1000;
	
	IdealElastic2DNode leftNodes[testNumberOfNodes];
	IdealElastic2DNode rightNodes[testNumberOfNodes];

	for (int k = 0; k < testNumberOfNodes; k++) {
		for (int i = 0; i < IdealElastic2DNode::M; i++) {
			leftNodes[k](i) = rightNodes[k](i) = rank;
		}
	}

	if (numberOfWorkers > 1) {
		if (rank == 0) {
			MPI_Sendrecv_replace(rightNodes, testNumberOfNodes, IdealElastic2DNode::MPI_NODE_TYPE,
			                     rank + 1, 1, rank + 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		} else if (rank == numberOfWorkers - 1) {
			MPI_Sendrecv_replace(leftNodes, testNumberOfNodes, IdealElastic2DNode::MPI_NODE_TYPE,
			                     rank - 1, 1, rank - 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		} else {
			MPI_Sendrecv_replace(rightNodes, testNumberOfNodes, IdealElastic2DNode::MPI_NODE_TYPE,
			                     rank + 1, 1, rank + 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv_replace(leftNodes, testNumberOfNodes, IdealElastic2DNode::MPI_NODE_TYPE,
			                     rank - 1, 1, rank - 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		}

		for (int k = 0; k < testNumberOfNodes; k++) {
			for (int i = 0; i < IdealElastic2DNode::M; i++) {
				if (rank == 0) {
					ASSERT_EQ(rightNodes[k](i), rank + 1);
				} else if (rank == numberOfWorkers - 1) {
					ASSERT_EQ(leftNodes[k](i), rank - 1);
				} else {
					ASSERT_EQ(leftNodes[k](i), rank - 1);
					ASSERT_EQ(rightNodes[k](i), rank + 1);
				}
			}
		}
	}
}


TEST(MPI, MPISolverVsSequenceSolver)
{
	Task task;
	task.accuracyOrder = 2;
	task.CourantNumber = 1.8;
	task.lambda0 = 2.0;
	task.mu0 = 0.5;
	task.rho0 = 4.0;
	task.X = 50;
	task.Y = 70;
	task.xLength = 4.0;
	task.yLength = 6.0;
	task.numberOfSnaps = 50;
	task.initialConditions = InitialConditions::Explosion;

	// calculate in sequence
	Mesh *meshSeq = new Mesh();
	Mesh *newMeshSeq = new Mesh();
	meshSeq->initialize(task, true);
	newMeshSeq->initialize(task, true);
	MPISolver sequenceSolver(meshSeq, newMeshSeq);
	sequenceSolver.calculate();

	// calculate in parallel
	Mesh *mesh = new Mesh();
	Mesh *newMesh = new Mesh();
	mesh->initialize(task);
	newMesh->initialize(task);
	MPISolver mpiSolver(mesh, newMesh);
	mpiSolver.calculate();

	// check that parallel result is equal to sequence result
	for (int y = 0; y < mesh->getYForTest(); y++) {
		for (int x = 0; x < mesh->getXForTest(); x++) {
			ASSERT_EQ(mesh->getNodeForTest(y, x).u,
			          meshSeq->getNodeForTest(y + mesh->getStartYForTest(), x).u);
		}
	}
}


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();

	testing::InitGoogleTest(&argc, argv);
	int allTestsResult = RUN_ALL_TESTS();

	MPI_Finalize();
	return allTestsResult;
}
