#include <gtest/gtest.h>

#include "lib/config.hpp"
#include "lib/util/DataBus.hpp"
#include "lib/Task.hpp"
#include "lib/grid/StructuredGrid.hpp"
#include "lib/solver/MpiStructuredSolver.hpp"
#include "lib/nodes/IdealElastic3DNode.hpp"
#include "lib/nodes/IdealElastic2DNode.hpp"
#include "lib/nodes/IdealElastic1DNode.hpp"
#include "lib/model/IdealElastic2DModel.hpp"

using namespace gcm;

template<class TNode>
class Test_MPI_Sendrecv_replace : public testing::Test {
protected:
	void test_MPI_Sendrecv_replace() {
		int rank = MPI::COMM_WORLD.Get_rank();
		int numberOfWorkers = MPI::COMM_WORLD.Get_size();

		const int testNumberOfNodes = 1000;

		TNode leftNodes[testNumberOfNodes];
		TNode rightNodes[testNumberOfNodes];

		for (int k = 0; k < testNumberOfNodes; k++) {
			for (int i = 0; i < TNode::M; i++) {
				leftNodes[k](i) = rightNodes[k](i) = rank;
			}
		}

		if (numberOfWorkers > 1) {
			if (rank == 0) {
				MPI_Sendrecv_replace(rightNodes, testNumberOfNodes, TNode::MPI_NODE_TYPE,
				                     rank + 1, 1, rank + 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
			} else if (rank == numberOfWorkers - 1) {
				MPI_Sendrecv_replace(leftNodes, testNumberOfNodes, TNode::MPI_NODE_TYPE,
				                     rank - 1, 1, rank - 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
			} else {
				MPI_Sendrecv_replace(rightNodes, testNumberOfNodes, TNode::MPI_NODE_TYPE,
				                     rank + 1, 1, rank + 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv_replace(leftNodes, testNumberOfNodes, TNode::MPI_NODE_TYPE,
				                     rank - 1, 1, rank - 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
			}

			for (int k = 0; k < testNumberOfNodes; k++) {
				for (int i = 0; i < TNode::M; i++) {
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
	};
};

/** Look at https://github.com/google/googletest/blob/master/googletest/samples/sample6_unittest.cc for explaination */
#if GTEST_HAS_TYPED_TEST_P
using testing::Types;
TYPED_TEST_CASE_P(Test_MPI_Sendrecv_replace);
TYPED_TEST_P(Test_MPI_Sendrecv_replace, MPI_NODE_TYPE) {
	this->test_MPI_Sendrecv_replace();
}
REGISTER_TYPED_TEST_CASE_P(Test_MPI_Sendrecv_replace, MPI_NODE_TYPE);

// write in generics all the Node implementations using in mpi connections
typedef Types<IdealElastic3DNode, IdealElastic2DNode, IdealElastic1DNode> AllImplementations;

INSTANTIATE_TYPED_TEST_CASE_P(AllNodeTypes, Test_MPI_Sendrecv_replace, AllImplementations);
#endif // GTEST_HAS_TYPED_TEST_P


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
	task.forceSequence = true;
	StructuredGrid<IdealElastic2DModel> *meshSeq = new StructuredGrid<IdealElastic2DModel>();
	StructuredGrid<IdealElastic2DModel> *newMeshSeq = new StructuredGrid<IdealElastic2DModel>();
	MpiStructuredSolver<IdealElastic2DModel> sequenceSolver;
	sequenceSolver.initialize(task, meshSeq, newMeshSeq);
	sequenceSolver.calculate();

	// calculate in parallel
	task.forceSequence = false;
	StructuredGrid<IdealElastic2DModel> *mesh = new StructuredGrid<IdealElastic2DModel>();
	StructuredGrid<IdealElastic2DModel> *newMesh = new StructuredGrid<IdealElastic2DModel>();
	MpiStructuredSolver<IdealElastic2DModel> mpiSolver;
	mpiSolver.initialize(task, mesh, newMesh);
	mpiSolver.calculate();

	// check that parallel result is equal to sequence result
	for (int y = 0; y < mesh->getYForTest(); y++) {
		for (int x = 0; x < mesh->getXForTest(); x++) {
			ASSERT_EQ(mesh->getNodeForTest(y, x), meshSeq->getNodeForTest(y + mesh->getStartYForTest(), x));
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
