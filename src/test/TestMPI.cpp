#include <gtest/gtest.h>

#include <lib/util/DataBus.hpp>
#include <lib/task/Task.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/solver/MpiStructuredSolver.hpp>
#include <lib/model/IdealElastic2DModel.hpp>
#include <lib/nodes/Node.hpp>

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
				leftNodes[k].u(i) = rightNodes[k].u(i) = rank;
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
	task.material = IsotropicMaterial(4.0, 2.0, 0.5);
	task.X = 50;
	task.Y = 70;
	task.xLength = 4.0;
	task.yLength = 6.0;
	task.numberOfSnaps = 50;
	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 2.0;
	pressure.area = std::make_shared<SphereArea>(0.3, linal::Vector3({2, 2, 0}));
	task.initialCondition.quantities.push_back(pressure);

//	task.borderConditions.at(CUBIC_BORDERS::Y_LEFT) = BorderCondition::CONDITION::FREE_BORDER;

	// calculate in sequence
	task.forceSequence = true;
	MpiStructuredSolver<IdealElastic2DModel> sequenceSolver;
	sequenceSolver.initialize(task);
	sequenceSolver.calculate();

	// calculate in parallel
	task.forceSequence = false;
	task.enableSnapshotting = true;
	MpiStructuredSolver<IdealElastic2DModel> mpiSolver;
	mpiSolver.initialize(task);
	mpiSolver.calculate();

	// check that parallel result is equal to sequence result
	for (int x = 0; x < mpiSolver.getMesh()->getXForTest(); x++) {
		for (int y = 0; y < mpiSolver.getMesh()->getYForTest(); y++) {
			ASSERT_EQ(mpiSolver.getMesh()->getNodeForTest(x, y, 0).u,
			          sequenceSolver.getMesh()->getNodeForTest(x + mpiSolver.getMesh()->getStartXForTest(), y, 0).u);
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
