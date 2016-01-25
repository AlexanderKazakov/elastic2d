#include <gtest/gtest.h>

#include <lib/util/DataBus.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/numeric/gcmethod/MpiStructuredSolver.hpp>

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

	task.sizes(0) = 20;
	task.sizes(1) = 10;
	task.lengthes = {2, 1, 1};
	task.numberOfSnaps = 5;

	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 2.0;

	pressure.area = std::make_shared<SphereArea>(0.2, linal::Vector3({1, 0.5, 0}));

	task.initialCondition.quantities.push_back(pressure);


//	task.borderConditions.at(CUBIC_BORDERS::X_LEFT) = BorderCondition::T::FREE_BORDER;

	// calculate in sequence
	task.forceSequence = true;
	MpiStructuredSolver<IdealElastic2DNode> sequenceSolver;
	sequenceSolver.initialize(task);
	sequenceSolver.calculate();

	// calculate in parallel
	task.forceSequence = false;
	task.enableSnapshotting = true;
	MpiStructuredSolver<IdealElastic2DNode> mpiSolver;
	mpiSolver.initialize(task);
	mpiSolver.calculate();

	// check that parallel result is equal to sequence result
	for (int x = 0; x < mpiSolver.getMesh()->getXForTest(); x++) {
		for (int y = 0; y < mpiSolver.getMesh()->getYForTest(); y++) {
			ASSERT_EQ(mpiSolver.getMesh()->getNodeForTest(x, y, 0).u,
			          sequenceSolver.getMesh()->getNodeForTest
					          (x + mpiSolver.getMesh()->getStartXForTest(), y, 0).u) <<
								"x = " << x <<
								" global x = " << x + mpiSolver.getMesh()->getStartXForTest()
			                      << " y = " << y << std::endl;
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
