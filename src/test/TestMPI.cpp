#include <gtest/gtest.h>
#include <test/wrappers/Wrappers.hpp>

#include <lib/util/DataBus.hpp>
#include <lib/util/task/Task.hpp>
#include <lib/rheology/models/Model.hpp>


using namespace gcm;

template<class TNode>
class TestMpiNodeTypes : public testing::Test {
protected:
	void testMpiNodeTypes() {
		int rank = MPI::COMM_WORLD.Get_rank();
		int numberOfWorkers = MPI::COMM_WORLD.Get_size();

		const int testNumberOfNodes = 1000;

		TNode leftNodes[testNumberOfNodes];
		TNode rightNodes[testNumberOfNodes];

		for (int k = 0; k < testNumberOfNodes; k++) {
			for (int i = 0; i < TNode::Vector::M; i++) {
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
				for (int i = 0; i < TNode::Vector::M; i++) {
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
TYPED_TEST_CASE_P(TestMpiNodeTypes);
TYPED_TEST_P(TestMpiNodeTypes, MPI_NODE_TYPE) {
	this->testMpiNodeTypes();
}
REGISTER_TYPED_TEST_CASE_P(TestMpiNodeTypes, MPI_NODE_TYPE);

// write in generics all the Node implementations using in mpi connections
typedef Types<
		Node<Elastic1DModel>,
		Node<Elastic2DModel>,
		Node<Elastic3DModel>,
		Node<OrthotropicElastic3DModel>,
		Node<ContinualDamageElastic2DModel>,
		Node<SuperDuperModel>> AllImplementations;

INSTANTIATE_TYPED_TEST_CASE_P(AllNodeTypes, TestMpiNodeTypes, AllImplementations);
#endif // GTEST_HAS_TYPED_TEST_P


TEST(MPI, MpiEngineVsSequenceEngine)
{
	Task task;
	task.accuracyOrder = 2;
	task.CourantNumber = 1.8;
	task.isotropicMaterial = IsotropicMaterial(4.0, 2.0, 0.5);

	task.sizes(0) = 20;
	task.sizes(1) = 10;
	task.lengthes = {2, 1, 1};
	task.numberOfSnaps = 5;

	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 2.0;
	pressure.area = std::make_shared<SphereArea>(0.2, linal::Vector3({1, 0.5, 0}));
	task.initialCondition.quantities.push_back(pressure);

	task.borderConditions.at(CUBIC_BORDERS::X_LEFT) = BorderCondition::T::FREE_BORDER;

	// calculate in sequence
	task.forceSequence = true;
	EngineWrapper<StructuredGrid<Elastic2DModel>> sequenceEngine;
	sequenceEngine.setSolver(new DefaultSolver<StructuredGrid<Elastic2DModel>>());
	sequenceEngine.initialize(task);
	sequenceEngine.run();

	// calculate in parallel
	task.forceSequence = false;
	task.enableSnapshotting = true;
	EngineWrapper<StructuredGrid<Elastic2DModel>> mpiEngine;
	mpiEngine.setSolver(new DefaultSolver<StructuredGrid<Elastic2DModel>>());
	mpiEngine.initialize(task);
	mpiEngine.run();

	// check that parallel result is equal to sequence result
	auto mpiMesh = mpiEngine.getSolverForTest()->getMesh();
	auto sequenceMesh = sequenceEngine.getSolverForTest()->getMesh();
	for (int x = 0; x < mpiMesh->getXForTest(); x++) {
		for (int y = 0; y < mpiMesh->getYForTest(); y++) {
			ASSERT_EQ(mpiMesh->getNodeForTest(x, y, 0).u,
			          sequenceMesh->getNodeForTest (x + mpiMesh->getStartXForTest(), y, 0).u)
					<< "x = " << x << " global x = " << x + mpiMesh->getStartXForTest() << " y = " << y << std::endl;
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
