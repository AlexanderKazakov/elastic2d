#include <gtest/gtest.h>
#include <test/wrappers/Wrappers.hpp>

#include <lib/util/task/Task.hpp>
#include <lib/rheology/models/Model.hpp>


using namespace gcm;

template<class TNode>
class TestMpiConnection : public testing::Test {
protected:
	void testMpiConnection() {
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
				MPI_Sendrecv_replace(rightNodes, (int) sizeof(TNode) * testNumberOfNodes, MPI_BYTE,
				                     rank + 1, 1, rank + 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
			} else if (rank == numberOfWorkers - 1) {
				MPI_Sendrecv_replace(leftNodes, (int) sizeof(TNode) * testNumberOfNodes, MPI_BYTE,
				                     rank - 1, 1, rank - 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
			} else {
				MPI_Sendrecv_replace(rightNodes, (int) sizeof(TNode) * testNumberOfNodes, MPI_BYTE,
				                     rank + 1, 1, rank + 1, 1, MPI::COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv_replace(leftNodes, (int) sizeof(TNode) * testNumberOfNodes, MPI_BYTE,
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
TYPED_TEST_CASE_P(TestMpiConnection);
TYPED_TEST_P(TestMpiConnection, MPI_NODE_TYPE) {
	this->testMpiConnection();
}
REGISTER_TYPED_TEST_CASE_P(TestMpiConnection, MPI_NODE_TYPE);

// write in generics all the Node implementations using in mpi connections
typedef Types<
		Elastic1DModel::/*Pde*/Vector,
		ContinualDamageElastic2DModel::/*Pde*/Vector,
		OrthotropicElastic3DModel::/*Pde*/Vector
> AllImplementations;

INSTANTIATE_TYPED_TEST_CASE_P(AllNodeTypes, TestMpiConnection, AllImplementations);
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
	EngineWrapper<DefaultGrid<Elastic2DModel, StructuredGrid>> sequenceEngine;
	sequenceEngine.setSolver(new DefaultSolver<DefaultGrid<Elastic2DModel, StructuredGrid>>());
	sequenceEngine.initialize(task);
	sequenceEngine.run();

	// calculate in parallel
	task.forceSequence = false;
	task.enableSnapshotting = true;
	EngineWrapper<DefaultGrid<Elastic2DModel, StructuredGrid>> mpiEngine;
	mpiEngine.setSolver(new DefaultSolver<DefaultGrid<Elastic2DModel, StructuredGrid>>());
	mpiEngine.initialize(task);
	mpiEngine.run();

	// check that parallel result is equal to sequence result
	auto mpiMesh = mpiEngine.getSolverForTest()->getMesh();
	auto sequenceMesh = sequenceEngine.getSolverForTest()->getMesh();

	int numberOfNodesAlongXPerOneCore = (int) std::round((real) task.sizes(0) / MPI::COMM_WORLD.Get_size());
	int startX = MPI::COMM_WORLD.Get_rank() * numberOfNodesAlongXPerOneCore;
	for (int x = 0; x < mpiMesh->getSizes()(0); x++) {
		for (int y = 0; y < mpiMesh->getSizes()(1); y++) {
			ASSERT_EQ(mpiMesh->getPde(x, y, 0), sequenceMesh->getPde(x + startX, y, 0))
					<< "x = " << x << " global x = " << x + startX << " y = " << y << std::endl;
		}
	}
}


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	testing::InitGoogleTest(&argc, argv);
	int allTestsResult = RUN_ALL_TESTS();

	MPI_Finalize();
	return allTestsResult;
}
