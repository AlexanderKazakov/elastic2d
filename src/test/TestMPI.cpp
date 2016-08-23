#include <gtest/gtest.h>

#include <lib/Engine.hpp>
#include <lib/mesh/DefaultMesh.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>

#include <lib/util/task/Task.hpp>
#include <lib/rheology/models/models.hpp>


using namespace gcm;

template<class Value>
class TestMpiConnection : public testing::Test {
protected:
	void testMpiConnection() {
		int rank = MPI::COMM_WORLD.Get_rank();
		int numberOfWorkers = MPI::COMM_WORLD.Get_size();

		const int testNumberOfNodes = 1000;

		Value lefValues[testNumberOfNodes];
		Value righValues[testNumberOfNodes];

		for (int k = 0; k < testNumberOfNodes; k++) {
			for (int i = 0; i < Value::M; i++) {
				lefValues[k](i) = righValues[k](i) = rank;
			}
		}

		if (numberOfWorkers > 1) {
			if (rank == 0) {
				MPI_Sendrecv_replace(righValues, (int) sizeof(Value) *
				                     testNumberOfNodes, MPI_BYTE,
				                     rank + 1, 1, rank + 1, 1, MPI::COMM_WORLD,
				                     MPI_STATUS_IGNORE);
			} else if (rank == numberOfWorkers - 1) {
				MPI_Sendrecv_replace(lefValues, (int) sizeof(Value) *
				                     testNumberOfNodes, MPI_BYTE,
				                     rank - 1, 1, rank - 1, 1, MPI::COMM_WORLD,
				                     MPI_STATUS_IGNORE);
			} else {
				MPI_Sendrecv_replace(righValues, (int) sizeof(Value) *
				                     testNumberOfNodes, MPI_BYTE,
				                     rank + 1, 1, rank + 1, 1, MPI::COMM_WORLD,
				                     MPI_STATUS_IGNORE);
				MPI_Sendrecv_replace(lefValues, (int) sizeof(Value) *
				                     testNumberOfNodes, MPI_BYTE,
				                     rank - 1, 1, rank - 1, 1, MPI::COMM_WORLD,
				                     MPI_STATUS_IGNORE);
			}

			for (int k = 0; k < testNumberOfNodes; k++) {
				for (int i = 0; i < Value::M; i++) {
					if (rank == 0) {
						ASSERT_EQ(righValues[k](i), rank + 1);
					} else if (rank == numberOfWorkers - 1) {
						ASSERT_EQ(lefValues[k](i), rank - 1);
					} else {
						ASSERT_EQ(lefValues[k](i), rank - 1);
						ASSERT_EQ(righValues[k](i), rank + 1);
					}
				}
			}
		}
	}

};


/** Look at https://github.com/google/googletest/blob/master/googletest/samples/sample6_unittest.cc
  for explaination */
#if GTEST_HAS_TYPED_TEST_P
using testing::Types;
TYPED_TEST_CASE_P(TestMpiConnection);
TYPED_TEST_P(TestMpiConnection, MPI_NODE_TYPE) {
	this->testMpiConnection();
}
REGISTER_TYPED_TEST_CASE_P(TestMpiConnection, MPI_NODE_TYPE);

// write in generics all the Node implementations using in mpi connections
typedef Types<
        ElasticModel<1>::PdeVector,
        ElasticModel<2>::PdeVector,
        ElasticModel<3>::PdeVector
        > AllImplementations;

INSTANTIATE_TYPED_TEST_CASE_P(AllNodeTypes, TestMpiConnection, AllImplementations);
#endif // GTEST_HAS_TYPED_TEST_P


TEST(MPI, MpiEngineVsSequenceEngine) {
	Task task;

	task.globalSettings.dimensionality = 2;
	task.globalSettings.gridId = Grids::T::CUBIC;

	task.bodies = {
		{0, {Materials::T::ISOTROPIC, Models::T::ELASTIC}}
	};

	task.cubicGrid.borderSize = 2;
	task.cubicGrid.sizes = {20, 10};
	task.cubicGrid.lengths = {2, 1};

	task.globalSettings.CourantNumber = 1.8;
	task.materialConditions.byAreas.defaultMaterial = 
			std::make_shared<IsotropicMaterial>(4, 2, 0.5);
	task.globalSettings.numberOfSnaps = 5;

	Task::InitialCondition::Quantity pressure;
	pressure.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	pressure.value = 2.0;
	pressure.area = std::make_shared<SphereArea>(0.2, Real3({1, 0.5, 0}));
	task.initialCondition.quantities.push_back(pressure);

	// calculate in sequence
	task.globalSettings.forceSequence = true;
	Engine sequenceEngine(task);
	sequenceEngine.run();

	// calculate in parallel
	task.globalSettings.forceSequence = false;
	Engine mpiEngine(task);
	mpiEngine.run();
	
	
	struct Wrapper {
		typedef DefaultMesh<ElasticModel<2>, CubicGrid<2>, IsotropicMaterial> Mesh;
		static const Mesh* getMesh(const Engine& engine) {
			const AbstractGrid* grid = engine.getSolver()->getAbstractMesh();
			const Mesh* mesh = dynamic_cast<const Mesh*>(grid);
			assert_true(mesh);
			return mesh;
		}
	};
	
	// check that parallel result is equal to sequence result
	auto mpiMesh = Wrapper::getMesh(mpiEngine);
	auto sequenceMesh = Wrapper::getMesh(sequenceEngine);

	int numberOfNodesAlongXPerOneCore = CubicGrid<2>::
			numberOfNodesAlongXPerOneCore(task.cubicGrid);
	int startX = Mpi::Rank() * numberOfNodesAlongXPerOneCore;
	
	for (int x = 0; x < mpiMesh->sizes(0); x++) {
		for (int y = 0; y < mpiMesh->sizes(1); y++) {
			ASSERT_EQ(mpiMesh->pde({x, y}), sequenceMesh->pde({x + startX, y})) 
					<< "x = " << x << " global x = " << x + startX << " y = " << y
					<< "\nMPI:\n" << mpiMesh->pde({x, y}) 
					<< "\nsequence:\n" << sequenceMesh->pde({x + startX, y});
		}
	}
}


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	testing::InitGoogleTest(&argc, argv);
	int allTestsResult = RUN_ALL_TESTS();

	MPI_Finalize();
	return allTestsResult;
}


