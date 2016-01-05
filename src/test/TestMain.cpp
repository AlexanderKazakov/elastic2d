#include <gtest/gtest.h>
#include <mpi.h>

#include "lib/util/DataBus.hpp"

#include "test/gcmMatrices.cpp"
#include "test/TestNode.cpp"
#include "test/TestInterpolator.cpp"
#include "test/TestMesh.cpp"
#include "test/TestSolver.cpp"
#include "test/assertions.cpp"
#include "test/linal.cpp"

using namespace gcm;

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();

	testing::InitGoogleTest(&argc, argv);
	int allTestsResult = RUN_ALL_TESTS();

	MPI_Finalize();
	return allTestsResult;
}
