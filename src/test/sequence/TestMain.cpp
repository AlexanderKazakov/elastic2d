#include <gtest/gtest.h>
#include <mpi.h>

#include "lib/util/DataBus.hpp"

using namespace gcm;

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	DataBus::createStaticTypes();

	testing::InitGoogleTest(&argc, argv);
	int allTestsResult = RUN_ALL_TESTS();

	MPI_Finalize();
	return allTestsResult;
}
