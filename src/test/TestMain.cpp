#include <gtest/gtest.h>

#include "test/TestMatrix.cpp"
#include "test/TestNode.cpp"
#include "test/TestInterpolator.cpp"
#include "test/TestMesh.cpp"
#include "test/TestSolver.cpp"


int main(int argc, char **argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
