#include <gtest/gtest.h>

#include "lib/grid/StructuredGrid.hpp"
#include "lib/model/IdealElastic2DModel.hpp"

using namespace gcm;

TEST(StructuredGrid, initialize) {
	Task task;
	task.accuracyOrder = 3;
	task.CourantNumber = 0.7;
	task.lambda0 = 2.0;
	task.mu0 = 0.5;
	task.rho0 = 4.0;
	task.X = 7;
	task.Y = 9;
	task.xLength = 20.0;
	task.yLength = 8.0;
	task.numberOfSnaps = 5;
	task.T = 100.0;
	task.initialConditions = InitialConditions::Zero;

	StructuredGrid<IdealElastic2DModel> structuredGrid;
	structuredGrid.initialize(task);
	ASSERT_NEAR(structuredGrid.getH0ForTest(), 3.333333333, EQUALITY_TOLERANCE);
	ASSERT_NEAR(structuredGrid.getH1ForTest(), 1.0, EQUALITY_TOLERANCE);
	ASSERT_NEAR(structuredGrid.getMaximalLambda(), 0.866025404, EQUALITY_TOLERANCE);
	ASSERT_NEAR(structuredGrid.getMinimalSpatialStep(), 1.0, EQUALITY_TOLERANCE);
	for (int x = 0; x < task.X; x++) {
		for (int y = 0; y < task.Y; y++) {
			for (int z = 0; z < task.Z; z++) {
				for (int i = 0; i < IdealElastic2DModel::Node::M; i++) {
					ASSERT_EQ(structuredGrid.getNodeForTest(x, y, z)(i), 0.0);
				}
			}
		}
	}
}


TEST(StructuredGrid, findSourcesForInterpolation)
{
	Task task;
	task.lambda0 = 2.0;
	task.mu0 = 1.0;
	task.rho0 = 2.0;

	task.X = task.Y = 3;
	task.accuracyOrder = 1;
	task.xLength = task.yLength = 2.0; // h_x = h_y = 1.0
	task.initialConditions = InitialConditions::TestExplosion;

	for (int stage = 0; stage <= 1; stage++) {

		StructuredGrid<IdealElastic2DModel> structuredGrid;
		structuredGrid.initialize(task);
		for (int x = 0; x < task.X; x++) {
			for (int y = 0; y < task.Y; y++) {
				// check that TestExplosion is set properly
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).Vx, 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).Vy, 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).Sxx, (x == 1 && y == 1) ? 1.0 : 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).Sxy, 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).Syy, (x == 1 && y == 1) ? 1.0 : 0.0);
			}
		}

		std::vector<IdealElastic2DModel::Node::Vector> src(task.accuracyOrder + 1);
		for (real dx = -1.0; dx <= 1.0; dx += 0.5) {
			structuredGrid.findSourcesForInterpolation(stage, 1, 1, 0, dx, src);
			for (int i = 0; i < IdealElastic2DModel::Node::M; i++) {
				ASSERT_EQ(src[0](i), (i == 2 || i == 4) ? 1.0 : 0.0);
				ASSERT_EQ(src[1](i), 0.0);
			}
		}

	}
}


TEST(StructuredGrid, interpolateValuesAround)
{
	Task task;
	task.lambda0 = 2.0;
	task.mu0 = 1.0;
	task.rho0 = 2.0;

	task.X = task.Y = 3;
	task.accuracyOrder = 1;
	task.xLength = task.yLength = 2.0; // h_x = h_y = 1.0
	task.initialConditions = InitialConditions::TestExplosion;

	for (int stage = 0; stage <= 1; stage++) {

		StructuredGrid<IdealElastic2DModel> structuredGrid;
		structuredGrid.initialize(task);
		for (int x = 0; x < task.X; x++) {
			for (int y = 0; y < task.Y; y++) {
				// check that TestExplosion is set properly
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).Vx, 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).Vy, 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).Sxx, (x == 1 && y == 1) ? 1.0 : 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).Sxy, 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).Syy, (x == 1 && y == 1) ? 1.0 : 0.0);
			}
		}

		IdealElastic2DModel::Node::Vector dx({-1, 1, -0.5, 0.5, 0});
		IdealElastic2DModel::Node::Matrix matrix = structuredGrid.interpolateValuesAround(stage, 1, 1, 0, dx);

		for (int i = 0; i < IdealElastic2DModel::Node::M; i++) {
			ASSERT_EQ(matrix(i, 0), 0.0) << "i = " << i; // Courant = 1
			ASSERT_EQ(matrix(i, 1), 0.0) << "i = " << i; // Courant = 1
			ASSERT_EQ(matrix(0, i), 0.0) << "i = " << i; // Vx
			ASSERT_EQ(matrix(1, i), 0.0) << "i = " << i; // Vy
			ASSERT_EQ(matrix(3, i), 0.0) << "i = " << i; // Sxy
		}

		ASSERT_EQ(matrix(2, 2), 0.5); // Courant = 0.5
		ASSERT_EQ(matrix(2, 3), 0.5); // Courant = 0.5
		ASSERT_EQ(matrix(4, 2), 0.5); // Courant = 0.5
		ASSERT_EQ(matrix(4, 3), 0.5); // Courant = 0.5

		ASSERT_EQ(matrix(2, 4), 1.0); // Courant = 0
		ASSERT_EQ(matrix(4, 4), 1.0); // Courant = 0
	}
}
