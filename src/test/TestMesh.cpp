#include <gtest/gtest.h>

#include "lib/Mesh.hpp"


TEST(Mesh, initialize) {
	Task task;
	task.accuracyOrder = 3;
	task.CourantNumber = 0.7;
	task.lambda0 = 2.0;
	task.mu0 = 0.5;
	task.rho0 = 4.0;
	task.X = 7;
	task.Y = 4;
	task.xLength = 20.0;
	task.yLength = 3.0;
	task.numberOfSnaps = 5;
	task.T = 100.0;
	task.initialConditions = InitialConditions::Zero;

	Mesh mesh;
	mesh.initialize(task);
	ASSERT_NEAR(mesh.getH0ForTest(), 3.333333333, EQUALITY_TOLERANCE);
	ASSERT_NEAR(mesh.getH1ForTest(), 1.0, EQUALITY_TOLERANCE);
	ASSERT_NEAR(mesh.getTauForTest(), 0.808290377, EQUALITY_TOLERANCE);
	ASSERT_NEAR(mesh.getTForTest(), 4.0414518843, EQUALITY_TOLERANCE);
	for (int x = 0; x < task.X; x++) {
		for (int y = 0; y < task.Y; y++) {
			for (int i = 0; i < N; i++) {
				ASSERT_EQ(mesh.getNodeForTest(y, x).u.get(i), 0.0);
			}
		}
	}
}


TEST(Mesh, findSourcesForInterpolation)
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

		Mesh mesh;
		mesh.initialize(task);
		for (int x = 0; x < task.X; x++) {
			for (int y = 0; y < task.Y; y++) {
				// check that TestExplosion is set properly
				ASSERT_EQ(mesh.getNodeForTest(y, x).get(NodeMap::Vx), 0.0);
				ASSERT_EQ(mesh.getNodeForTest(y, x).get(NodeMap::Vy), 0.0);
				ASSERT_EQ(mesh.getNodeForTest(y, x).get(NodeMap::Sxx), (x == 1 && y == 1) ? 1.0 : 0.0);
				ASSERT_EQ(mesh.getNodeForTest(y, x).get(NodeMap::Sxy), 0.0);
				ASSERT_EQ(mesh.getNodeForTest(y, x).get(NodeMap::Syy), (x == 1 && y == 1) ? 1.0 : 0.0);
			}
		}

		std::vector<Vector> src(task.accuracyOrder + 1);
		for (real dx = -1.0; dx <= 1.0; dx += 0.5) {
			mesh.findSourcesForInterpolation(stage, 1, 1, dx, src);
			for (int i = 0; i < N; i++) {
				ASSERT_EQ(src[0].get(i), (i == 2 || i == 4) ? 1.0 : 0.0);
				ASSERT_EQ(src[1].get(i), 0.0);
			}
		}

	}
}


TEST(Mesh, interpolateValuesAround)
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

		Mesh mesh;
		mesh.initialize(task);
		for (int x = 0; x < task.X; x++) {
			for (int y = 0; y < task.Y; y++) {
				// check that TestExplosion is set properly
				ASSERT_EQ(mesh.getNodeForTest(y, x).get(NodeMap::Vx), 0.0);
				ASSERT_EQ(mesh.getNodeForTest(y, x).get(NodeMap::Vy), 0.0);
				ASSERT_EQ(mesh.getNodeForTest(y, x).get(NodeMap::Sxx), (x == 1 && y == 1) ? 1.0 : 0.0);
				ASSERT_EQ(mesh.getNodeForTest(y, x).get(NodeMap::Sxy), 0.0);
				ASSERT_EQ(mesh.getNodeForTest(y, x).get(NodeMap::Syy), (x == 1 && y == 1) ? 1.0 : 0.0);
			}
		}

		Vector dx;
		dx.createVector({-1, 1, -0.5, 0.5, 0});
		Matrix matrix = mesh.interpolateValuesAround(stage, 1, 1, dx);

		for (int i = 0; i < N; i++) {
			ASSERT_EQ(matrix.get(i, 0), 0.0) << "i = " << i; // Courant = 1
			ASSERT_EQ(matrix.get(i, 1), 0.0) << "i = " << i; // Courant = 1
			ASSERT_EQ(matrix.get(0, i), 0.0) << "i = " << i; // Vx
			ASSERT_EQ(matrix.get(1, i), 0.0) << "i = " << i; // Vy
			ASSERT_EQ(matrix.get(3, i), 0.0) << "i = " << i; // Sxy
		}

		ASSERT_EQ(matrix.get(2, 2), 0.5); // Courant = 0.5
		ASSERT_EQ(matrix.get(2, 3), 0.5); // Courant = 0.5
		ASSERT_EQ(matrix.get(4, 2), 0.5); // Courant = 0.5
		ASSERT_EQ(matrix.get(4, 3), 0.5); // Courant = 0.5

		ASSERT_EQ(matrix.get(2, 4), 1.0); // Courant = 0
		ASSERT_EQ(matrix.get(4, 4), 1.0); // Courant = 0

	}
}
