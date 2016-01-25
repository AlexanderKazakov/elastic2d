#include <gtest/gtest.h>

#include <lib/util/areas/SphereArea.hpp>
#include <lib/grid/StructuredGrid.hpp>

using namespace gcm;

TEST(StructuredGrid, initialize) {
	Task task;
	task.accuracyOrder = 3;
	task.CourantNumber = 0.7;
	task.material = IsotropicMaterial(4.0, 2.0, 0.5);
	task.sizes(0) = 7;
	task.sizes(1) = 9;
	task.lengthes = {20, 8, 1};
	task.numberOfSnaps = 5;
	task.T = 100.0;

	StructuredGrid<IdealElastic2DNode> structuredGrid;
	structuredGrid.initialize(task);
	ASSERT_NEAR(structuredGrid.getH0ForTest(), 3.333333333, EQUALITY_TOLERANCE);
	ASSERT_NEAR(structuredGrid.getH1ForTest(), 1.0, EQUALITY_TOLERANCE);
	ASSERT_NEAR(structuredGrid.getMaximalLambda(), 0.866025404, EQUALITY_TOLERANCE);
	ASSERT_NEAR(structuredGrid.getMinimalSpatialStep(), 1.0, EQUALITY_TOLERANCE);
	for (int x = 0; x < task.sizes(0); x++) {
		for (int y = 0; y < task.sizes(1); y++) {
			for (int z = 0; z < task.sizes(2); z++) {
				for (int i = 0; i < IdealElastic2DNode::M; i++) {
					ASSERT_EQ(structuredGrid.getNodeForTest(x, y, z).u(i), 0.0);
				}
			}
		}
	}
}


TEST(StructuredGrid, findSourcesForInterpolation)
{
	Task task;
	task.material = IsotropicMaterial(2.0, 2.0, 1.0);

	task.sizes(0) = task.sizes(1) = 3;
	task.accuracyOrder = 1;
	task.lengthes = {2.0, 2.0, 1.0}; // h_x = h_y = 1.0

	Task::InitialCondition::Quantity quantity;
	quantity.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	quantity.value = -1.0;
	quantity.area = std::make_shared<SphereArea>(0.1, linal::Vector3({1,1,0}));
	task.initialCondition.quantities.push_back(quantity);

	for (int stage = 0; stage <= 1; stage++) {

		StructuredGrid<IdealElastic2DNode> structuredGrid;
		structuredGrid.initialize(task);
		for (int x = 0; x < task.sizes(0); x++) {
			for (int y = 0; y < task.sizes(1); y++) {
				// check that values is set properly
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).u.V[0], 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).u.V[1], 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).u.sigma(0, 0), (x == 1 && y == 1) ? 1.0 : 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).u.sigma(0, 1), 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).u.sigma(1, 1), (x == 1 && y == 1) ? 1.0 : 0.0);
			}
		}

		std::vector<IdealElastic2DNode::Vector> src((unsigned long)task.accuracyOrder + 1);
		for (real dx = -1.0; dx <= 1.0; dx += 0.5) {
			structuredGrid.findSourcesForInterpolation(stage, 1, 1, 0, dx, src);
			for (int i = 0; i < IdealElastic2DNode::M; i++) {
				ASSERT_EQ(src[0](i), (i == 2 || i == 4) ? 1.0 : 0.0);
				ASSERT_EQ(src[1](i), 0.0);
			}
		}

	}
}


TEST(StructuredGrid, interpolateValuesAround)
{
	Task task;
	task.material = IsotropicMaterial(2.0, 2.0, 1.0);

	task.sizes(0) = task.sizes(1) = 3;
	task.accuracyOrder = 1;
	task.lengthes = {2, 2, 2}; // h_x = h_y = 1.0

	Task::InitialCondition::Quantity quantity;
	quantity.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	quantity.value = -1.0;
	quantity.area = std::make_shared<SphereArea>(0.1, linal::Vector3({1,1,0}));
	task.initialCondition.quantities.push_back(quantity);

	for (int stage = 0; stage <= 1; stage++) {

		StructuredGrid<IdealElastic2DNode> structuredGrid;
		structuredGrid.initialize(task);
		for (int x = 0; x < task.sizes(0); x++) {
			for (int y = 0; y < task.sizes(1); y++) {
				// check that values is set properly
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).u.V[0], 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).u.V[1], 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).u.sigma(0, 0), (x == 1 && y == 1) ? 1.0 : 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).u.sigma(0, 1), 0.0);
				ASSERT_EQ(structuredGrid.getNodeForTest(x, y, 0).u.sigma(1, 1), (x == 1 && y == 1) ? 1.0 : 0.0);
			}
		}

		IdealElastic2DNode::Vector dx({-1, 1, -0.5, 0.5, 0});
		IdealElastic2DNode::Matrix matrix = structuredGrid.interpolateValuesAround(stage, 1, 1, 0, dx);

		for (int i = 0; i < IdealElastic2DNode::M; i++) {
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
