#include <gtest/gtest.h>

#include <lib/util/areas/SphereArea.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/numeric/gcmethod/GcmHandler.hpp>
#include <test/wrappers/Wrappers.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

TEST(StructuredGrid, initialize) {
	Task task;
	task.accuracyOrder = 3;
	task.CourantNumber = 0.7;
	task.isotropicMaterial = IsotropicMaterial(4.0, 2.0, 0.5);
	task.sizes(0) = 7;
	task.sizes(1) = 9;
	task.lengthes = {20, 8, 1};
	task.numberOfSnaps = 5;
	task.T = 100.0;

	MeshWrapper<DefaultGrid<Elastic2DModel, StructuredGrid>> grid;
	grid.initialize(task);
	ASSERT_NEAR(grid.getH()(0), 3.333333333, EQUALITY_TOLERANCE);
	ASSERT_NEAR(grid.getH()(1), 1.0, EQUALITY_TOLERANCE);
	ASSERT_NEAR(grid.getMaximalLambda(), 0.866025404, EQUALITY_TOLERANCE);
	ASSERT_NEAR(grid.getMinimalSpatialStep(), 1.0, EQUALITY_TOLERANCE);
	for (int x = 0; x < task.sizes(0); x++) {
		for (int y = 0; y < task.sizes(1); y++) {
			for (int z = 0; z < task.sizes(2); z++) {
				for (int i = 0; i < DefaultGrid<Elastic2DModel, StructuredGrid>::PdeVector::M; i++) {
					ASSERT_EQ(grid.getPde(x, y, z)(i), 0.0);
				}
			}
		}
	}
}


TEST(StructuredGrid, interpolateValuesAround)
{
	Task task;
	task.isotropicMaterial = IsotropicMaterial(2.0, 2.0, 1.0);

	task.sizes(0) = task.sizes(1) = 3;
	task.accuracyOrder = 1;
	task.lengthes = {2, 2, 2}; // h_x = h_y = 1.0

	Task::InitialCondition::Quantity quantity;
	quantity.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	quantity.value = -1.0;
	quantity.area = std::make_shared<SphereArea>(0.1, linal::Vector3({1,1,0}));
	task.initialCondition.quantities.push_back(quantity);

	for (int stage = 0; stage <= 1; stage++) {

		MeshWrapper<DefaultGrid<Elastic2DModel, StructuredGrid>> grid;
		grid.initialize(task);
		for (int x = 0; x < task.sizes(0); x++) {
			for (int y = 0; y < task.sizes(1); y++) {
				// check that values is set properly
				ASSERT_EQ(grid.getPde(x, y, 0).V[0], 0.0);
				ASSERT_EQ(grid.getPde(x, y, 0).V[1], 0.0);
				ASSERT_EQ(grid.getPde(x, y, 0).sigma(0, 0), (x == 1 && y == 1) ? 1.0 : 0.0);
				ASSERT_EQ(grid.getPde(x, y, 0).sigma(0, 1), 0.0);
				ASSERT_EQ(grid.getPde(x, y, 0).sigma(1, 1), (x == 1 && y == 1) ? 1.0 : 0.0);
			}
		}

		DefaultGrid<Elastic2DModel, StructuredGrid>::PdeVector dx({-1, 1, -0.5, 0.5, 0});
		DefaultGrid<Elastic2DModel, StructuredGrid>::Iterator it({1, 1, 0}, grid.getSizes());
		DefaultGrid<Elastic2DModel, StructuredGrid>::Matrix matrix =
				GcmHandler<Elastic2DModel, StructuredGrid>::interpolateValuesAround(grid, stage, it, dx);

		for (int i = 0; i < DefaultGrid<Elastic2DModel, StructuredGrid>::PdeVector::M; i++) {
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
