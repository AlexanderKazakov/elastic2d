#include <gtest/gtest.h>

#include <lib/util/areas/areas.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/numeric/gcm/GcmHandler.hpp>
#include <test/wrappers/Wrappers.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

TEST(CubicGrid, initialize) {
	Task task;
	Statement statement;
	task.borderSize = 3;
	task.dimensionality = 2;
	statement.CourantNumber = 0.7;
	statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(4, 2, 0.5);
	task.sizes(0) = 7;
	task.sizes(1) = 9;
	task.lengthes = {20, 8, 1};
	statement.numberOfSnaps = 5;
	statement.T = 100.0;
	
	task.statements.push_back(statement);

	MeshWrapper<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>> grid(task);
	grid.beforeStatement(statement);
	ASSERT_NEAR(grid.h(0), 3.333333333, EQUALITY_TOLERANCE);
	ASSERT_NEAR(grid.h(1), 1.0, EQUALITY_TOLERANCE);
	ASSERT_NEAR(grid.getMaximalLambda(), 0.866025404, EQUALITY_TOLERANCE);
	ASSERT_NEAR(grid.getMinimalSpatialStep(), 1.0, EQUALITY_TOLERANCE);
	for (int x = 0; x < task.sizes(0); x++) {
		for (int y = 0; y < task.sizes(1); y++) {
			for (int z = 0; z < task.sizes(2); z++) {
				for (int i = 0; i < DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>::PdeVector::M; i++) {
					ASSERT_EQ(grid.pde({x, y, z})(i), 0.0);
				}
			}
		}
	}
}

TEST(CubicGrid, PartIterator) {
	Task task;

	task.dimensionality = 3;
	task.borderSize = 2;
	int X = 5, Y = 7, Z = 6;
	task.sizes = {X, Y, Z};
	task.lengthes = {20, 8, 1};
	
	MeshWrapper<DefaultMesh<Elastic3DModel, CubicGrid, IsotropicMaterial>> grid(task);
	
	int counter = 0;
	for (auto it = grid.slice((int)DIRECTION::Z, 3); it != it.end(); ++it) {
		counter++;
	}
	ASSERT_EQ(X * Y, counter);
	counter = 0;
	CubicGrid::Int3 min = {2, 2, 3}, max = {4, 3, 6};
	for (auto it = grid.box(min, max); it != it.end(); ++it) {
		counter++;
	}
	ASSERT_EQ(linal::directProduct(max-min), counter);
	auto partIter = grid.box({0, 0, 0}, grid.sizes);
	for (auto it : grid) {
		ASSERT_EQ(it, *partIter);
		++partIter;
	}
	for (int a = 1; a <= grid.borderSize; a++) {
		int direction = 0;
		int d = direction, d1 = (d+1)%3, d2 = (d+2)%3;
		auto realIter = grid.slice(direction,  a);
		auto virtIter = grid.slice(direction, -a);
		while (realIter != realIter.end()) {
			ASSERT_EQ(realIter(d1),  virtIter(d1));
			ASSERT_EQ(realIter(d2),  virtIter(d2));
			ASSERT_EQ(realIter(d), - virtIter(d));
			++realIter; ++virtIter;
		}
		ASSERT_EQ(virtIter.end(), virtIter);
	}
}


TEST(CubicGrid, interpolateValuesAround)
{
	Task task;
	Statement statement;
	statement.materialConditions.defaultMaterial = std::make_shared<IsotropicMaterial>(2, 2, 1);

	task.dimensionality = 2;
	task.sizes(0) = task.sizes(1) = 3;
	task.borderSize = 1;
	task.lengthes = {2, 2, 2}; // h_x = h_y = 1.0

	Statement::InitialCondition::Quantity quantity;
	quantity.physicalQuantity = PhysicalQuantities::T::PRESSURE;
	quantity.value = -1.0;
	quantity.area = std::make_shared<SphereArea>(0.1, linal::Vector3({1,1,0}));
	statement.initialCondition.quantities.push_back(quantity);

	task.statements.push_back(statement);
	
	for (int stage = 0; stage <= 1; stage++) {

		MeshWrapper<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>> grid(task);
		grid.beforeStatement(statement);
		for (int x = 0; x < task.sizes(0); x++) {
			for (int y = 0; y < task.sizes(1); y++) {
				// check that values is set properly
				ASSERT_EQ(grid.pde({x, y, 0}).V[0], 0.0);
				ASSERT_EQ(grid.pde({x, y, 0}).V[1], 0.0);
				ASSERT_EQ(grid.pde({x, y, 0}).sigma(0, 0), (x == 1 && y == 1) ? 1.0 : 0.0);
				ASSERT_EQ(grid.pde({x, y, 0}).sigma(0, 1), 0.0);
				ASSERT_EQ(grid.pde({x, y, 0}).sigma(1, 1), (x == 1 && y == 1) ? 1.0 : 0.0);
			}
		}

		DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>::PdeVector dx({-1, 1, -0.5, 0.5, 0});
		DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>::Iterator it({1, 1, 0}, grid.sizes);
		auto m = GcmHandler<Elastic2DModel, CubicGrid, IsotropicMaterial>::
				interpolateValuesAround(grid, stage, it, dx);

		for (int i = 0; i < DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial>::PdeVector::M; i++) {
			ASSERT_EQ(m(i, 0), 0.0) << "i = " << i; // Courant = 1
			ASSERT_EQ(m(i, 1), 0.0) << "i = " << i; // Courant = 1
			ASSERT_EQ(m(0, i), 0.0) << "i = " << i; // Vx
			ASSERT_EQ(m(1, i), 0.0) << "i = " << i; // Vy
			ASSERT_EQ(m(3, i), 0.0) << "i = " << i; // Sxy
		}

		ASSERT_EQ(m(2, 2), 0.5); // Courant = 0.5
		ASSERT_EQ(m(2, 3), 0.5); // Courant = 0.5
		ASSERT_EQ(m(4, 2), 0.5); // Courant = 0.5
		ASSERT_EQ(m(4, 3), 0.5); // Courant = 0.5

		ASSERT_EQ(m(2, 4), 1.0); // Courant = 0
		ASSERT_EQ(m(4, 4), 1.0); // Courant = 0
	}
}
