#include <gtest/gtest.h>

#include <lib/util/Area.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <test/wrappers/Wrappers.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

TEST(CubicGrid, initialize) {
	Task task;
	Statement statement;
	task.cubicGrid.borderSize = 3;
	task.cubicGrid.dimensionality = 2;
	task.cubicGrid.sizes = {7, 9, 1};
	task.cubicGrid.lengths = {20, 8, 1};
	statement.materialConditions.defaultMaterial = 
			std::make_shared<IsotropicMaterial>(4, 2, 0.5);
	task.statements.push_back(statement);

	CubicGrid::preprocessTask(task.cubicGrid);
	MeshWrapper<DefaultMesh<Elastic2DModel, CubicGrid, IsotropicMaterial> > grid(task);
	grid.beforeStatementForTest(statement);
	
	ASSERT_NEAR(grid.h(0), 3.333333333, EQUALITY_TOLERANCE);
	ASSERT_NEAR(grid.h(1), 1.0, EQUALITY_TOLERANCE);
	ASSERT_NEAR(grid.getMaximalEigenvalue(), 0.866025404, EQUALITY_TOLERANCE);
	ASSERT_NEAR(grid.getMinimalSpatialStep(), 1.0, EQUALITY_TOLERANCE);
	
	for (int x = 0; x < task.cubicGrid.sizes(0); x++) {
		for (int y = 0; y < task.cubicGrid.sizes(1); y++) {
			for (int z = 0; z < task.cubicGrid.sizes(2); z++) {
				for (int i = 0; i < 5; i++) {
					ASSERT_EQ(grid.pde({x, y, z}) (i), 0.0);
				}
			}
		}
	}
}


TEST(CubicGrid, PartIterator) {
	Task task;

	task.cubicGrid.dimensionality = 3;
	task.cubicGrid.borderSize = 2;
	int X = 5, Y = 7, Z = 6;
	task.cubicGrid.sizes = {X, Y, Z};
	task.cubicGrid.lengths = {20, 8, 1};

	CubicGrid::preprocessTask(task.cubicGrid);
	MeshWrapper<DefaultMesh<Elastic3DModel, CubicGrid, IsotropicMaterial> > grid(task);

	int counter = 0;
	for (auto it = grid.slice((int)DIRECTION::Z, 3); it != it.end(); ++it) {
		counter++;
	}
	ASSERT_EQ(X * Y, counter);
	counter = 0;
	Int3 min = {2, 2, 3}, max = {4, 3, 6};
	for (auto it = grid.box(min, max); it != it.end(); ++it) {
		counter++;
	}
	ASSERT_EQ(linal::directProduct(max - min), counter);
	auto partIter = grid.box({0, 0, 0}, grid.sizes);
	for (auto it : grid) {
		ASSERT_EQ(it, *partIter);
		++partIter;
	}
	for (int a = 1; a <= grid.borderSize; a++) {
		int direction = 0;
		int d = direction, d1 = (d + 1) % 3, d2 = (d + 2) % 3;
		auto realIter = grid.slice(direction, a);
		auto virtIter = grid.slice(direction, -a);
		while (realIter != realIter.end()) {
			ASSERT_EQ(realIter(d1), virtIter(d1));
			ASSERT_EQ(realIter(d2), virtIter(d2));
			ASSERT_EQ(realIter(d), -virtIter(d));
			++realIter; ++virtIter;
		}
		ASSERT_EQ(virtIter.end(), virtIter);
	}
}


