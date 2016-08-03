#include <gtest/gtest.h>

#include <lib/util/Area.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/mesh/DefaultMesh.hpp>
#include <lib/rheology/models/models.hpp>

using namespace gcm;

TEST(CubicGrid, initialize) {
	Task task;
	Statement statement;
	task.cubicGrid.borderSize = 3;
	task.cubicGrid.sizes = {7, 9};
	task.cubicGrid.lengths = {20, 8};
	statement.materialConditions.byAreas.defaultMaterial = 
			std::make_shared<IsotropicMaterial>(4, 2, 0.5);
	task.statements.push_back(statement);
	
	typedef CubicGrid<2> Grid;
	typedef typename Grid::GlobalScene GS;
	std::shared_ptr<GS> gs(new GS(task));
	DefaultMesh<ElasticModel<2>, Grid, IsotropicMaterial> mesh(task, gs.get(), 0);
	mesh.beforeStatement(statement);
	
	ASSERT_NEAR(mesh.h(0), 3.333333333, EQUALITY_TOLERANCE);
	ASSERT_NEAR(mesh.h(1), 1.0, EQUALITY_TOLERANCE);
	ASSERT_NEAR(mesh.getMaximalEigenvalue(), 0.866025404, EQUALITY_TOLERANCE);
	ASSERT_NEAR(mesh.getMinimalSpatialStep(), 1.0, EQUALITY_TOLERANCE);
	
	for (int x = 0; x < task.cubicGrid.sizes.at(0); x++) {
		for (int y = 0; y < task.cubicGrid.sizes.at(1); y++) {
			for (int i = 0; i < 5; i++) {
				ASSERT_EQ(mesh.pde({x, y})(i), 0.0);
			}
		}
	}
}


TEST(CubicGrid, PartIterator) {
	Task task;
	
	task.cubicGrid.borderSize = 2;
	int X = 5, Y = 7, Z = 6;
	task.cubicGrid.sizes = {X, Y, Z};
	task.cubicGrid.lengths = {20, 8, 1};
	
	typedef CubicGrid<3> Grid;
	typedef typename Grid::GlobalScene GS;
	std::shared_ptr<GS> gs(new GS(task));
	DefaultMesh<ElasticModel<3>, Grid, IsotropicMaterial> mesh(task, gs.get(), 0);
	
	int counter = 0;
	for (auto it = mesh.slice((int)DIRECTION::Z, 3); it != it.end(); ++it) {
		counter++;
	}
	ASSERT_EQ(X * Y, counter);
	counter = 0;
	Int3 min = {2, 2, 3}, max = {4, 3, 6};
	for (auto it = mesh.box(min, max); it != it.end(); ++it) {
		counter++;
	}
	ASSERT_EQ(linal::directProduct(max - min), counter);
	auto partIter = mesh.box({0, 0, 0}, mesh.sizes);
	for (auto it : mesh) {
		ASSERT_EQ(it, partIter);
		++partIter;
	}
	for (int a = 1; a <= mesh.borderSize; a++) {
		int direction = 0;
		int d = direction, d1 = (d + 1) % 3, d2 = (d + 2) % 3;
		auto realIter = mesh.slice(direction, a);
		auto virtIter = mesh.slice(direction, -a);
		while (realIter != realIter.end()) {
			ASSERT_EQ(realIter(d1), virtIter(d1));
			ASSERT_EQ(realIter(d2), virtIter(d2));
			ASSERT_EQ(realIter(d), -virtIter(d));
			++realIter; ++virtIter;
		}
		ASSERT_EQ(virtIter.end(), virtIter);
	}
}


