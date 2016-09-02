#include <gtest/gtest.h>

#include <libgcm/util/math/Area.hpp>
#include <libgcm/grid/cubic/CubicGrid.hpp>
#include <libgcm/engine/mesh/DefaultMesh.hpp>
#include <libgcm/rheology/models/models.hpp>


using namespace gcm;


TEST(CubicGrid, initialize) {
	Task task;
	task.materialConditions.byAreas.defaultMaterial = 
			std::make_shared<IsotropicMaterial>(4, 2, 0.5);
	
	typedef CubicGrid<2> Grid;
	typedef typename Grid::ConstructionPack ConstructionPack;
	ConstructionPack cp;
	cp.borderSize = 3;
	cp.sizes = {7, 9};
	cp.h = {1, 1};
	
	DefaultMesh<ElasticModel<2>, Grid, IsotropicMaterial> mesh(task, 0, cp);
	
	ASSERT_NEAR(mesh.getMaximalEigenvalue(), 0.866025404, EQUALITY_TOLERANCE);
	ASSERT_NEAR(mesh.getMinimalSpatialStep(), 1.0, EQUALITY_TOLERANCE);
	
	for (int x = 0; x < cp.sizes(0); x++) {
		for (int y = 0; y < cp.sizes(1); y++) {
			for (int i = 0; i < 5; i++) {
				ASSERT_EQ(mesh.pde({x, y})(i), 0.0);
			}
		}
	}
}


TEST(CubicGrid, PartIterator) {
	typedef CubicGrid<3> Grid;
	typedef typename Grid::ConstructionPack ConstructionPack;
	ConstructionPack cp;
	cp.borderSize = 2;
	int X = 5, Y = 7, Z = 6;
	cp.sizes = {X, Y, Z};
	cp.h = {1, 1, 1};
	Grid grid(0, cp);
	
	int counter = 0;
	for (auto it = grid.slice(2, 3); it != it.end(); ++it) {
		counter++;
	}
	ASSERT_EQ(X * Y, counter);
	ASSERT_EQ(grid.slice(2, 3).size(), counter);
	
	counter = 0;
	Int3 min = {2, 2, 3}, max = {4, 3, 6};
	for (auto it = grid.box(min, max); it != it.end(); ++it) {
		counter++;
	}
	ASSERT_EQ(linal::directProduct(max - min), counter);
	ASSERT_EQ(grid.box(min, max).size(), counter);
	
	auto partIter = grid.box({0, 0, 0}, grid.sizes);
	for (auto it : grid) {
		ASSERT_EQ(it, partIter);
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


