#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>
#include <libgcm/grid/simplex/SimplexGrid.hpp>
#include <libgcm/util/snapshot/VtkSnapshotter.hpp>

#include <gtest/gtest.h>

using namespace gcm;

typedef SimplexGrid<2, CgalTriangulation> Grid;
typedef typename Grid::Cell               Cell;
typedef typename Grid::Triangulation      Triangulation;
typedef typename Grid::Iterator           Iterator;
typedef typename Grid::RealD              RealD;


inline void testContains(const Grid& grid, const Cell& cell,
		const Iterator& it, const RealD& shift) {
	const RealD start = grid.coordsD(it);
	const RealD q = start + shift;
	if (cell.n == 3) {
		const RealD a = grid.coordsD(cell(0));
		const RealD b = grid.coordsD(cell(1));
		const RealD c = grid.coordsD(cell(2));
		ASSERT_TRUE(linal::triangleContains(a, b, c, q, grid.localEqualityTolerance()));
	} else if (!grid.isInner(it)) {
		ASSERT_EQ(0, cell.n);
	} else if (cell.n == 2) {
		const RealD a = grid.coordsD(cell(0));
		const RealD b = grid.coordsD(cell(1));
		RealD intersection = linal::linesIntersection(start, q, a, b);
		ASSERT_TRUE(linal::segmentContains(a, b, intersection, grid.localEqualityTolerance()));
	} else if (cell.n == 1) {
		const RealD a = grid.coordsD(cell(0));
		ASSERT_TRUE(linal::segmentContains(start, q, a, grid.localEqualityTolerance()));
	}
}


inline int testWholeGridOneDirection(
		const Grid& grid, const RealD& shift, const int m) {
	int hitCounter = 0;
	for (int i = 1; i < m; i++) {
		for (Iterator it : grid) {
			Cell c = grid.findCellCrossedByTheRay(it, shift * i);
			if (c.n > 0) { ++hitCounter; }
			testContains(grid, c, it, shift * i);
		}
	}
	return hitCounter;
}


TEST(LineWalkSearch, 2D) {
	Task task;
	real h = 0.5;
	task.simplexGrid.spatialStep = h;
	Task::SimplexGrid::Body::Border outer = {
		{3, 3}, {-3, 3}, {-3, -3}, {3, -3}, {2, 2}
	};
	Task::SimplexGrid::Body::Border inner = {
		{-2, -1}, {-1, 0}, {0, -1}, {-1, -2}
	};
	task.simplexGrid.bodies = {
		Task::SimplexGrid::Body({0, outer, {inner} }),
		Task::SimplexGrid::Body({0, {{-2, 5}, {2, 5}, {0, 7}}, {} })
	};
	Triangulation triangulation(task);
	Grid grid(0, {&triangulation});
	VtkUtils::dumpGridToVtk(grid);
	
	real step = h / 3;
	for (int i = 0; i < 16; i++) {
		real phi = i * M_PI / 8;
		RealD shift = step * RealD({cos(phi), sin(phi)});
		testWholeGridOneDirection(grid, shift, 10);
	}
}





