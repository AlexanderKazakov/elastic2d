#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>
#include <libgcm/grid/simplex/SimplexGrid.hpp>
#include <libgcm/util/snapshot/VtkSnapshotter.hpp>

#include <gtest/gtest.h>

using namespace gcm;

typedef SimplexGrid<2, CgalTriangulation> Grid;
typedef typename Grid::Cell               Cell;
typedef typename Grid::Triangulation      Triangulation;
typedef typename Grid::Iterator           Iterator;
typedef typename Grid::BorderIterator     BorderIter;
typedef typename Grid::RealD              RealD;
typedef elements::Element<Real2, 3>       RealCell;

#define CellIterToCellRealD(gridName) [&](Iterator iter) {return  gridName.coordsD(iter);}


inline void testContains(const Grid& grid, const RealCell& cell,
		const Iterator& it, const RealD& shift, int& hitCounter) {
	const RealD start = grid.coordsD(it);
	const RealD q = start + shift;
	if (cell.n == 3) {
		ASSERT_TRUE(linal::triangleContains(
				cell(0), cell(1), cell(2), q, EQUALITY_TOLERANCE));
		++hitCounter;
	} else if (!grid.isInner(it)) {
		ASSERT_EQ(0, cell.n);
	} else if (cell.n == 2) {
		RealD intersection = linal::linesIntersection(start, q, cell(0), cell(1));
		ASSERT_TRUE(linal::segmentContains(cell(0), cell(1), intersection,
				EQUALITY_TOLERANCE, grid.localEqualityTolerance()));
		++hitCounter;
	} else if (cell.n == 1) {
		ASSERT_TRUE(linal::segmentContains(start, q, cell(0), EQUALITY_TOLERANCE));
		++hitCounter;
	}
}


inline void testContains(const Grid& grid, const Cell& cell,
		const Iterator& it, const RealD& shift, int& hitCounter) {
	RealCell rc(cell, CellIterToCellRealD(grid));
	testContains(grid, rc, it, shift, hitCounter);
}


inline void checkBothCellsContainQueryPoint(
		const RealCell a, const RealCell b, const RealD q, const real eps) {
	ASSERT_EQ(a.n, b.n);
	std::set<RealD> common = a.equalPoints(b);
	if (common.size() == 3) {
		ASSERT_TRUE(linal::triangleContains(a(0), a(1), a(2), q, EQUALITY_TOLERANCE));
	} else if (common.size() == 2) {
		ASSERT_TRUE(linal::segmentContains(
				*common.begin(), *common.rbegin(), q, EQUALITY_TOLERANCE, eps));
	} else if (common.size() == 1) {
		ASSERT_LT(linal::length(*common.begin() - q), eps);
	}
}


inline void matchSearchResults(const Grid& grid,
		const Cell& byLineWalk, const Cell& byCgal, 
		const Iterator& it, const RealD& shift,
		const bool gridIsConvex, int& hitCounter) {
	const RealD start = grid.coordsD(it);
	const RealD q = start + shift;
	if (byLineWalk.n == 3 && byCgal.n == 3) {
		RealCell a(byLineWalk, CellIterToCellRealD(grid));
		RealCell b(byCgal,     CellIterToCellRealD(grid));
		checkBothCellsContainQueryPoint(a, b, q, grid.localEqualityTolerance());
	} else if (gridIsConvex) {
		ASSERT_EQ(0, byCgal.n) << "byLineWalk.n == " << byLineWalk.n;
	}
	testContains(grid, byLineWalk, it, shift, hitCounter);
}


template<typename Predicate>
inline void testWholeGridOneDirection(
		const Grid& grid, const RealD& shift, const int m,
		const Predicate test, int& hitCounter) {
	hitCounter = 0;
	for (int i = 1; i < m; i++) {
		for (Iterator it : grid) {
			test(it, shift * i, hitCounter);
		}
	}
}


TEST(LineWalkSearch2D, VersusLinal) {
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
	
	real step = h / 3;
	for (int i = 0; i < 16; i++) {
		real phi = i * M_PI / 8;
		RealD direction = step * RealD({cos(phi), sin(phi)});
		int hitCount = 0;
		testWholeGridOneDirection(grid, direction, 10,
				[&](Iterator it, RealD shift, int& hitCounter) {
					Cell c = grid.findCellCrossedByTheRay(it, shift);
					testContains(grid, c, it, shift, hitCounter);
				}, hitCount);
		ASSERT_NEAR(3100, hitCount, 300);
	}
}


TEST(LineWalkSearch2D, VersusCgalSearchFunction) {
	Task task;
	real h = 0.5;
	task.simplexGrid.spatialStep = h;
	// For convex cases CGAL search and ours line walk search must match each other
	Task::SimplexGrid::Body::Border outer = {
		{3, 3}, {-3, 3}, {-3, -3}, {3, -3}, {5, 1}
	};
	task.simplexGrid.bodies = {
		Task::SimplexGrid::Body({ 0, outer, { } })
	};
	Triangulation triangulation(task);
	Grid grid(0, {&triangulation});
	
	real step = h / 3;
	for (int i = 0; i < 16; i++) {
		real phi = i * M_PI / 8;
		RealD direction = step * RealD({cos(phi), sin(phi)});
		int hitCount = 0;
		testWholeGridOneDirection(grid, direction, 10,
				[&](Iterator it, RealD shift, int& hitCounter) {
					Cell lw = grid.findCellCrossedByTheRay(it, shift);
					Cell cgal = grid.locateOwnerCell(it, shift);
					matchSearchResults(grid, lw, cgal, it, shift, true, hitCounter);
				}, hitCount);
		ASSERT_NEAR(3750, hitCount, 300);
	}
}


TEST(LineWalkSearch2D, findCellCrossedByTheRayNotConvex) {
	Task task;
	real h = 5;
	task.simplexGrid.spatialStep = h;
	// triangulation of the first body will be the same in "one" and "aFew"
	task.simplexGrid.bodies = {
		Task::SimplexGrid::Body({ 0, {{3, 3}, {-3, 3}, {-3, -3}, {3, -3}}, { } }),
	};
	Triangulation triangulation(task);
	Grid one(0, {&triangulation});
	
	task.simplexGrid.bodies.push_back(Task::SimplexGrid::Body(
			{ 0, {{-2, 5}, {2, 5}, {0, 7}}, { } }));
	task.simplexGrid.bodies.push_back(Task::SimplexGrid::Body(
			{ 0, {{-7, -1}, {-5, -4},  {-2, -5},  {2, -5}, {0, -7}, {-7, -7}}, { } }));
	task.simplexGrid.bodies.push_back(Task::SimplexGrid::Body(
			{ 0, {{-10, -10}, {-10, 10}, {10, 10}, {10, -10}},
					{ {{-9, -9}, {-9, 9}, {9, 9}, {9, -9}} } }));
	Triangulation triangulation2(task);
	Grid aFew(0, {&triangulation2});
	
	real step = h / 3;
	for (int i = 0; i < 16; i++) {
		real phi = i * M_PI / 8;
		RealD direction = step * RealD({cos(phi), sin(phi)});
		int hitCount = 0;
		testWholeGridOneDirection(one, direction, 10,
				[&](Iterator it, RealD shift, int& hitCounter) {
					RealD start = one.coordsD(it);
					RealD query = start + shift;
					Iterator it2 = aFew.findVertexByCoordinates(start);
					Cell a =  one.findCellCrossedByTheRay(it,  shift);
					Cell b = aFew.findCellCrossedByTheRay(it2, shift);
					RealCell ar(a, CellIterToCellRealD(one));
					RealCell br(b, CellIterToCellRealD(aFew));
					RealD q = query;
					if (ar.n == 2) { q = linal::linesIntersection(ar(0), ar(1), start, query); }
					if (ar.n == 1) { q = ar(0); }
					checkBothCellsContainQueryPoint(ar, br, q, one.localEqualityTolerance());
					testContains(one, a, it, shift, hitCounter);
				}, hitCount);
		ASSERT_NEAR(18, hitCount, 4);
	}
}


TEST(LineWalkSearch2D, CasesAlongBorder) {
	Task task;
	real h = 0.1, step = h / 3;
	task.simplexGrid.spatialStep = h;
	Task::SimplexGrid::Body::Border border = {{-1, -1}, {-1, 1}, {1, 1}, {1, -1}};
	task.simplexGrid.bodies = {Task::SimplexGrid::Body({0, border, {}})};
	Triangulation triangulation(task);
	Grid grid(0, {&triangulation});
	
	auto check = [&](RealD a, RealD b, RealD shift, int& hitCounter) {
		for (Iterator it : grid) {
			RealD start = grid.coordsD(it);
			if (!linal::segmentContains(a, b, start,
					EQUALITY_TOLERANCE, grid.localEqualityTolerance())) { continue; }
			RealD query = start + shift;
			Cell c = grid.findCellCrossedByTheRay(it, shift);
			if (linal::segmentContains(a, b, query,
					EQUALITY_TOLERANCE, grid.localEqualityTolerance())) {
				ASSERT_EQ(3, c.n);
				testContains(grid, c, it, shift, hitCounter);
			} else {
				ASSERT_EQ(0, c.n);
			}
		}
	};
	
	int hitCountPrev = 0;
	for (size_t i = 0; i < border.size(); i++) {
		int hitCount = 0;
		Task::SimplexGrid::Body::Point pA = border[i];
		Task::SimplexGrid::Body::Point pB = border[(i + 1) % border.size()];
		RealD a = {pA[0], pA[1]}, b = {pB[0], pB[1]};
		RealD direction = linal::normalize(b - a);
		for (int j = 0; j < 10; j++) {
			RealD shift = direction * step * j;
			check(a, b,  shift, hitCount);
			check(a, b, -shift, hitCount);
		}
		if (i == 0) { hitCountPrev = hitCount; }
		else { ASSERT_EQ(hitCountPrev, hitCount); }
	}
	ASSERT_NEAR(600, hitCountPrev, 10);
}


#undef CellIterToCellRealD
