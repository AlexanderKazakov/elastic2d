#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>
#include <libgcm/grid/simplex/SimplexGrid.hpp>
#include <libgcm/util/snapshot/VtkSnapshotter.hpp>

#include <gtest/gtest.h>

using namespace gcm;

typedef SimplexGrid<3, CgalTriangulation> Grid;
typedef typename Grid::Cell               Cell;
typedef typename Grid::Triangulation      Triangulation;
typedef typename Grid::Iterator           Iterator;
typedef typename Grid::BorderIterator     BorderIter;
typedef typename Grid::RealD              RealD;
typedef elements::Element<Real3, 4>       RealCell;

#define CellIterToCellReal(gridName) [&](Iterator iter) {return  gridName.coordsD(iter);}


inline void testContains(const Grid& grid, const RealCell& cell,
		const Iterator& it, const RealD& shift, int& hitCounter) {
	const RealD start = grid.coordsD(it);
	const RealD q = start + shift;
	if (cell.n == 4) {
		ASSERT_TRUE(linal::tetrahedronContains(
				cell(0), cell(1), cell(2), cell(3), q, EQUALITY_TOLERANCE));
		++hitCounter;
	} else if (!grid.isInner(it)) {
		ASSERT_EQ(0, cell.n);
	} else if (cell.n == 3) {
		RealD intersection = linal::lineWithFlatIntersection(
				cell(0), cell(1), cell(2), start, q);
		ASSERT_TRUE(linal::triangleContains(cell(0), cell(1), cell(2),
				intersection, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
		++hitCounter;
	} else if (cell.n == 2) {
		RealD intersection = linal::linesIntersection(cell(0), cell(1), start, q);
		ASSERT_TRUE(linal::segmentContains(cell(0), cell(1),
				intersection, EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	} else if (cell.n == 1) {
		ASSERT_TRUE(linal::segmentContains(
				start, q, cell(0), EQUALITY_TOLERANCE, EQUALITY_TOLERANCE));
	} else if (cell.n == 0) {
		ASSERT_FALSE(grid.isInner(it));
	} else {
		FAIL() << "Impossible value for cell.n: " << cell.n;
	}
}


inline void testContains(const Grid& grid, const Cell& cell,
		const Iterator& it, const RealD& shift, int& hitCounter) {
	RealCell rc(cell, CellIterToCellReal(grid));
	testContains(grid, rc, it, shift, hitCounter);
}


inline void checkBothCellsContainQueryPoint(
		const RealCell a, const RealCell b, const RealD q, const real eps) {
	ASSERT_EQ(a.n, b.n);
	std::set<RealD> common = a.equalPoints(b);
	std::vector<RealD> p(common.begin(), common.end());
	ASSERT_GT(p.size(), 0);
	if (p.size() == 4) {
		ASSERT_TRUE(linal::tetrahedronContains(
				p[0], p[1], p[2], p[3], q, EQUALITY_TOLERANCE));
	} else if (p.size() == 3) {
		if (linal::triangleContains(
				p[0], p[1], p[2], q, EQUALITY_TOLERANCE, eps)) { return; }
		if (linal::segmentContains(
				p[0], p[1], q, EQUALITY_TOLERANCE, eps)) { return; }
		if (linal::segmentContains(
				p[0], p[2], q, EQUALITY_TOLERANCE, eps)) { return; }
		ASSERT_TRUE(linal::segmentContains(
				p[1], p[2], q, EQUALITY_TOLERANCE, eps));
	} else if (p.size() == 2) {
		ASSERT_TRUE(linal::segmentContains(
				*p.begin(), *p.rbegin(), q, EQUALITY_TOLERANCE, eps));
	} else if (p.size() == 1) {
		ASSERT_LT(linal::length(q - *p.begin()), eps);
	}
}


inline void matchSearchResults(const Grid& grid,
		const Cell& byLineWalk, const Cell& byCgal, 
		const Iterator& it, const RealD& shift) {
	const RealD start = grid.coordsD(it);
	const RealD q = start + shift;
	if (byLineWalk.n == 4 && byCgal.n == 4) {
		RealCell a(byLineWalk, CellIterToCellReal(grid));
		RealCell b(byCgal,     CellIterToCellReal(grid));
		checkBothCellsContainQueryPoint(a, b, q, EQUALITY_TOLERANCE);
	}
}


template<typename Predicate>
inline void testWholeGridOneDirection(
		const Grid& grid, const RealD& shift, const int m,
		const Predicate test, int& hitCounter) {
	hitCounter = 0;
	for (int i = 1; i < m; i++) {
		for (Iterator it : grid) {
			RealD dir = shift * i;
			try {
				test(it, dir, hitCounter);
			} catch (Exception& e) {
				RealD start = grid.coordsD(it);
				std::cout << e.what() << std::endl 
						<< "start = " << start << "query = " << start + dir;
				throw;
			}
		}
	}
}


inline void test3DFigure(
		const std::string filename, const real h, const bool sharpEdges,
		const int hitCountMin) {
	std::cout << "Start testing grid from file " << filename
			<< " with h == " << h << std::endl;
	Task task;
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::CGAL_MESHER;
	task.simplexGrid.spatialStep = h;
	task.simplexGrid.detectSharpEdges = sharpEdges;
	task.simplexGrid.fileName = filename;
	Triangulation triangulation(task);
	Grid grid(0, {&triangulation});
	real step = task.simplexGrid.spatialStep / 3;
	for (int i = 0; i < 16; i++) {
		real phi = i * M_PI / 8;
		for (int j = 0; j < 16; j++) {
			real teta = j * M_PI / 8;
			RealD direction = step * RealD({
					cos(phi) * cos(teta), sin(phi) * cos(teta), sin(teta)});
			int hitCount = 0;
			testWholeGridOneDirection(grid, direction, 10,
					[&](Iterator it, RealD shift, int& hitCounter) {
						Cell c = grid.findCellCrossedByTheRay(it, shift);
						testContains(grid, c, it, shift, hitCounter);
					}, hitCount);
			std::cout << "hitCount == " << hitCount << std::endl;
			ASSERT_GT(hitCount, hitCountMin);
			testWholeGridOneDirection(grid, direction, 10,
					[&](Iterator it, RealD shift, int&) {
						Cell lw = grid.findCellCrossedByTheRay(it, shift);
						Cell cgal = grid.locateOwnerCell(it, shift);
						matchSearchResults(grid, lw, cgal, it, shift);
					}, hitCount);
		}
	}
}
/*
TEST(LineWalkSearch3D, VersusLinalAndCgal) {
	test3DFigure("meshes/tetrahedron.off", 0.6, true, 20);
	test3DFigure("meshes/cube.off", 0.4, true, 50);
	test3DFigure("meshes/icosahedron.off", 1.0, false, 60);
	
	test3DFigure("meshes/tetrahedron.off", 0.4, true, 150);
	test3DFigure("meshes/cube.off", 0.2, true, 700);
	test3DFigure("meshes/icosahedron.off", 0.4, false, 1700);
	
//	test3DFigure("meshes/cube.off", 0.025, true, 400000);
}


inline bool quadrateContains(
		RealD a, RealD b, RealD c, RealD d, RealD q, const real eps) {
	assert_true(linal::isDegenerate(a, b, c, d, eps));
	return  linal::triangleContains(a, b, c, q, EQUALITY_TOLERANCE, eps) ||
			linal::triangleContains(a, b, d, q, EQUALITY_TOLERANCE, eps) ||
			linal::triangleContains(a, c, d, q, EQUALITY_TOLERANCE, eps) ||
			linal::triangleContains(c, b, d, q, EQUALITY_TOLERANCE, eps);
}

TEST(LineWalkSearch3D, CasesAlongBorder) {
	Task task;
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::CGAL_MESHER;
//	real h = 0.1, step = h / 3;
	real h = 0.025, step = h / 3;
	task.simplexGrid.spatialStep = h;
	task.simplexGrid.detectSharpEdges = true;
	task.simplexGrid.fileName = "meshes/cube.off";
	
	try {
	
	Triangulation triangulation(task);
	Grid grid(0, {&triangulation});
	// check that cube geometry is not changed in the file
	ASSERT_NO_THROW(grid.findVertexByCoordinates({0, 0, 0}));
	ASSERT_NO_THROW(grid.findVertexByCoordinates({1, 1, 1}));
	
	auto check = [&](RealD a, RealD b, RealD c, RealD d,
			Iterator it, RealD shift, int& hitCounter) {
		RealD start = grid.coordsD(it);
		RealD query = start + shift;
		Cell iterCell;
		if (!quadrateContains(
				a, b, c, d, start, EQUALITY_TOLERANCE)) { return; }
		try {
			iterCell = grid.findCellCrossedByTheRay(it, shift);
		} catch (Exception& e) {
			std::cout << e.what(); throw;
		}
		RealCell cell(iterCell, CellIterToCellReal(grid));
		if (quadrateContains(a, b, c, d, query, EQUALITY_TOLERANCE)) {
			ASSERT_EQ(4, cell.n) << "start == " << start << "query == " << query << std::endl;
			ASSERT_TRUE(linal::tetrahedronContains(
					cell(0), cell(1), cell(2), cell(3), query, EQUALITY_TOLERANCE));
			++hitCounter;
		} else {
			ASSERT_EQ(0, cell.n);
		}
	};
	
	std::vector<std::vector<RealD>> facets = {
		{RealD({0,0,0}), RealD({1,0,0}), RealD({0,1,0}), RealD({1,1,0})},
		{RealD({0,0,0}), RealD({0,0,1}), RealD({0,1,0}), RealD({0,1,1})},
		{RealD({0,0,0}), RealD({0,0,1}), RealD({1,0,0}), RealD({1,0,1})},
		{RealD({1,1,1}), RealD({0,1,1}), RealD({1,0,1}), RealD({0,0,1})},
		{RealD({1,1,1}), RealD({1,1,0}), RealD({1,0,1}), RealD({1,0,0})},
		{RealD({1,1,1}), RealD({1,1,0}), RealD({0,1,1}), RealD({0,1,0})},
	};
	
	for (std::vector<RealD> facet : facets) {
		RealD a = facet[0], b = facet[1], c = facet[2], d = facet[3];
		int hitCount = 0;
		for (size_t i = 0; i < 4; i++) {
			for (size_t j = 0; j < 4; j++) if (j != i) {
				RealD direction = linal::normalize(facet[i] - facet[j]);
				for (int m = 0; m < 10; m++) {
					RealD shift = direction * step * m;
					for (BorderIter it = grid.borderBegin(); it != grid.borderEnd(); ++it) {
						check(a, b, c, d, *it, shift, hitCount);
					}
				}
			}
		}
		ASSERT_GT(hitCount, 10000);
	}
	
	} catch (Exception& e) {
		std::cout << e.what(); throw;
	}
}
*/

TEST(LineWalkSearch3D, Skull) {
	Task task;
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::INM_MESHER;
	task.simplexGrid.fileName = "meshes/coarse/mesh-aneurysm.out";
	task.simplexGrid.scale = 10;
	try {
	
	Triangulation triangulation(task);
	Grid grid(5, {&triangulation});
	vtk_utils::dumpGridToVtk(grid);
	real h = grid.getAverageHeight(), step = h / 3;
	
	for (int i = 0; i < 16; i++) {
		real phi = i * M_PI / 8;
		for (int j = 0; j < 16; j++) {
			real teta = j * M_PI / 8;
			RealD direction = step * RealD({
					cos(phi) * cos(teta), sin(phi) * cos(teta), sin(teta)});
			int hitCount = 0;
			testWholeGridOneDirection(grid, direction, 10,
					[&](Iterator it, RealD shift, int& hitCounter) {
						Cell c = grid.findCellCrossedByTheRay(it, shift);
						testContains(grid, c, it, shift, hitCounter);
					}, hitCount);
			std::cout << "i == " << i << "j == " << j << "hitCount == " << hitCount << std::endl;
			ASSERT_GT(hitCount, 130000);
		}
	}
	
	} catch (Exception& e) {
		std::cout << e.what(); throw;
	}
}


#undef CellIterToCellReal
