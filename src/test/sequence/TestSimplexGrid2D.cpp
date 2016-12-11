#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>
#include <libgcm/grid/simplex/SimplexGrid.hpp>
#include <libgcm/util/snapshot/VtkSnapshotter.hpp>

#include <gtest/gtest.h>


using namespace gcm;


TEST(SimplexGrid2D, ownerTriangleVsBarycentric) {
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
	
	typedef SimplexGrid<2, CgalTriangulation> Grid;
	typedef typename Grid::Triangulation Triangulation;
	Triangulation triangulation(task);
	Grid grid(0, {&triangulation});
	VtkUtils::dumpGridToVtk(grid);
	
	Utils::seedRand();
	for (int multiplier = 1; multiplier < 10; multiplier++) {
		int cntXY = 0/*, cntX = 0, cntY = 0*/;
		real step = h / 3 * multiplier;
		for (auto& it : grid) {
//			if (linal::length(grid.coordsD(it) - Real2({1.93934, 3})) > 0.01) { continue; } //!
			for (int test_counter = 0; test_counter < 8; test_counter++) {
				Real2 shift;
				auto checkInnerTriangle = [&](int& counter) {
					typename Grid::Cell triangleF;
					typename Grid::Cell triangleL;
//					try {
						triangleF = grid.findCellCrossedByTheRay(it, shift);
						triangleL = grid.locateOwnerCell(it, shift);
//					} catch (Exception e) {
//						std::cout << "shift == " << shift << "step = " << step << std::endl
//								<< e.what() << std::endl;
//					}
					if (triangleF.n == 3) {
						auto lambda = linal::barycentricCoordinates(
								grid.coordsD(triangleF(0)), 
								grid.coordsD(triangleF(1)), 
								grid.coordsD(triangleF(2)), 
								grid.coordsD(it) + shift);
						for (int i = 0; i < 3; i++) {
							ASSERT_TRUE(lambda(i) > -EQUALITY_TOLERANCE) << lambda(i);
						}
						counter++;
						
					} else if (triangleF.n == 2) {
						Real2 cross = linal::linesIntersection(
								grid.coordsD(triangleF(0)), 
								grid.coordsD(triangleF(1)),
								grid.coordsD(it),
								grid.coordsD(it) + shift);
						ASSERT_LT(fabs(linal::orientedArea(
								grid.coordsD(triangleF(0)), 
								grid.coordsD(triangleF(1)),
								cross)), EQUALITY_TOLERANCE);
					}
					
					if (triangleL.n == 3) {
						auto lambda = linal::barycentricCoordinates(
								grid.coordsD(triangleL(0)), 
								grid.coordsD(triangleL(1)), 
								grid.coordsD(triangleL(2)), 
								grid.coordsD(it) + shift);
						for (int i = 0; i < 3; i++) {
							ASSERT_TRUE(lambda(i) > -EQUALITY_TOLERANCE) << lambda(i);
						}
					}
				};
				
				shift = {Utils::randomReal(-step, step), Utils::randomReal(-step, step)};
//				shift = { 0.499485, -0.586246 }; //!
				checkInnerTriangle(cntXY);
				
//				shift = {Utils::randomReal(-step, step), 0};
//				checkInnerTriangle(cntX);
				
//				shift = {0, Utils::randomReal(-step, step)};
//				checkInnerTriangle(cntY);
			}
		}
		
//		ASSERT_TRUE(Utils::approximatelyEqual(2400, cntXY, 0.3)) << "Actual: " << cntXY;
//		ASSERT_TRUE(Utils::approximatelyEqual(2400, cntX, 0.3)) << "Actual: " << cntX;
//		ASSERT_TRUE(Utils::approximatelyEqual(2400, cntY, 0.3)) << "Actual: " << cntY;
	}
}


TEST(SimplexGrid2D, locateVsFindOwnerTriangle) {
	Task task;
	real h = 0.5;
	task.simplexGrid.spatialStep = h;
	
	Task::SimplexGrid::Body::Border outer = {
		{3, 3}, {-3, 3}, {-3, -3}, {3, -3}, {5, 1}
	};
	task.simplexGrid.bodies = {
		Task::SimplexGrid::Body({ 0, outer, { } })
	};
	
	typedef SimplexGrid<2, CgalTriangulation> Grid;
	typedef typename Grid::Triangulation Triangulation;
	Triangulation triangulation(task);
	Grid grid(0, {&triangulation});
	
	
	Utils::seedRand();
	for (int multiplier = 1; multiplier < 10; multiplier++) {
		real step = h / 3 * multiplier;
		for (auto& it : grid) {
			for (int test_counter = 0; test_counter < 8; test_counter++) {
			
				Real2 shift;
				auto checkTriangles = [&]() {
					auto triangleF = grid.findCellCrossedByTheRay(it, shift);
					auto triangleL = grid.locateOwnerCell(it, shift);
					
					if (triangleF.n == 3 && triangleL.n == 3) {
					// both are inner
						auto common = triangleF.equalPoints(triangleL);
						bool correct = common.size() == 3;
						if (common.size() == 2) {
							Real2 a = grid.coordsD(*common.begin());
							Real2 q = grid.coordsD(it) + shift;
							Real2 c = grid.coordsD(*common.rbegin());
							correct = linal::isDegenerate(a, q, c, EQUALITY_TOLERANCE * h) && 
									  linal::dotProduct(a - q, c - q) < 0;
						} else if (common.size() == 1) {
							correct = (shift == Real2::Zeros());
						}
						ASSERT_TRUE(correct)
							<< "it = " << grid.coordsD(it)
							<< "query = " << grid.coordsD(it) + shift
							<< "find: " << grid.coordsD(triangleF(0)) 
								<< grid.coordsD(triangleF(1)) << grid.coordsD(triangleF(2))
							<< "locate: " << grid.coordsD(triangleL(0))
								<< grid.coordsD(triangleL(1)) << grid.coordsD(triangleL(2))
							<< std::endl;
					
					} else {
					// both are outer
						ASSERT_TRUE((triangleF.n != 3) && (triangleL.n != 3))
							<< "F: " << triangleF.n << " L: " << triangleL.n << "\n"
							<< "it = " << grid.getIndex(it) << grid.coordsD(it) 
							<< "query = " << grid.coordsD(it) + shift
							<< "locate: " << grid.coordsD(triangleL(0))
							<< grid.coordsD(triangleL(1)) << grid.coordsD(triangleL(2))
							<< std::endl;;
					}
					
				};
				
				shift = {Utils::randomReal(-step, step), Utils::randomReal(-step, step)};
				checkTriangles();
				
//				shift = {Utils::randomReal(-step, step), 0};
//				checkTriangles();
				
//				shift = {0, Utils::randomReal(-step, step)};
//				checkTriangles();
			}
		}
	}
}


TEST(SimplexGrid2D, findCellCrossedByTheRayTwoBodies) {
	Task task;
	real h = 5;
	task.simplexGrid.spatialStep = h;
	
	task.simplexGrid.bodies = {
		Task::SimplexGrid::Body({ 0, {{3, 3}, {-3, 3}, {-3, -3}, {3, -3}}, { } }),
	};
	
	typedef SimplexGrid<2, CgalTriangulation> Grid;
	typedef typename Grid::Triangulation Triangulation;
	Triangulation triangulation(task);
	Grid oneBody(0, {&triangulation});
	
	
	task.simplexGrid.bodies.push_back(
			Task::SimplexGrid::Body({ 0, {{-2, 5}, {2, 5}, {0, 7}}, { } }));
	task.simplexGrid.bodies.push_back(
			Task::SimplexGrid::Body({ 0, {{-7, -1}, {-5, -4},  {-2, -5},  {2, -5}, {0, -7}, {-7, -7}}, { } }));
	task.simplexGrid.bodies.push_back(
			Task::SimplexGrid::Body({ 0, {{-10, -10}, {-10, 10}, {10, 10}, {10, -10}},
					{ {{-9, -9}, {-9, 9}, {9, 9}, {9, -9}} } }));
	
	Triangulation triangulation2(task);
	Grid twoBodies(0, {&triangulation2});
	
	// fortunately, triangulations of the first bodies are equal in
	// both triangulations with such parameters
	Utils::seedRand();
	for (int multiplier = 1; multiplier < 10; multiplier++) {
		real step = h / 3 * multiplier;
		for (auto& it1 : oneBody) {
			auto it2 = twoBodies.findVertexByCoordinates(oneBody.coordsD(it1));
			
			for (int test_counter = 0; test_counter < 8; test_counter++) {
				Real2 shift;
				auto checkTriangles = [&]() {
					auto triangle1 = oneBody.findCellCrossedByTheRay(it1, shift);
					auto triangle2 = twoBodies.findCellCrossedByTheRay(it2, shift);
					
					if (triangle1.n == 3 && triangle2.n == 3) {
					// both are inner
						elements::Triangle<Real2> t1(triangle1, 
								std::function<Real2(const Grid::Iterator&)>(
										[&](const Grid::Iterator& iterator) {
											return oneBody.coordsD(iterator); 
										}));
						elements::Triangle<Real2> t2(triangle2, 
								std::function<Real2(const Grid::Iterator&)>(
										[&](const Grid::Iterator& iterator) {
											return twoBodies.coordsD(iterator); 
										}));

						auto common = t1.equalPoints(t2);
						bool correct = common.size() == 3;
						if (common.size() == 2) {
							Real2 a = *common.begin();
							Real2 q = oneBody.coordsD(it1) + shift;
							Real2 c = *common.rbegin();
							correct = linal::isDegenerate(a, q, c, EQUALITY_TOLERANCE * h) && 
									  linal::dotProduct(a - q, c - q) < 0;
						}
						ASSERT_TRUE(correct)
							<< "it = " << oneBody.coordsD(it1)
							<< "query = " << oneBody.coordsD(it1) + shift
							<< "one: " << t1 << "two: " << t2
							<< std::endl;
					
					} else {
					// both are outer
						ASSERT_TRUE((triangle1.n != 3) && (triangle2.n != 3))
							<< "one: " << triangle1.n << " two: " << triangle2.n << "\n"
							<< "it1 = " << oneBody.getIndex(it1) << oneBody.coordsD(it1) 
							<< "query = " << oneBody.coordsD(it1) + shift
							<< std::endl;
					}
					
				};
				
				shift = {Utils::randomReal(-step, step), Utils::randomReal(-step, step)};
				checkTriangles();
				
//				shift = {Utils::randomReal(-step, step), 0};
//				checkTriangles();
				
//				shift = {0, Utils::randomReal(-step, step)};
//				checkTriangles();
			}
		}
	}
}


TEST(SimplexGrid2D, miscellaneous) {
	Task task;
	task.simplexGrid.spatialStep = 0.5;
	task.simplexGrid.bodies = {
		Task::SimplexGrid::Body({ 0, { {0, 0}, {0, 1}, {1, 1}, {1, 0} }, { } })
	};
	
	typedef SimplexGrid<2, CgalTriangulation> Grid;
	typedef typename Grid::Triangulation Triangulation;
	Triangulation triangulation(task);
	Grid grid(0, {&triangulation});
	
	ASSERT_EQ(21, grid.sizeOfAllNodes());
	ASSERT_EQ(grid.sizeOfRealNodes(), grid.sizeOfAllNodes());
	ASSERT_NEAR(0.2210, grid.getAverageHeight(), 1e-4);
	
	
	ASSERT_EQ(Real2::Zeros(), grid.borderNormal(*(grid.innerBegin())));
	
	for (auto borderIt  = grid.borderBegin();
	          borderIt != grid.borderEnd(); ++borderIt) {
		
		Real2 coords = grid.coordsD(*borderIt);
		Real2 normal = Real2::Zeros();
		for (int i = 0; i < 2; i++) {
			if      (coords(i) == 0) { normal(i) = -1; }
			else if (coords(i) == 1) { normal(i) =  1; }
		}
		normal = linal::normalize(normal);
		
		Real2 error = grid.borderNormal(*borderIt) - normal;
		ASSERT_LT(linal::length(error), EQUALITY_TOLERANCE)
				<< "coords:" << coords 
				<< "correct normal:" << normal 
				<< "actual normal:" << grid.borderNormal(*borderIt);
	}
	
	
	ASSERT_EQ(2, grid.findNeighborVertices(grid.findVertexByCoordinates({0, 0})).size());
	ASSERT_EQ(4, grid.findNeighborVertices(grid.findVertexByCoordinates({0.5, 0.5})).size());
}





