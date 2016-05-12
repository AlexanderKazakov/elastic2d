#include <gtest/gtest.h>

#include <lib/mesh/grid/Cgal2DGrid.hpp>

#include <lib/util/snapshot/VtkSnapshotter.hpp>

using namespace gcm;


TEST(Cgal2DGrid, ownerTriangleVsBarycentric) {
	Task task;
	real h = 0.5;
	task.cgal2DGrid.spatialStep = h;
	
	Task::Cgal2DGrid::Body::Border outer = {
		{3, 3}, {-3, 3}, {-3, -3}, {3, -3}, {2, 2}
	};
	Task::Cgal2DGrid::Body::Border inner = {
		{-2, -1}, {-1, 0}, {0, -1}, {-1, -2}
	};
	task.cgal2DGrid.bodies = {
		Task::Cgal2DGrid::Body(outer, {inner}),
	};
	
	Cgal2DGrid grid(task);
	
	Utils::seedRand();
	int cntXY = 0, cntX = 0, cntY = 0;
	for (int crnt = 1; crnt < 10; crnt++) {
		real step = h / 3 * crnt;
		for (auto& it : grid) {
			for (int test_counter = 0; test_counter < 8; test_counter++) {
				Real2 shift;
				auto checkInnerTriangle = [&](int& counter) {
					auto triangleF = grid.findOwnerTriangle(it, shift);
					auto triangleL = grid.locateOwnerTriangle(it, shift);
					
					if (triangleF.inner) {
						auto lambda = linal::barycentricCoordinates(
								grid.coords2d(triangleF.p[0]), 
								grid.coords2d(triangleF.p[1]), 
								grid.coords2d(triangleF.p[2]), 
								grid.coords2d(it) + shift);
						for (int i = 0; i < 3; i++) {
							ASSERT_TRUE(lambda(i) > -EQUALITY_TOLERANCE) << lambda(i);
						}
						counter++;
					}
					
					if (triangleL.inner) {
						auto lambda = linal::barycentricCoordinates(
								grid.coords2d(triangleL.p[0]), 
								grid.coords2d(triangleL.p[1]), 
								grid.coords2d(triangleL.p[2]), 
								grid.coords2d(it) + shift);
						for (int i = 0; i < 3; i++) {
							ASSERT_TRUE(lambda(i) > -EQUALITY_TOLERANCE) << lambda(i);
						}
					}
				};
				
				shift = {Utils::randomReal(-step, step), Utils::randomReal(-step, step)};
				checkInnerTriangle(cntXY);
				
				shift = {Utils::randomReal(-step, step), 0};
				checkInnerTriangle(cntX);
				
				shift = {0, Utils::randomReal(-step, step)};
				checkInnerTriangle(cntY);
			}
		}
		
		ASSERT_TRUE(Utils::approximatelyEqual(2400, cntXY, 100)) << "Actual: " << cntXY;
		ASSERT_TRUE(Utils::approximatelyEqual(2400, cntX, 100)) << "Actual: " << cntX;
		ASSERT_TRUE(Utils::approximatelyEqual(2400, cntY, 100)) << "Actual: " << cntY;
	}
}


TEST(Cgal2DGrid, locateVsFindOwnerTriangle) {
	Task task;
	real h = 0.5;
	task.cgal2DGrid.spatialStep = h;
	
	Task::Cgal2DGrid::Body::Border outer = {
		{3, 3}, {-3, 3}, {-3, -3}, {3, -3}, {5, 1}
	};
	task.cgal2DGrid.bodies = {
		Task::Cgal2DGrid::Body(outer, { })
	};
	
	Cgal2DGrid grid(task);
	
	Utils::seedRand();
	for (int crnt = 1; crnt < 10; crnt++) {
		real step = h / 3 * crnt;
		for (auto& it : grid) {
			for (int test_counter = 0; test_counter < 8; test_counter++) {
			
				Real2 shift;
				auto checkTriangles = [&]() {
					auto triangleF = grid.findOwnerTriangle(it, shift);
					auto triangleL = grid.locateOwnerTriangle(it, shift);
					
					if (triangleF.inner && triangleL.inner) {
					// both are inner
						auto common = triangleF.equalPoints(triangleL);
						bool correct = common.size() == 3;
						if (common.size() == 2) {
							Real2 a = grid.coords2d(*common.begin());
							Real2 q = grid.coords2d(it) + shift;
							Real2 c = grid.coords2d(*common.rbegin());
							correct = linal::isDegenerate(a, q, c) && 
									  linal::dotProduct(a - q, c - q) < 0;
						}
						ASSERT_TRUE(correct)
							<< "it = " << grid.coords2d(it)
							<< "query = " << grid.coords2d(it) + shift
							<< "find: " << grid.coords2d(triangleF.p[0]) 
								<< grid.coords2d(triangleF.p[1]) << grid.coords2d(triangleF.p[2])
							<< "locate: " << grid.coords2d(triangleL.p[0])
								<< grid.coords2d(triangleL.p[1]) << grid.coords2d(triangleL.p[2])
							<< std::endl;
					
					} else {
					// both are outer
						ASSERT_TRUE(!triangleF.inner && !triangleL.inner) 
							<< "F: " << triangleF.inner << " L: " << triangleL.inner << "\n"
							<< "it = " << it.iter << grid.coords2d(it) 
							<< "query = " << grid.coords2d(it) + shift
							<< "locate: " << grid.coords2d(triangleL.p[0])
							<< grid.coords2d(triangleL.p[1]) << grid.coords2d(triangleL.p[2])
							<< std::endl;;
					}
					
				};
				
				shift = {Utils::randomReal(-step, step), Utils::randomReal(-step, step)};
				checkTriangles();
				
				shift = {Utils::randomReal(-step, step), 0};
				checkTriangles();
				
				shift = {0, Utils::randomReal(-step, step)};
				checkTriangles();
			}
		}
	}
}

