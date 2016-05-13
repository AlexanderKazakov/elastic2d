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
		Task::Cgal2DGrid::Body({{-2, 5}, {2, 5}, {0, 7}}, {})
	};
	
	Cgal2DGrid grid(task);
	
	Utils::seedRand();
	for (int multiplier = 1; multiplier < 10; multiplier++) {
		int cntXY = 0, cntX = 0, cntY = 0;
		real step = h / 3 * multiplier;
		for (auto& it : grid) {
			for (int test_counter = 0; test_counter < 8; test_counter++) {
				Real2 shift;
				auto checkInnerTriangle = [&](int& counter) {
					auto triangleF = grid.findOwnerTriangle(it, shift);
					auto triangleL = grid.locateOwnerTriangle(it, shift);
					
					if (triangleF.valid) {
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
					
					if (triangleL.valid) {
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
		
		ASSERT_TRUE(Utils::approximatelyEqual(2400, cntXY, 0.3)) << "Actual: " << cntXY;
		ASSERT_TRUE(Utils::approximatelyEqual(2400, cntX, 0.3)) << "Actual: " << cntX;
		ASSERT_TRUE(Utils::approximatelyEqual(2400, cntY, 0.3)) << "Actual: " << cntY;
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
	for (int multiplier = 1; multiplier < 10; multiplier++) {
		real step = h / 3 * multiplier;
		for (auto& it : grid) {
			for (int test_counter = 0; test_counter < 8; test_counter++) {
			
				Real2 shift;
				auto checkTriangles = [&]() {
					auto triangleF = grid.findOwnerTriangle(it, shift);
					auto triangleL = grid.locateOwnerTriangle(it, shift);
					
					if (triangleF.valid && triangleL.valid) {
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
						ASSERT_TRUE(!triangleF.valid && !triangleL.valid) 
							<< "F: " << triangleF.valid << " L: " << triangleL.valid << "\n"
							<< "it = " << grid.getIndex(it) << grid.coords2d(it) 
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


TEST(Cgal2DGrid, findOwnerTriangleTwoBodies) {
	Task task;
	real h = 5;
	task.cgal2DGrid.spatialStep = h;
	
	task.cgal2DGrid.bodies = {
		Task::Cgal2DGrid::Body({{3, 3}, {-3, 3}, {-3, -3}, {3, -3}}, { }),
	};
	Cgal2DGrid oneBody(task);
	VtkUtils::dumpGridToVtk(oneBody, "one");
	
	task.cgal2DGrid.bodies.push_back(
			Task::Cgal2DGrid::Body({{-2, 5}, {2, 5}, {0, 7}}, { }));
	Cgal2DGrid twoBodies(task);
	VtkUtils::dumpGridToVtk(twoBodies, "two");

	
	// fortunately, triangulations of the first bodies are equal in
	// both triangulations with such parameters
	Utils::seedRand();
	/*for (int multiplier = 1; multiplier < 10; multiplier++)*/ {
//		real step = h / 3 * multiplier;
		/*for (auto& it1 : oneBody)*/ {
			Cgal2DGrid::Iterator it1 = 7;
			auto it2 = twoBodies.findVertexByCoordinates(oneBody.coords2d(it1));
			
			/*for (int test_counter = 0; test_counter < 8; test_counter++)*/ {
				Real2 shift;
				auto checkTriangles = [&]() {
					auto triangle1 = oneBody.findOwnerTriangle(it1, shift);
					auto triangle2 = twoBodies.findOwnerTriangle(it2, shift);
					
					if (triangle1.valid && triangle2.valid) {
					// both are inner
						elements::Triangle<Real2> t1(triangle1, 
								std::function<Real2(const Cgal2DGrid::Iterator&)>(
										[&](const Cgal2DGrid::Iterator& iterator) {
											return oneBody.coords2d(iterator); 
										}));
						elements::Triangle<Real2> t2(triangle2, 
								std::function<Real2(const Cgal2DGrid::Iterator&)>(
										[&](const Cgal2DGrid::Iterator& iterator) {
											return twoBodies.coords2d(iterator); 
										}));

						auto common = t1.equalPoints(t2);
						bool correct = common.size() == 3;
						if (common.size() == 2) {
							Real2 a = *common.begin();
							Real2 q = oneBody.coords2d(it1) + shift;
							Real2 c = *common.rbegin();
							correct = linal::isDegenerate(a, q, c) && 
									  linal::dotProduct(a - q, c - q) < 0;
						}
						ASSERT_TRUE(correct)
							<< "it = " << oneBody.coords2d(it1)
							<< "query = " << oneBody.coords2d(it1) + shift
							<< "one: " << t1 << "two: " << t2
							<< std::endl;
					
					} else {
					// both are outer
						ASSERT_TRUE(!triangle1.valid && !triangle2.valid)
							<< "one: " << triangle1.valid << " two: " << triangle2.valid << "\n"
							<< "it1 = " << oneBody.getIndex(it1) << oneBody.coords2d(it1) 
							<< "query = " << oneBody.coords2d(it1) + shift
							<< std::endl;
					}
					
				};
				
//				shift = {Utils::randomReal(-step, step), Utils::randomReal(-step, step)};
//				checkTriangles();
				
//				shift = {Utils::randomReal(-step, step), 0};
//				checkTriangles();
				
//				shift = {0, Utils::randomReal(-step, step)};
				shift = {2, 0};
				checkTriangles();
			}
		}
	}
}
