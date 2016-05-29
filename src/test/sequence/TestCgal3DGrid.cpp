#include <lib/mesh/grid/Cgal3DGrid.hpp>
#include <lib/util/snapshot/VtkSnapshotter.hpp>

#include <gtest/gtest.h>

using namespace gcm;



//TEST(Cgal3DGrid, miscellaneous) {
//	Task task;
//	task.cgal3DGrid.spatialStep = 0.2;
//	task.cgal3DGrid.detectSharpEdges = true;
//	task.cgal3DGrid.polyhedronFileName = "meshes/cube.off";
//	Cgal3DGrid grid(task);
//	ASSERT_EQ(190, grid.sizeOfAllNodes());
//	ASSERT_EQ(grid.sizeOfRealNodes(), grid.sizeOfAllNodes());
//	ASSERT_EQ(task.cgal3DGrid.spatialStep, grid.getMinimalSpatialStep());
	
//	ASSERT_THROW(grid.normal(*(grid.innerBegin())), Exception);
	
//	for (auto borderIt  = grid.borderBegin();
//	          borderIt != grid.borderEnd(); ++borderIt) {
		
//		Real3 coords = grid.coords(*borderIt);
//		Real3 normal = Real3::Zeros();
//		for (int i = 0; i < 3; i++) {
//			if      (coords(i) == 0) { normal(i) = -1; }
//			else if (coords(i) == 1) { normal(i) =  1; }
//		}
//		normal = linal::normalize(normal);
		
//		Real3 error = grid.normal(*borderIt) - normal;
//		ASSERT_LT(linal::length(error), 0.3)
//				<< "coords:" << coords 
//				<< "correct normal:" << normal 
//				<< "actual normal:" << grid.normal(*borderIt);
//	}
	
//}


//TEST(Cgal3DGrid, ownerTetrahedronVsBarycentric) {
//	Task task;
//	real h = 0.1;
//	task.cgal3DGrid.spatialStep = h;
//	task.cgal3DGrid.detectSharpEdges = true;
//	task.cgal3DGrid.polyhedronFileName = "meshes/cube.off";
//	Cgal3DGrid grid(task);
	
//	Utils::seedRand();
//	int cntXYZ = 0, cntX = 0, cntY = 0, cntZ = 0;
//	for (int multiplier = 0; multiplier < 10; multiplier++) {
//		real step = h / 3 * multiplier;
//		for (auto& it : grid) {
//			for (int test_counter = 0; test_counter < 8; test_counter++) {
//				Real3 shift;
				
//				auto checkInnerTriangle = [&](int& counter) {
//					Cgal3DGrid::Cell tetrF;
//					Cgal3DGrid::Cell tetrL;
//					try {
//						tetrF = grid.findOwnerCell(it, shift);
//						tetrL = grid.locateOwnerCell(it, shift);
//					} catch (Exception e) {
//						std::cerr << e.what() << std::endl;
//					}
					
//					if (tetrF.valid) {
//						Real4 lambda = linal::barycentricCoordinates(
//								grid.coords(tetrF(0)), 
//								grid.coords(tetrF(1)), 
//								grid.coords(tetrF(2)), 
//								grid.coords(tetrF(3)), 
//								grid.coords(it) + shift);
//						for (int i = 0; i < 4; i++) {
//							ASSERT_TRUE(lambda(i) > -EQUALITY_TOLERANCE) << lambda(i);
//						}
//						counter++;
//					}
					
//					if (tetrL.valid) {
//						Real4 lambda = linal::barycentricCoordinates(
//								grid.coords(tetrL(0)), 
//								grid.coords(tetrL(1)), 
//								grid.coords(tetrL(2)), 
//								grid.coords(tetrL(3)), 
//								grid.coords(it) + shift);
//						for (int i = 0; i < 4; i++) {
//							ASSERT_TRUE(lambda(i) > -EQUALITY_TOLERANCE) << lambda(i);
//						}
//					}
//				};
				
//				shift = { Utils::randomReal(-step, step), 
//				          Utils::randomReal(-step, step),
//				          Utils::randomReal(-step, step)};
//				checkInnerTriangle(cntXYZ);
				
//				shift = {Utils::randomReal(-step, step), 0, 0};
//				checkInnerTriangle(cntX);
				
//				shift = {0, Utils::randomReal(-step, step), 0};
//				checkInnerTriangle(cntY);
				
//				shift = {0, 0, Utils::randomReal(-step, step)};
//				checkInnerTriangle(cntZ);
//			}
//		}
//	}
	
//	ASSERT_TRUE(Utils::approximatelyEqual(70000, cntX, 0.1)) << "Actual: " << cntX;
//	ASSERT_TRUE(Utils::approximatelyEqual(70000, cntY, 0.1)) << "Actual: " << cntY;
//	ASSERT_TRUE(Utils::approximatelyEqual(70000, cntZ, 0.1)) << "Actual: " << cntZ;	
//	ASSERT_TRUE(Utils::approximatelyEqual(50000, cntXYZ, 0.1)) << "Actual: " << cntXYZ;
//}


TEST(Cgal3DGrid, locateVsFindOwnerCell) {
	Task task;
	real h = 0.2;
	task.cgal3DGrid.spatialStep = h;
	task.cgal3DGrid.detectSharpEdges = true;
	task.cgal3DGrid.polyhedronFileName = "meshes/tetrahedron.off";
	Cgal3DGrid grid(task);
	VtkUtils::dumpGridToVtk(grid);
	
	std::function<Real3(const Cgal3DGrid::Iterator&)> coordinates = 
			[&](const Cgal3DGrid::Iterator& iterator) {
				return grid.coords(iterator); 
			};
	
	Utils::seedRand();
	for (int multiplier = 0; multiplier < 10; multiplier++) {
		std::cout << "multiplier = " << multiplier << std::endl;
		real step = h / 3 * multiplier;
		for (auto& it : grid) {
//			Cgal3DGrid::Iterator it = 242;///////
			for (int test_counter = 0; test_counter < 8; test_counter++) {
				Real3 shift;
				
				auto checkTetrs = [&]() {
					auto tetrL = grid.locateOwnerCell(it, shift);
					auto tetrF = grid.findOwnerCell(it, shift);
					
					if (tetrF.valid && tetrL.valid) {
					// both are inner
						auto common = tetrF.equalPoints(tetrL);
						bool correct = common.size() == 4;
						if (common.size() == 3) {
							Real3 a = grid.coords(*(common.begin()));
							Real3 b = grid.coords(*(std::next(common.begin())));
							Real3 c = grid.coords(*(common.rbegin()));
							Real3 q = grid.coords(it) + shift;
							if (fabs(linal::orientedVolume(a, b, c, q)) < EQUALITY_TOLERANCE) {
								Real3 lambda = linal::barycentricCoordinates(a, b, c, q);
								correct = lambda(0) > -EQUALITY_TOLERANCE && 
								          lambda(1) > -EQUALITY_TOLERANCE &&
								          lambda(2) > -EQUALITY_TOLERANCE;
							}
						} else if (common.size() == 2) {
							Real3 a = grid.coords(*(common.begin()));
							Real3 b = grid.coords(*(common.rbegin()));
							Real3 q = grid.coords(it) + shift;
							
							correct = linal::length(linal::crossProduct(
									a - b, q - b)) < EQUALITY_TOLERANCE && 
									linal::dotProduct(a - q, b - q) < 0;
						} else if (common.size() == 1) {
							correct = (shift == Real3::Zeros());
						}
						ASSERT_TRUE(correct)
							<< "common.size() == " << common.size()
							<< "\nit = " << grid.getIndex(it) << grid.coords(it)
							<< "query = " << grid.coords(it) + shift
							<< "find:" << elements::Tetrahedron<Real3>(tetrF, coordinates)
							<< "locate:" << elements::Tetrahedron<Real3>(tetrL, coordinates);
					
					} else {
					// both are outer
						ASSERT_TRUE(!tetrF.valid && !tetrL.valid) 
							<< "it = " << grid.getIndex(it) << grid.coords(it)
							<< "query = " << grid.coords(it) + shift
							<< "find:" << elements::Tetrahedron<Real3>(tetrF, coordinates)
							<< "locate:" << elements::Tetrahedron<Real3>(tetrL, coordinates);
					}
					
				};
				
				shift = { Utils::randomReal(-step, step), 
				          Utils::randomReal(-step, step),
				          Utils::randomReal(-step, step)};
				checkTetrs();
				
				shift = {Utils::randomReal(-step, step), 0, 0};
				checkTetrs();
				
				shift = {0, Utils::randomReal(-step, step), 0};
				checkTetrs();
				
				shift = {0, 0, Utils::randomReal(-step, step)};
//				shift = {0.05, 0, 0};///////
				checkTetrs();
			}
		}
	}
}





