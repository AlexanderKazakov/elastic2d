#include <lib/mesh/grid/Cgal3DGrid.hpp>
#include <lib/util/snapshot/VtkSnapshotter.hpp>

#include <gtest/gtest.h>

using namespace gcm;



TEST(Cgal3DGrid, miscellaneous) {
	Task task;
	task.cgal3DGrid.spatialStep = 0.2;
	task.cgal3DGrid.detectSharpEdges = true;
	task.cgal3DGrid.polyhedronFileName = "meshes/cube.off";
	Cgal3DGrid grid(task);
	ASSERT_EQ(190, grid.sizeOfAllNodes());
	ASSERT_EQ(grid.sizeOfRealNodes(), grid.sizeOfAllNodes());
	ASSERT_EQ(task.cgal3DGrid.spatialStep, grid.getMinimalSpatialStep());
	
	ASSERT_THROW(grid.normal(*(grid.innerBegin())), Exception);
	
	for (auto borderIt  = grid.borderBegin();
	          borderIt != grid.borderEnd(); ++borderIt) {
		
		Real3 coords = grid.coords(*borderIt);
		Real3 normal = Real3::Zeros();
		for (int i = 0; i < 3; i++) {
			if      (coords(i) == 0) { normal(i) = -1; }
			else if (coords(i) == 1) { normal(i) =  1; }
		}
		normal = linal::normalize(normal);
		
		Real3 error = grid.normal(*borderIt) - normal;
		ASSERT_LT(linal::length(error), 0.3)
				<< "coords:" << coords 
				<< "correct normal:" << normal 
				<< "actual normal:" << grid.normal(*borderIt);
	}
	
}


inline void checkOwnerCellVsBarycentric(const Cgal3DGrid* grid, const real step,
                                        int& cntFind, int& cntLocate) {
	cntFind = cntLocate = 0;
	for (const auto& it : *grid) {
		for (int test_counter = 0; test_counter < 30; test_counter++) {
			Real3 shift = { Utils::randomReal(-step, step), 
			                Utils::randomReal(-step, step),
			                Utils::randomReal(-step, step) };
			
			Cgal3DGrid::Cell tetrF;
			Cgal3DGrid::Cell tetrL;
			tetrF = grid->findOwnerCell(it, shift);
			tetrL = grid->locateOwnerCell(it, shift);
			
			if (tetrF.valid) {
				Real4 lambda = linal::barycentricCoordinates(
						grid->coords(tetrF(0)), 
						grid->coords(tetrF(1)), 
						grid->coords(tetrF(2)), 
						grid->coords(tetrF(3)), 
						grid->coords(it) + shift);
				for (int i = 0; i < 4; i++) {
					ASSERT_TRUE(lambda(i) > -EQUALITY_TOLERANCE) << lambda(i);
				}
				cntFind++;
			}
			
			if (tetrL.valid) {
				Real4 lambda = linal::barycentricCoordinates(
						grid->coords(tetrL(0)), 
						grid->coords(tetrL(1)), 
						grid->coords(tetrL(2)), 
						grid->coords(tetrL(3)), 
						grid->coords(it) + shift);
				for (int i = 0; i < 4; i++) {
					ASSERT_TRUE(lambda(i) > -EQUALITY_TOLERANCE) << lambda(i);
				}
				cntLocate++;
			}
			
		}
	}
}


TEST(Cgal3DGrid, ownerTetrahedronVsBarycentric) {
	Task task;
	real h = 0.1;
	task.cgal3DGrid.spatialStep = h;
	task.cgal3DGrid.detectSharpEdges = true;
	
	task.cgal3DGrid.polyhedronFileName = "meshes/tetrahedron.off";
	Cgal3DGrid tetrGrid(task);
	Utils::seedRand();
	for (int multiplier = 1; multiplier < 30; multiplier++) {
//		std::cout << "tetr: multiplier == " << multiplier << std::endl;
		int cntFind = 0, cntLocate = 0;
		checkOwnerCellVsBarycentric(&tetrGrid, h / 3 * multiplier, cntFind, cntLocate);
		ASSERT_EQ(cntLocate, cntFind)
				<< "Find: " << cntFind << " Locate: " << cntLocate;
	}
	
	task.cgal3DGrid.polyhedronFileName = "meshes/cube.off";
	Cgal3DGrid cubeGrid(task);
	for (int multiplier = 1; multiplier < 30; multiplier++) {
//		std::cout << "cube: multiplier == " << multiplier << std::endl;
		int cntFind = 0, cntLocate = 0;
		checkOwnerCellVsBarycentric(&cubeGrid, h / 3 * multiplier, cntFind, cntLocate);
		ASSERT_EQ(cntLocate, cntFind)
				<< "Find: " << cntFind << " Locate: " << cntLocate;
	}
}



