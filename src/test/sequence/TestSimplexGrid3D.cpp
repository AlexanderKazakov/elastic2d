#include <lib/mesh/grid/cgal/CgalTriangulation.hpp>
#include <lib/mesh/grid/SimplexGrid.hpp>
#include <lib/util/snapshot/VtkSnapshotter.hpp>

#include <gtest/gtest.h>

using namespace gcm;



TEST(SimplexGrid3D, miscellaneous) {
	Task task;
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::CGAL_MESHER;
	task.simplexGrid.spatialStep = 0.2;
	task.simplexGrid.detectSharpEdges = true;
	task.simplexGrid.fileName = "meshes/cube.off";
	
//	try {
	typedef SimplexGrid<3, CgalTriangulation> Grid;
	typedef typename Grid::GlobalScene GS;
	std::shared_ptr<GS> gs(new GS(task));
	Grid grid(task, gs.get(), 0);
	
	
	ASSERT_EQ(190, grid.sizeOfAllNodes());
	ASSERT_EQ(grid.sizeOfRealNodes(), grid.sizeOfAllNodes());
	ASSERT_NEAR(0.1384, grid.getMinimalSpatialStep(), 1e-4);
	
	ASSERT_THROW(grid.borderNormal(*(grid.innerBegin())), Exception);
	
	for (auto borderIt  = grid.borderBegin();
	          borderIt != grid.borderEnd(); ++borderIt) {
		
		Real3 coords = grid.coords(*borderIt);
		Real3 normal = Real3::Zeros();
		for (int i = 0; i < 3; i++) {
			if      (coords(i) == 0) { normal(i) = -1; }
			else if (coords(i) == 1) { normal(i) =  1; }
		}
		normal = linal::normalize(normal);
		
		Real3 error = grid.borderNormal(*borderIt) - normal;
		ASSERT_LT(linal::length(error), 0.3)
				<< "coords:" << coords 
				<< "correct normal:" << normal 
				<< "actual normal:" << grid.borderNormal(*borderIt);
	}
	
//	} catch (Exception e) {
//		std::cout << e.what();
//	}
}


inline void checkOwnerCellVsBarycentric(
		const SimplexGrid<3, CgalTriangulation>* grid, const real step,
		int& cntFind, int& cntLocate) {
	
	cntFind = cntLocate = 0;
	for (const auto& it : *grid) {
		for (int test_counter = 0; test_counter < 30; test_counter++) {
			Real3 shift = { Utils::randomReal(-step, step), 
			                Utils::randomReal(-step, step),
			                Utils::randomReal(-step, step) };
			
			typename SimplexGrid<3, CgalTriangulation>::Cell tetrF;
			typename SimplexGrid<3, CgalTriangulation>::Cell tetrL;
			tetrF = grid->findOwnerCell(it, shift);
			tetrL = grid->locateOwnerCell(it, shift);
			
			if (tetrF.n == 4) {
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
			
			} else if (tetrF.n == 3) {
				Real3 cross = linal::lineWithFlatIntersection(
						grid->coords(tetrF(0)), 
						grid->coords(tetrF(1)), 
						grid->coords(tetrF(2)),
						grid->coords(it),
						grid->coords(it) + shift);
				Real3 lambda = linal::barycentricCoordinates(
						grid->coords(tetrF(0)), 
						grid->coords(tetrF(1)), 
						grid->coords(tetrF(2)), 
						cross);
				for (int i = 0; i < 3; i++) {
					ASSERT_TRUE(lambda(i) > -EQUALITY_TOLERANCE) << lambda(i);
				}
			}
			
			if (tetrL.n == 4) {
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


TEST(SimplexGrid3D, ownerTetrahedronVsBarycentric) {
	Task task;
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::CGAL_MESHER;
	task.simplexGrid.spatialStep = 0.4;
	task.simplexGrid.detectSharpEdges = true;
	
	task.simplexGrid.fileName = "meshes/tetrahedron.off";
	
	typedef SimplexGrid<3, CgalTriangulation> Grid;
	typedef typename Grid::GlobalScene GS;
	std::shared_ptr<GS> gs(new GS(task));
	Grid tetrGrid(task, gs.get(), 0);
	
	Utils::seedRand();
	for (int multiplier = 1; multiplier < 30; multiplier++) {
//		std::cout << "tetr: multiplier == " << multiplier << std::endl;
		int cntFind = 0, cntLocate = 0;
		checkOwnerCellVsBarycentric(&tetrGrid, 
				task.simplexGrid.spatialStep / 3 * multiplier, cntFind, cntLocate);
		ASSERT_EQ(cntLocate, cntFind)
				<< "Find: " << cntFind << " Locate: " << cntLocate;
	}
	
	
	
	task.simplexGrid.spatialStep = 0.2;
	task.simplexGrid.fileName = "meshes/cube.off";
	
	std::shared_ptr<GS> gs2(new GS(task));
	Grid cubeGrid(task, gs2.get(), 0);
	
	for (int multiplier = 1; multiplier < 30; multiplier++) {
//		std::cout << "cube: multiplier == " << multiplier << std::endl;
		int cntFind = 0, cntLocate = 0;
		checkOwnerCellVsBarycentric(&cubeGrid,
				task.simplexGrid.spatialStep / 3 * multiplier, cntFind, cntLocate);
		ASSERT_TRUE(abs(cntLocate - cntFind) < cntLocate / 20)
				<< "Find: " << cntFind << " Locate: " << cntLocate;
	}
	
	
	task.simplexGrid.spatialStep = 0.4;
	task.simplexGrid.fileName = "meshes/icosahedron.off";
	task.simplexGrid.detectSharpEdges = false;
	
	std::shared_ptr<GS> gs3(new GS(task));
	Grid icosGrid(task, gs3.get(), 0);
	
	for (int multiplier = 1; multiplier < 30; multiplier++) {
//		std::cout << "icosahedron: multiplier == " << multiplier << std::endl;
		int cntFind = 0, cntLocate = 0;
		checkOwnerCellVsBarycentric(&icosGrid,
				task.simplexGrid.spatialStep / 3 * multiplier, cntFind, cntLocate);
		ASSERT_GE(cntLocate, cntFind) // non-convex case - find less or equal than locate
				<< "multiplier == " << multiplier 
				<< " Find: " << cntFind << " Locate: " << cntLocate << std::endl;
	}
}



