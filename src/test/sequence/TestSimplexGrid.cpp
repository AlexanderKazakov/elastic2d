#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>
#include <libgcm/grid/simplex/SimplexGrid.hpp>
#include <libgcm/util/snapshot/VtkSnapshotter.hpp>

#include <gtest/gtest.h>

using namespace gcm;


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
	ASSERT_NEAR(0.22, grid.getAverageHeight(), 0.01);
	
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


TEST(SimplexGrid3D, miscellaneous) {
	Task task;
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::CGAL_MESHER;
	task.simplexGrid.spatialStep = 0.2;
	task.simplexGrid.detectSharpEdges = true;
	task.simplexGrid.fileName = "meshes/cube.off";
	
//	try {
	typedef SimplexGrid<3, CgalTriangulation> Grid;
	typedef typename Grid::Triangulation Triangulation;
	Triangulation triangulation(task);
	Grid grid(0, {&triangulation});
	
	ASSERT_NEAR(190, (real)grid.sizeOfAllNodes(), 10);
	ASSERT_EQ(grid.sizeOfRealNodes(), grid.sizeOfAllNodes());
	ASSERT_NEAR(0.1384, grid.getAverageHeight(), 0.03);
	
	ASSERT_EQ(Real3::Zeros(), grid.borderNormal(*(grid.innerBegin())));
	
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


