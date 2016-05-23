#include <lib/mesh/grid/Cgal3DGrid.hpp>
#include <lib/util/snapshot/VtkSnapshotter.hpp>

#include <gtest/gtest.h>

using namespace gcm;



TEST(Cgal3DGrid, miscellaneous) {
	Task task;

//	task.cgal3DGrid.spatialStep = 1.0;
//	task.cgal3DGrid.polyhedronFileName = "meshes/tetrahedron.off";
	task.cgal3DGrid.spatialStep = 0.05;
	task.cgal3DGrid.polyhedronFileName = "meshes/elephant.off";

	Cgal3DGrid grid(task);
//	ASSERT_EQ(21, grid.sizeOfAllNodes());
	ASSERT_EQ(grid.sizeOfRealNodes(), grid.sizeOfAllNodes());
	
	ASSERT_THROW(grid.normal(*(grid.innerBegin())), Exception);
	
	ASSERT_EQ(task.cgal3DGrid.spatialStep, grid.getMinimalSpatialStep());
	
//	ASSERT_EQ(2, grid.findNeighborVertices(grid.findVertexByCoordinates({0, 0})).size());
//	ASSERT_EQ(4, grid.findNeighborVertices(grid.findVertexByCoordinates({0.5, 0.5})).size());
	
	VtkUtils::dumpGridToVtk(grid);
}





