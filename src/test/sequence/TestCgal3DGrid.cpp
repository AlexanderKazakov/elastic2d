#include <lib/mesh/grid/Cgal3DGrid.hpp>
#include <lib/util/snapshot/VtkSnapshotter.hpp>

#include <gtest/gtest.h>

using namespace gcm;



TEST(Cgal3DGrid, miscellaneous) {
	Task task;

	try {
		task.cgal3DGrid.spatialStep = 0.2;
		task.cgal3DGrid.detectSharpEdges = true;
		task.cgal3DGrid.polyhedronFileName = "meshes/cube.off";
//		task.cgal3DGrid.spatialStep = 0.03;
//		task.cgal3DGrid.polyhedronFileName = "meshes/elephant.off";
	
		Cgal3DGrid grid(task);
//		ASSERT_EQ(21, grid.sizeOfAllNodes());
		ASSERT_EQ(grid.sizeOfRealNodes(), grid.sizeOfAllNodes());
		
//		ASSERT_THROW(grid.normal(*(grid.innerBegin())), Exception);
		
		ASSERT_EQ(task.cgal3DGrid.spatialStep, grid.getMinimalSpatialStep());
			
		VtkUtils::dumpGridToVtk(grid);
	
	} catch (Exception e) {
		std::cerr << e.what() << std::endl;
	}
}





