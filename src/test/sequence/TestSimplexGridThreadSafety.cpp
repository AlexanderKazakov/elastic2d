#include <lib/mesh/grid/cgal/CgalTriangulation.hpp>
#include <lib/mesh/grid/SimplexGrid.hpp>
#include <lib/mesh/grid/SimplexGlobalScene.hpp>

#include <gtest/gtest.h>

using namespace gcm;


TEST(SimplexGrid3D, ThreadSafety) {
	Task task;
	task.simplexGrid.mesher = Task::SimplexGrid::Mesher::CGAL_MESHER;
	task.simplexGrid.spatialStep = 0.1;
	task.simplexGrid.detectSharpEdges = true;
	task.simplexGrid.fileName = "meshes/tetrahedron.off";
	
	typedef SimplexGrid<3, CgalTriangulation> Grid;
	typedef typename Grid::GlobalScene GS;
	std::shared_ptr<GS> gs(new GS(task));
	Grid grid(task, gs.get(), 0);
	
	std::vector<size_t> sequence(grid.sizeOfAllNodes());
	for (size_t it = 0; it < grid.sizeOfRealNodes(); ++it) {
		const auto neighbors = grid.findNeighborVertices(it);
		sequence[it] = neighbors.size();
	}
	
	for (int cntr = 0; cntr < 100; cntr++) {
		std::vector<size_t> parallel(grid.sizeOfAllNodes());
		#pragma omp parallel for
		for (size_t it = 0; it < grid.sizeOfRealNodes(); ++it) {
			const auto neighbors = grid.findNeighborVertices(it);
			parallel[it] = neighbors.size();
		}
		ASSERT_EQ(sequence, parallel);
	}
}


TEST(SimplexGrid2D, ThreadSafety) {
	Task task;
	real h = 0.1;
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
	typedef typename Grid::GlobalScene GS;
	std::shared_ptr<GS> gs(new GS(task));
	Grid grid(task, gs.get(), 0);
	
	std::vector<size_t> sequence(grid.sizeOfAllNodes());
	for (size_t it = 0; it < grid.sizeOfRealNodes(); ++it) {
		const auto neighbors = grid.findNeighborVertices(it);
		sequence[it] = neighbors.size();
	}
	
	for (int cntr = 0; cntr < 100; cntr++) {
		std::vector<size_t> parallel(grid.sizeOfAllNodes());
		#pragma omp parallel for
		for (size_t it = 0; it < grid.sizeOfRealNodes(); ++it) {
			const auto neighbors = grid.findNeighborVertices(it);
			parallel[it] = neighbors.size();
		}
		ASSERT_EQ(sequence, parallel);
	}
}


