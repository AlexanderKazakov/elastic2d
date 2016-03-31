#include <lib/mesh/grid/Cgal2DGrid.hpp>

using namespace gcm;


Cgal2DGrid::
Cgal2DGrid(const Task& task) :
	UnstructuredGrid(task) {
	effectiveSpatialStep = task.cgal2DGrid.spatialStep;
	triangulate();
	vertexHandles.resize(triangulation.number_of_vertices());
	size_t vertexIndex = 0;
	for (auto it = triangulation.finite_vertices_begin();
	     it != triangulation.finite_vertices_end(); it++) {
		auto handle = it->handle();
		vertexHandles[vertexIndex] = handle;
		verticesIndices.insert({handle, vertexIndex});
		vertexIndex++;
	}
}


void Cgal2DGrid::
triangulate() {
	VertexHandle va = triangulation.insert(Point(-2, -1));
	VertexHandle vb = triangulation.insert(Point(-1, 0));
	VertexHandle vc = triangulation.insert(Point(0, -1));
	VertexHandle vd = triangulation.insert(Point(-1, -2));
	triangulation.insert_constraint(va, vb);
	triangulation.insert_constraint(vb, vc);
	triangulation.insert_constraint(vc, vd);
	triangulation.insert_constraint(vd, va);

	va = triangulation.insert(Point(3, 3));
	vb = triangulation.insert(Point(-3, 3));
	vc = triangulation.insert(Point(-3, -3));
	vd = triangulation.insert(Point(3, -3));
	VertexHandle ve = triangulation.insert(Point(2, 2));
	triangulation.insert_constraint(va, vb);
	triangulation.insert_constraint(vb, vc);
	triangulation.insert_constraint(vc, vd);
	triangulation.insert_constraint(vd, ve);
	triangulation.insert_constraint(ve, va);

	std::list<Point> listOfSeeds;
	listOfSeeds.push_back(Point(-1, -1));

	LOG_DEBUG("Number of vertices: " << triangulation.number_of_vertices());
	LOG_DEBUG("Number of faces: " << triangulation.number_of_faces());
	LOG_DEBUG("Meshing the triangulation...");
	Mesher mesher(triangulation);
	mesher.set_seeds(listOfSeeds.begin(), listOfSeeds.end());
	Criteria meshingCriteria;
	assert_gt(effectiveSpatialStep, 0);
	meshingCriteria.set_size_bound(effectiveSpatialStep);
	mesher.set_criteria(meshingCriteria);
	mesher.refine_mesh();
	LOG_DEBUG("Number of vertices: " << triangulation.number_of_vertices());
	LOG_DEBUG("Number of faces: " << triangulation.number_of_faces());
}


