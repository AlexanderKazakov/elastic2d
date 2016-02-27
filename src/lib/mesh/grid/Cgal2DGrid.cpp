#include <lib/mesh/grid/Cgal2DGrid.hpp>

using namespace gcm;


void Cgal2DGrid::initializeImpl(const Task& task) {
	LOG_INFO("Start initialization");
	minimalSpatialStep = task.spatialStep;
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

	recalculateMinimalSpatialStep();
	initializeImplImpl(task);
};

void Cgal2DGrid::triangulate() {
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

	std::cout << "Number of vertices: " << triangulation.number_of_vertices() << std::endl;
	std::cout << "Number of faces: " << triangulation.number_of_faces() << std::endl;
	std::cout << "Meshing the triangulation..." << std::endl;
	Mesher mesher(triangulation);
	mesher.set_seeds(listOfSeeds.begin(), listOfSeeds.end());
	Criteria meshingCriteria; 
	meshingCriteria.set_size_bound(minimalSpatialStep);
	mesher.set_criteria(meshingCriteria);
	mesher.refine_mesh();
	std::cout << "Number of vertices: " << triangulation.number_of_vertices() << std::endl;
	std::cout << "Number of faces: " << triangulation.number_of_faces() << std::endl;
};

