#include <lib/mesh/grid/Cgal2DGrid.hpp>

using namespace gcm;


Cgal2DGrid::
Cgal2DGrid(const Task& task) :
	UnstructuredGrid(task),
	effectiveSpatialStep(task.cgal2DGrid.spatialStep), 
	movable(task.cgal2DGrid.movable) {

	triangulate(task.cgal2DGrid);
	
	vertexHandles.resize(triangulation.number_of_vertices());
	size_t vertexIndex = 0;
	for (auto it = triangulation.finite_vertices_begin();
	     it != triangulation.finite_vertices_end(); it++) {
		auto handle = it->handle();
		handle->info() = vertexIndex;
		vertexHandles[vertexIndex] = handle;
		vertexIndex++;
	}

	markBorders();
}


Real2 Cgal2DGrid::
normal(const Iterator& it) const {
	// only for border nodes
	assert_true(borderIndices.find(it.iter) != borderIndices.end());

	VertexHandle v = vertexHandles[it.iter];
	
	// find border edges as two neighbor different domain faces 
	// moving counterclockwise around the vertex
	auto firstFace = triangulation.incident_faces(v);
	auto secondFace = firstFace; ++secondFace;
	
	while ( !(  firstFace->is_in_domain() && 
	           !secondFace->is_in_domain() ) ) {
		++firstFace; ++secondFace;
	}
	
	VertexHandle firstBorderNeighbour = 
			firstFace->vertex(firstFace->cw(firstFace->index(v)));
	assert_true(secondFace->has_vertex(firstBorderNeighbour));
	
	while ( !( !firstFace->is_in_domain() && 
	            secondFace->is_in_domain() ) ) {
		++firstFace; ++secondFace;
	}
	
	VertexHandle secondBorderNeighbour = 
			firstFace->vertex(firstFace->cw(firstFace->index(v)));
	assert_true(secondFace->has_vertex(secondBorderNeighbour));
	
	Real2 borderVector = PointToReal2(secondBorderNeighbour->point()) -
						 PointToReal2(firstBorderNeighbour->point());
	
	return linal::normalize(
	       linal::perpendicularClockwise(borderVector));
}


void Cgal2DGrid::
triangulate(const Task::Cgal2DGrid& task) {
	// Special "seeds" in the interior of inner cavities
	// to tell CGAL do not mesh these cavities
	std::list<Point> listOfSeeds;
	
	// insert all task outer borders and inner cavities
	for (const auto& body : task.bodies) {
		// TODO - check bodies and cavities intersections
		insertPolygon(makePolygon(body.outer));
		for (const auto& innerCavity : body.inner) {
			insertPolygon(makePolygon(innerCavity));
			listOfSeeds.push_back(findInnerPoint(makePolygon(innerCavity)));
		}
	}

	LOG_DEBUG("Number of vertices before meshing: " << triangulation.number_of_vertices());
	LOG_DEBUG("Number of faces before meshing: " << triangulation.number_of_faces());
	LOG_DEBUG("Meshing the triangulation...");
	
	Mesher mesher(triangulation);
	mesher.set_seeds(listOfSeeds.begin(), listOfSeeds.end());
	
	Criteria meshingCriteria;
	assert_gt(effectiveSpatialStep, 0);
	meshingCriteria.set_size_bound(effectiveSpatialStep);
	mesher.set_criteria(meshingCriteria);
	
	mesher.refine_mesh(); // meshing
	
	LOG_DEBUG("Number of vertices after meshing: " << triangulation.number_of_vertices());
	LOG_DEBUG("Number of faces after meshing: " << triangulation.number_of_faces());
}


Cgal2DGrid::Polygon Cgal2DGrid::
makePolygon(const std::vector<Real2>& points) {
	/// convert point set to simple CGAL polygon
	assert_ge(points.size(), 3); // polygon is a closed line
	
	Polygon polygon;
	for (const auto& p : points) {
		polygon.push_back(Real2ToPoint(p));
	}
	assert_true(polygon.is_simple()); // has not intersections
	
	return polygon;
}


void Cgal2DGrid::
insertPolygon(const Polygon& polygon) {
	/// insert points and constraints - lines between them - to triangulation
	auto point = polygon.vertices_begin();
	VertexHandle first = triangulation.insert(*point);
	VertexHandle last = first;
	point++;
	
	while(point != polygon.vertices_end()) {
		VertexHandle current = triangulation.insert(*point);
		triangulation.insert_constraint(last, current);
		last = current;
		point++;
	}
	
	triangulation.insert_constraint(last, first);
}


Cgal2DGrid::Point Cgal2DGrid::
findInnerPoint(const Polygon& polygon) {
	/// find inner point in polygon
	Real2 a = PointToReal2(polygon.vertex(0));
	Real2 b = PointToReal2(polygon.vertex(1));
	Real2 middle = (a + b) / 2.0;
	Real2 cross = linal::perpendicularClockwise(b - a);
	Real2 innerPoint = middle + cross;
	int n = 1; // iteration number
	while (!polygon.has_on_bounded_side(Real2ToPoint(innerPoint))) {
	// watching on the both sides of the edge, getting closer on each iteration
		innerPoint = middle + cross / pow(-2, n);
		if (n > 20) {
			THROW_BAD_MESH("Something is wrong with polygon inner point search");
		}
		n++;
	}
	return Real2ToPoint(innerPoint);
}


void Cgal2DGrid::markBorders() {
	/// insert indices of border vertices into borderIndices
	borderIndices.clear();

	for (auto cell = triangulation.all_faces_begin();
		 cell != triangulation.all_faces_end(); ++cell) {
		
		if (!cell->is_in_domain()) {
			for (int i = 0; i < 3; i++) {
				VertexHandle v = cell->vertex(i);
				if (!triangulation.is_infinite(v)) {
					borderIndices.insert(v->info());
				}
			}
		}
	}
	LOG_DEBUG("Number of border vertices: " << borderIndices.size());
}





