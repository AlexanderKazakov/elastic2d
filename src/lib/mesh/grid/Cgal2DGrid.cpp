#include <lib/mesh/grid/Cgal2DGrid.hpp>

#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Delaunay_mesher_2.h>


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

	markInnersAndBorders();
}


std::set<Cgal2DGrid::Iterator> Cgal2DGrid::
findNeighborVertices(const Iterator& it) const {
	VertexHandle v = vertexHandles[getIndex(it)];
	std::set<Iterator> ans;
	
	auto beginFace = triangulation.incident_faces(v);
	auto faceCirculator = beginFace;
	do {
		if (faceCirculator->is_in_domain()) {
			int vLocalIndex = faceCirculator->index(v);
			ans.insert(getIterator(
					faceCirculator->vertex(faceCirculator->cw(vLocalIndex))));
			ans.insert(getIterator(
					faceCirculator->vertex(faceCirculator->ccw(vLocalIndex))));
		}
		++faceCirculator;
	} while (faceCirculator != beginFace);
	
	return ans;
}


Cgal2DGrid::Cell Cgal2DGrid::
findOwnerCell(const Iterator& it, const Real2& shift) const {
	if (shift == Real2({0, 0})) {
		auto startFace = vertexHandles[getIndex(it)]->incident_faces();
		while (!startFace->is_in_domain()) { ++startFace; }
		return createCell(it, startFace, startFace);
	}
	
	LineWalker lineWalker(this, it, shift);
	Real2 query = coordsD(it) + shift;
	FaceHandle previousFace = lineWalker.faceHandle();
	FaceHandle currentFace = previousFace;
	while (lineWalker.isValid() && !contains(currentFace, query)) {		
		previousFace = currentFace;
		currentFace = lineWalker.next();
	}
	
	return createCell(it, currentFace, previousFace);
}


Cgal2DGrid::Cell Cgal2DGrid::
locateOwnerCell(const Iterator& it, const Real2& shift) const {
	VertexHandle beginVertex = vertexHandles[getIndex(it)];
	CgalPoint2 query = beginVertex->point() + cgalVector2(shift);
	FaceHandle ownerFace = triangulation.locate(query, beginVertex->incident_faces());
	return createCell(it, ownerFace, ownerFace);
}


std::pair<Cgal2DGrid::Iterator, Cgal2DGrid::Iterator> Cgal2DGrid::
findBorderNeighbors(const Iterator& it) const {
	/// only for border nodes
	assert_true(isBorder(it));
	VertexHandle v = vertexHandles[getIndex(it)];
	
	/// find border edges as two neighbor different domain faces 
	/// moving counterclockwise around the vertex
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
	
	return {getIterator(secondBorderNeighbour), getIterator(firstBorderNeighbour)};
}


Cgal2DGrid::Iterator Cgal2DGrid::
findVertexByCoordinates(const Real2& coordinates) const {
	for (const auto& it : *this) {
		if (coordsD(it) == coordinates) {
			return it;
		}
	}
	THROW_INVALID_ARG("There isn't a vertex with such coordinates");
}


void Cgal2DGrid::
triangulate(const Task::Cgal2DGrid& task) {
	// Special "seeds" in the interior of inner cavities
	// to tell CGAL do not mesh these cavities
	std::list<CgalPoint2> listOfSeeds;
	
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
		polygon.push_back(cgalPoint2(p));
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
	++point;
	
	while(point != polygon.vertices_end()) {
		VertexHandle current = triangulation.insert(*point);
		triangulation.insert_constraint(last, current);
		last = current;
		++point;
	}
	
	triangulation.insert_constraint(last, first);
}


Cgal2DGrid::CgalPoint2 Cgal2DGrid::
findInnerPoint(const Polygon& polygon) {
	/// find inner point in polygon
	Real2 a = real2(polygon.vertex(0));
	Real2 b = real2(polygon.vertex(1));
	Real2 middle = (a + b) / 2.0;
	Real2 cross = linal::perpendicularClockwise(b - a);
	Real2 innerPoint = middle + cross;
	int n = 1; // iteration number
	while (!polygon.has_on_bounded_side(cgalPoint2(innerPoint))) {
	// watching on the both sides of the edge, getting closer on each iteration
		innerPoint = middle + cross / pow(-2, n);
		if (n > 20) {
			THROW_BAD_MESH("Exceed number of iterations in polygon inner point search");
		}
		n++;
	}
	return cgalPoint2(innerPoint);
}


void Cgal2DGrid::markInnersAndBorders() {
	/// insert indices of border vertices into borderIndices
	/// and indices of inner vertices into innerIndices
	borderIndices.clear();
	innerIndices.clear();
	
	for (const auto it : *this) {
		if (isBorder(it)) {
			borderIndices.push_back(getIndex(it));
		} else {
			innerIndices.push_back(getIndex(it));
		}
	}
	
	assert_eq(borderIndices.size() + innerIndices.size(), sizeOfAllNodes());
	LOG_DEBUG("Number of border vertices: " << borderIndices.size());
	LOG_DEBUG("Number of inner vertices: " << innerIndices.size());
}


std::vector<Cgal2DGrid::VertexHandle> Cgal2DGrid::
commonVertices(const FaceHandle& a, const FaceHandle& b) const {
	/// return common vertices of given faces
	std::vector<VertexHandle> ans;
	for (int i = 0; i < 3; i++) {
		VertexHandle candidate = a->vertex(i);
		if (b->has_vertex(candidate)) {
			ans.push_back(candidate);
		}
	}
	return ans;
}


void Cgal2DGrid::
printFace(const FaceHandle& f, const std::string& name) const {
	/// debugging helper
	SUPPRESS_WUNUSED(name);
	LOG_DEBUG("Face " << name << ":");
	for (int i = 0; i < 3; i++) {
		if (triangulation.is_infinite(f->vertex(i))) {
			LOG_DEBUG("\nINFINITE\n");
		} else {
			LOG_DEBUG(real2(f->vertex(i)->point()));
		}
	}
}








