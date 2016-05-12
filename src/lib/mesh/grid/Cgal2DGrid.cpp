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


Cgal2DGrid::Triangle Cgal2DGrid::
findOwnerTriangle(const Iterator& it, const Real2& shift) const {
	/// zero shift not handled to reduce if-else hell
	assert_true(shift != Real2({0, 0}));
	
	Triangle ans;
	ans.inner = false;
		
	VertexHandle startVertex = vertexHandles[getIndex(it)];
	CgalPoint2 query = startVertex->point() + cgalVector2(shift);
	auto startFace = startVertex->incident_faces();
	auto lineWalker = triangulation.line_walk(startVertex->point(), query, startFace);
	
	if (!lineWalker.is_empty()) {
	// normal situation; go along the line until found or come out of the body
		while (lineWalker->is_in_domain() &&
		       triangulation.triangle(lineWalker).has_on_unbounded_side(query)) {
			++lineWalker;
		}
		
	} else {
	// LineFaceCirculator recognizes tangent faces if they are at the left 
	// of the line only; try to change line direction
		lineWalker = triangulation.line_walk(query, startVertex->point());
		
		if (!lineWalker.is_empty()) {
		// the line is tangent to border yet; go along the line in the inverse direction
			auto lineWalkerBegin = lineWalker;
			while ( !(lineWalker->has_vertex(startVertex) && lineWalker->is_in_domain()) ) {
				--lineWalker;
				if (lineWalker == lineWalkerBegin) {
				// some convex corner
					return ans;
				}
			}
			while (lineWalker->is_in_domain() &&
			       triangulation.triangle(lineWalker).has_on_unbounded_side(query)) {
				--lineWalker;
			}
			
		} else {
		// really empty (some convex corner too)
			return ans;
		}
		
	}
	
	if (lineWalker->is_in_domain()) {
		ans.inner = true;
		for (int i = 0; i < 3; i++) {
			ans.p[i] = getIterator(lineWalker->vertex(i));
		}
	}
	
	// TODO - handle the case of exact hit to the border edge, when outer face is chosen
	return ans;
}


Cgal2DGrid::Triangle Cgal2DGrid::
locateOwnerTriangle(const Iterator& it, const Real2& shift) const {
	Triangle ans;
	ans.inner = false;

	VertexHandle beginVertex = vertexHandles[getIndex(it)];
	CgalPoint2 query = beginVertex->point() + cgalVector2(shift);
	FaceHandle ownerFace = triangulation.locate(query, beginVertex->incident_faces());
	
	if (ownerFace->is_in_domain()) {
		ans.inner = true;
		for (int i = 0; i < 3; i++) {
			ans.p[i] = getIterator(ownerFace->vertex(i));
		}
	}
	// TODO - handle the case of exact hit to the border edge, when outer face is chosen
	
	return ans;
}


std::pair<Cgal2DGrid::Iterator, Cgal2DGrid::Iterator> Cgal2DGrid::
findCrossingBorder(const Iterator& start, const Real2& direction) const {
	VertexHandle v = vertexHandles[getIndex(start)];
	FaceHandle startFace = findCrossedFace(v, direction);
	return findCrossingBorder(v, startFace, direction);
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
findBorderFlexion(Iterator first, Iterator second) const {
	
	auto findThird = [&]() {
		auto neighbors = findBorderNeighbors(second);
		return (first != neighbors.first) ? neighbors.first
		                                  : neighbors.second;
	};
	
	Iterator third = findThird();
	while (linal::isDegenerate(coords2d(first), coords2d(second), coords2d(third))) {
		first = second;
		second = third;
		third = findThird();
	}
	
	return second;
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
	
	for (auto v  = triangulation.finite_vertices_begin();
	          v != triangulation.finite_vertices_end(); ++v) {

		bool isBorderNode = false;
		auto beginFace = triangulation.incident_faces(v);
		auto faceCirculator = beginFace;
		do {
			if ( !faceCirculator->is_in_domain() ) {
				isBorderNode = true;
				break;
			}
			++faceCirculator;
		} while (faceCirculator != beginFace);
		
		if (isBorderNode) {
			borderIndices.insert(getIndex(getIterator(v)));
		} else {
			innerIndices.insert(getIndex(getIterator(v)));
		}
		
	}
	
	assert_eq(borderIndices.size() + innerIndices.size(), sizeOfAllNodes());
	LOG_DEBUG("Number of border vertices: " << borderIndices.size());
	LOG_DEBUG("Number of inner vertices: " << innerIndices.size());
}


Cgal2DGrid::FaceHandle Cgal2DGrid::
findCrossedFace(const VertexHandle start, const Real2& direction) const {
	/// Choose among incident faces of the given point that one which is
	/// crossed by the line from that point in specified direction.

	auto faceCirculatorBegin = triangulation.incident_faces(start);
	while (triangulation.is_infinite(faceCirculatorBegin)) { ++faceCirculatorBegin; }
	int n = 1; // iteration number
	FaceHandle ans = NULL;
	while (ans == NULL) {
		auto faceCirculator = faceCirculatorBegin;
		do {
			CgalPoint2 q = start->point() + cgalVector2(direction / pow(2, n));
			if (!triangulation.triangle(faceCirculator).has_on_unbounded_side(q)) {
				ans = faceCirculator;
				break;
			}
			do { ++faceCirculator; } while (triangulation.is_infinite(faceCirculator));
		} while (faceCirculator != faceCirculatorBegin);
		
		++n; 
		if (n > 20) {
			THROW_BAD_MESH("Exceed number of iterations in findCrossedFace");
		}
	}
	
	return ans;	
}


Cgal2DGrid::VertexHandle Cgal2DGrid::
commonVertex(const FaceHandle& a, const FaceHandle& b) const {
	/// return some common vertex of given faces
	/// or throw Exception if there isn't such
	for (int i = 0; i < 3; i++) {
		VertexHandle candidate = a->vertex(i);
		if (b->has_vertex(candidate)) {
			return candidate;
		}
	}
	THROW_INVALID_ARG("There aren't common vertices");
}


std::pair<Cgal2DGrid::Iterator, Cgal2DGrid::Iterator> Cgal2DGrid::
findCrossingBorder(const VertexHandle& v, const FaceHandle& f,
		const Real2& direction) const {
	/// given a start point and direction, go along this line until
	/// intersection with border;
	/// @note face is incident for given vertex,
	/// given line must go across given face
	/// @return found border as pair of its vertices
	
	assert_true(f->is_in_domain());
	
	auto lineWalker = triangulation.line_walk(
			v->point(), v->point() + cgalVector2(direction), f);
	while (lineWalker->is_in_domain()) {
	// go along the line until come out of the body
		++lineWalker;
	}
	
	FaceHandle outerFace = lineWalker;
	FaceHandle innerFace = --lineWalker;
	
	Iterator first;
	Iterator second;
	int neighborIndex = -1;
	if (innerFace->has_neighbor(outerFace, neighborIndex)) {
	// line is crossing common edge
		first = getIterator(innerFace->vertex(innerFace->cw(neighborIndex)));
		second = getIterator(innerFace->vertex(innerFace->ccw(neighborIndex)));
	} else {
	// line hits exact to unique common vertex
		first = getIterator(commonVertex(outerFace, innerFace));
		second = findBorderNeighbors(first).first;
	}
	
	return {first, second};
}


void Cgal2DGrid::
printFace(const FaceHandle& f, const std::string& name) const {
	/// debugging helper
	LOG_DEBUG("Face " << name << ":");
	for (int i = 0; i < 3; i++) {
		if (triangulation.is_infinite(f->vertex(i))) {
			LOG_DEBUG("\nINFINITE\n");
		} else {
			LOG_DEBUG(real2(f->vertex(i)->point()));
		}
	}
}








