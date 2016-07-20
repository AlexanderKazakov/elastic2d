#include <lib/mesh/grid/cgal/Cgal2DGrid.hpp>
#include <lib/mesh/grid/cgal/Cgal2DLineWalker.hpp>

#include <libcgalmesh/Cgal2DMesher.hpp>

using namespace gcm;


Cgal2DGrid::
Cgal2DGrid(const Task& task) :
	UnstructuredGrid(task),
	effectiveSpatialStep(task.cgal2DGrid.spatialStep),
	movable(task.cgal2DGrid.movable) {
	
	// convert task to cgal mesher format
	typedef cgalmesh::Cgal2DMesher::Body Body;
	std::vector<Body> bodies;
	for (const auto& b : task.cgal2DGrid.bodies) {
		bodies.push_back({b.outer, b.inner});
	}
	
	LOG_DEBUG("Call Cgal2DMesher");
	cgalmesh::Cgal2DMesher::triangulate(
			task.cgal2DGrid.spatialStep, bodies, triangulation);
	LOG_DEBUG("Number of vertices after meshing: " << triangulation.number_of_vertices());
	LOG_DEBUG("Number of faces after meshing: " << triangulation.number_of_faces());
	
	
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
	VertexHandle v = vertexHandle(it);
	std::set<Iterator> ans;
	
	auto beginFace = triangulation.incident_faces(v);
	auto faceCirculator = beginFace;
	do {
		if (isInDomain(faceCirculator)) {
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
		auto startFace = vertexHandle(it)->incident_faces();
		while (!isInDomain(startFace)) { ++startFace; }
		return createCell(it, startFace, startFace);
	}
	
	Cgal2DLineWalker lineWalker(this, it, shift);
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
	VertexHandle beginVertex = vertexHandle(it);
	CgalPoint2 query = beginVertex->point() + cgalVector2(shift);
	FaceHandle ownerFace = triangulation.locate(query, beginVertex->incident_faces());
	return createCell(it, ownerFace, ownerFace);
}


std::pair<Cgal2DGrid::Iterator, Cgal2DGrid::Iterator> Cgal2DGrid::
findBorderNeighbors(const Iterator& it) const {
	/// only for border nodes
	assert_true(isBorder(it));
	VertexHandle v = vertexHandle(it);
	
	/// find border edges as two neighbor different domain faces 
	/// moving counterclockwise around the vertex
	auto firstFace = triangulation.incident_faces(v);
	auto secondFace = firstFace; ++secondFace;
	
	while ( !(  isInDomain(firstFace) && 
	           !isInDomain(secondFace) ) ) {
		++firstFace; ++secondFace;
	}
	VertexHandle firstBorderNeighbour = 
			firstFace->vertex(firstFace->cw(firstFace->index(v)));
	assert_true(secondFace->has_vertex(firstBorderNeighbour));
	
	while ( !( !isInDomain(firstFace) && 
	            isInDomain(secondFace) ) ) {
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




