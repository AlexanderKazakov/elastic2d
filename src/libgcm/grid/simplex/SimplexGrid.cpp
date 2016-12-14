#include <libgcm/grid/simplex/SimplexGrid.hpp>

#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>
#include <libgcm/grid/simplex/cgal/LineWalker.hpp>


using namespace gcm;


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
SimplexGrid<Dimensionality, TriangulationT>::
SimplexGrid(const GridId id_, const ConstructionPack& constructionPack) :
		UnstructuredGrid(id_),
		triangulation(constructionPack.triangulation) {
	
	assert_ne(id, EmptySpaceFlag);
	LOG_INFO("Start construction of the grid " << id << " ...");
	
	/// find local cells and vertices in global triangulation
	std::set<VertexHandle> localVertices;
	for (auto cellIter  = triangulation->allCellsBegin();
	          cellIter != triangulation->allCellsEnd(); ++cellIter) {
		if (cellIter->info().getGridId() == id) {
			
			// check on mesher artifacts
			if (Triangulation::minimalCellHeight(cellIter) < EQUALITY_TOLERANCE) {
				// FIXME - EQUALITY_TOLERANCE is bad solution, move to mesher?
				cellIter->info().setGridId(EmptySpaceFlag);
				triangulation->printCell(cellIter, std::string(
						"replaced as degenerate with minimalHeight == ") +
						std::to_string(Triangulation::minimalCellHeight(cellIter)));
//				THROW_BAD_MESH("Fix degenerate cells");
				continue;
			}
			
			cellHandles.push_back(cellIter);
			for (int i = 0; i < CELL_POINTS_NUMBER; i++) {
				localVertices.insert(cellIter->vertex(i));
			}
		}
	}
	vertexHandles.assign(localVertices.begin(), localVertices.end());
	
	/// write local vertices indices to vertices info
	/// note (!) later, it will be invalidated by other grids constructors,
	/// so it is just temporary auxiliary info
	for (size_t i = 0; i < vertexHandles.size(); i++) {
		vertexHandles[i]->info() = i;
	}
	
	/// write local vertices indices to cells info
	/// this is not temporary, because cell belongs to the only one grid
	for (CellIterator cell = cellBegin(); cell != cellEnd(); ++cell) {
		for (int i = 0; i < CELL_POINTS_NUMBER; i++) {
			(*cell)->info().localVertexIndices[i] = (*cell)->vertex(i)->info();
		}
	}
	// TODO - write border/contact/inner state into vertices info 
	// because this information is equal for all grids in contact
	markInnersAndBorders();
	calculateMinimalSpatialStep();
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
typename SimplexGrid<Dimensionality, TriangulationT>::Cell
SimplexGrid<Dimensionality, TriangulationT>::
findCellCrossedByTheRay(const Iterator& it, const RealD& shift) const {
	/// Crucial requirements for LineWalker:
	/// possible return values of LineWalker::cellsAlongSegment:
	/// - empty (no cells found)
	/// - all cells belong to the grid (needed cell seem to be found)
	/// - the number of cells greater than 1 and (only!) the last one
	///   is outside the grid -- the ray seems to cross the border
	typedef LineWalker<Triangulation, DIMENSIONALITY> LINE_WALKER;
	RealD start = coordsD(it);
	RealD query = start + shift;
	const auto isLocalCell = [=](const CellHandle c) {
		return belongsToTheGrid(c);
	};
	
//	if (linal::length(start - RealD({0.87696,0.172084,0.87696})) < 0.01 &&
//			linal::length(query - RealD({1.11266,0.172084,1.11266})) < 0.01) {
//		std::cout << "Here!";
//	}
	
	/// Try to go along the ray from it to it+shift cell-by-cell
	std::vector<CellHandle> cellsAlong = LINE_WALKER::cellsAlongSegment(
			triangulation, isLocalCell, vertexHandle(it), query);
//	int i = 0;
//	LOG_INFO("Start search from vertex");
//	for (CellHandle ch : cellsAlong) {
//		triangulation->printCell(ch, std::to_string(i++));
//	}
	Cell foundCell = checkLineWalkFoundCell(it, cellsAlong, start, query);
	if (foundCell.n > 0) { return foundCell; }
	
	/// Now, two situations are possible:
	/// - the ray is going out of the grid from a border (contact) node,
	/// - we have some numerical inexactness in geometrical operations.
	/// In order to avoid inexactness try to start from the center of
	/// a cell incident to vertexHandle(it) instead of the vertex itself
	CellHandle startCell = triangulation->findCrossedIncidentCell(
			isLocalCell, vertexHandle(it), query, 0);
	if (startCell == NULL) {
		startCell = triangulation->findCrossedIncidentCell(
				isLocalCell, vertexHandle(it), query, EQUALITY_TOLERANCE);
	}
	if (startCell == NULL) {
		startCell = localIncidentCells(it).front();
	}
	cellsAlong = LINE_WALKER::cellsAlongSegment(
			triangulation, isLocalCell, startCell, query);
//	i = 0;
//	LOG_INFO("Start search from cell center");
//	for (CellHandle ch : cellsAlong) {
//		triangulation->printCell(ch, std::to_string(i++));
//	}
	foundCell = checkLineWalkFoundCell(it, cellsAlong, start, query);
	if (foundCell.n > 0) { return foundCell; }
	
	/// Now, the majority (unlikely all) of inexactness cases are handled.
	/// Also, the case when the ray starts from the border node and 
	/// then goes outside the grid through another border 
	/// (firstly going through the inner cells) is not handled.
	/// However, by now, it's not handled by the gcm-method too.
	/// If somewhen fix it, remember the false cases -- immediate outgoing 
	/// from a border node (see old commits before 01.12.2016 for algo).
	assert_false(isInner(it));
	return createCell();
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
typename SimplexGrid<Dimensionality, TriangulationT>::Cell
SimplexGrid<Dimensionality, TriangulationT>::
checkLineWalkFoundCell(const Iterator& it,
		const std::vector<CellHandle>& cellsAlong,
		const RealD& start, const RealD& query) const {
	if (cellsAlong.empty()) { return createCell(); }
	
	/// Handle cases when the ray seems to hit into the inner space of the grid.
	/// Even when some cell belongs to the grid we have to check if that cell
	/// contains query point (due to numerical inexactness in line walk search)
	CellHandle last = cellsAlong.back();
	if (belongsToTheGrid(last) &&
			Triangulation::contains(last, query, EQUALITY_TOLERANCE)) {
		return createCell(last);
	}
	
	if (cellsAlong.size() == 1) {
		assert_false(isInner(it));
		return createCell();
	}
	
	/// Handle cases when the ray seems to hit into the inner space near
	/// the border, but due to numerical inexactness the last cell is outside
	CellHandle prev = *std::next(cellsAlong.rbegin());
	if (Triangulation::contains(prev, query, EQUALITY_TOLERANCE)) {
		assert_true(belongsToTheGrid(prev));
		return createCell(prev);
	}
	
	if (!isInner(it)) {
		return createCell();
	}
	
	/// Handle cases when the ray seems to hit outside the grid 
	/// going from an inner node through a border facet.
	/// Again, we must explicitly check that the border facet contains
	/// intersection of the ray and the border facet
	assert_true(belongsToTheGrid(prev));
	if (!belongsToTheGrid(last)) {
		std::vector<VertexHandle> common = Triangulation::filterFaceNotCrossedByTheRay(
				Triangulation::commonVertices(prev, last),
				start, query, localEqualityTolerance());
		return Cell(common.begin(), common.end(), [=](const VertexHandle vh) {
				return localVertexIndex(vh, prev); });
	}
	
	return createCell();
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
typename SimplexGrid<Dimensionality, TriangulationT>::Iterator
SimplexGrid<Dimensionality, TriangulationT>::
findVertexByCoordinates(const RealD& coordinates) const {
	for (const auto& it : *this) {
		if (coordsD(it) == coordinates) {
			return it;
		}
	}
	THROW_INVALID_ARG("There isn't a node with such coordinates");
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void
SimplexGrid<Dimensionality, TriangulationT>::
markInnersAndBorders() {
/// insert indices of contact vertices into contactIndices
/// insert indices of border vertices into borderIndices
/// and indices of inner vertices into innerIndices
	contactIndices.clear();
	borderIndices.clear();
	innerIndices.clear();
	size_t multicontactCounter = 0;
	
	for (const auto it : *this) {
		switch (borderState(it)) {
			case BorderState::CONTACT:
				contactIndices.push_back(it);
				break;
				
			case BorderState::BORDER:
				borderIndices.push_back(it);
				break;
				
			case BorderState::INNER:
				innerIndices.push_back(it);
				break;
				
			case BorderState::MULTICONTACT:
				++multicontactCounter;
				break;
				
			default:
				THROW_BAD_MESH("Unknown border state");
				break;
		}
	}
	
	assert_eq(contactIndices.size() + borderIndices.size() + innerIndices.size() +
			multicontactCounter, sizeOfAllNodes());
	
	LOG_INFO("Number of multicontact vertices: " << multicontactCounter);
	LOG_INFO("Number of contact vertices: " << contactIndices.size());
	LOG_INFO("Number of border vertices: " << borderIndices.size());
	LOG_INFO("Number of inner vertices: " << innerIndices.size());
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void
SimplexGrid<Dimensionality, TriangulationT>::
calculateMinimalSpatialStep() {
	real summHeight = 0;
	size_t cellsNumber = 0;
	minimalSpatialStep = std::numeric_limits<real>::max();
	
	for (auto cell = cellBegin(); cell != cellEnd(); ++cell) {
		real h = Triangulation::minimalCellHeight(*cell);
		if (minimalSpatialStep > h) {
			minimalSpatialStep = h;
		}
		
		real summPrevious = summHeight;
		summHeight += h;
		assert_ne(summPrevious, summHeight); // check overflow and zero
		++cellsNumber;
	}
	
	averageSpatialStep = summHeight / (real)cellsNumber;
	LOG_INFO("Minimal height: " << minimalSpatialStep
	    << ", Average height: " << averageSpatialStep);
}



template class SimplexGrid<2, CgalTriangulation>;
template class SimplexGrid<3, CgalTriangulation>;

