#include <libgcm/grid/simplex/SimplexGrid.hpp>

#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>
#include <libgcm/grid/simplex/cgal/LineWalker.hpp>
#include <libgcm/util/snapshot/VtkUtils.hpp>


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
//				cellIter->info().setGridId(EmptySpaceFlag);
//				THROW_BAD_MESH("Fix degenerate cells");
//				continue;
				triangulation->printCell(cellIter, std::string(
						"is degenerate with minimalHeight == ") +
						std::to_string(Triangulation::minimalCellHeight(cellIter)));
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
	
	bool here = false;
	if (Dimensionality == 3) {
		here =
			linal::length(start - RealD({-1.67878, -3.05063, 143.62})) < getAverageHeight() / 10 &&
			linal::length(query - RealD({-1.74938, -3.05063, 143.62})) < getAverageHeight() / 10;
	}
	if (here) {
		std::cout << "Here!" << std::endl;
	}
	
	/// Try to go along the ray from it to it+shift cell-by-cell
	std::vector<CellHandle> cellsAlong = LINE_WALKER::cellsAlongSegment(
			triangulation, isLocalCell, vertexHandle(it), query);
	if (here) {
		int i = 0;
		LOG_INFO("Start search from vertex");
		vtk_utils::drawSegmentToVtk(start, query, "lineFromVertex");
		writeCellsToVtk(cellsAlong, "fromVertex");
		for (CellHandle ch : cellsAlong) {
			triangulation->printCell(ch, std::to_string(i++));
			std::cout << "isLocalCell(ch) == " << isLocalCell(ch) << std::endl;
			std::cout << "height(ch) == " << Triangulation::minimalCellHeight(ch) << std::endl;
		}
	}
	Cell foundCell = checkLineWalkFoundCell(it, cellsAlong, start, query);
	
	if (foundCell.n > 0) { return foundCell; }
	
	/// Now, two situations are possible:
	/// - the ray is going out of the grid from a border (contact) node,
	/// - we have some numerical inexactness in geometrical operations.
	/// In order to avoid inexactness try to start from inside of
	/// the incident cell which is crossed by the search direction
	CellHandle startCell = triangulation->findCrossedIncidentCell(
			isLocalCell, vertexHandle(it), query, 0);
	if (startCell == NULL) {
		startCell = triangulation->findCrossedIncidentCell(
				isLocalCell, vertexHandle(it), query, EQUALITY_TOLERANCE);
	}
	if (startCell != NULL) {
		const int ITERATION_MAX = (int)(-log2(EQUALITY_TOLERANCE));
		for (int iteration = 0; iteration < ITERATION_MAX; iteration++) {
			real w = std::pow(2, -iteration);
			RealD startPoint = Triangulation::center(startCell) * w + start * (1 - w);
			cellsAlong = LINE_WALKER::cellsAlongSegment(
					triangulation, isLocalCell, startCell, startPoint, query);
			if (here) {
				vtk_utils::drawSegmentToVtk(startPoint, query,
						"lineFromCenter" + std::to_string(iteration));
				writeCellsToVtk(cellsAlong, "fromCenter" + std::to_string(iteration));
			}
			foundCell = checkLineWalkFoundCell(it, cellsAlong, start, query);
			if (foundCell.n > 0) { return foundCell; }
		}
	}
	
	/// Now, try to start from centers of all local incident cells
	for (CellHandle incidentCell : localIncidentCells(it)) {
		cellsAlong = LINE_WALKER::cellsAlongSegment(
				triangulation, isLocalCell, incidentCell,
						Triangulation::center(incidentCell), query);
		if (here) {
			int i = 0;
			LOG_INFO("Continue search all cells centers");
			for (CellHandle ch : cellsAlong) {
				triangulation->printCell(ch, std::to_string(i++));
				std::cout << "isLocalCell(ch) == " << isLocalCell(ch) << std::endl;
			}
		}
		foundCell = checkLineWalkFoundCell(it, cellsAlong, start, query);
		if (foundCell.n > 0) { return foundCell; }
	}
	
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
				start, query, EQUALITY_TOLERANCE);
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
//		real summPrevious = summHeight;
		summHeight += h;
		// by now, zero cells are allowed
//		assert_ne(summPrevious, summHeight); // check overflow and zero cells
		++cellsNumber;
	}
	
	averageSpatialStep = summHeight / (real)cellsNumber;
	LOG_INFO("Minimal height: " << minimalSpatialStep
	    << ", Average height: " << averageSpatialStep);
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void SimplexGrid<Dimensionality, TriangulationT>::
writeCellsToVtk(const std::vector<CellHandle>& cells, const std::string filename) const {
	std::vector<elements::Element<RealD, CELL_POINTS_NUMBER>> elems;
	for (const CellHandle ch : cells) {
		assert_false(triangulation->isInfinite(ch));
		elements::Element<RealD, CELL_POINTS_NUMBER> elem;
		elem.n = CELL_POINTS_NUMBER;
		for (int i = 0; i < CELL_POINTS_NUMBER; i++) {
			elem(i) = Triangulation::realD(ch->vertex(i));
		}
		elems.push_back(elem);
	}
	vtk_utils::drawCellsToVtk(elems, filename);
}

template class SimplexGrid<2, CgalTriangulation>;
template class SimplexGrid<3, CgalTriangulation>;

