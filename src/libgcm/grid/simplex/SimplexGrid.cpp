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
//				triangulation->printCell(cellIter, std::string(
//						"replaced as degenerate with minimalHeight == ") +
//						std::to_string(Triangulation::minimalCellHeight(cellIter)));
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
	
	markInnersAndBorders();
	calculateMinimalSpatialStep();
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
typename SimplexGrid<Dimensionality, TriangulationT>::RealD
SimplexGrid<Dimensionality, TriangulationT>::
contactNormal(const Iterator& it, const GridId neighborId) const {

	std::list<RealD> facesNormals;
	
	const std::list<CellHandle> localCells = localIncidentCells(it);
	for (const CellHandle localCell : localCells) {
		
		for (int i = 0; i < CELL_POINTS_NUMBER; i++) {
			CellHandle outerCell = localCell->neighbor(i);
			
			if (
			/// localCell's neighbor with neighborId ..
				outerCell->info().getGridId() == neighborId &&
			/// .. which also has our vertex
				Utils::has(Triangulation::commonVertices(localCell, outerCell),
						vertexHandle(it)) ) {
				
				facesNormals.push_back(
						Triangulation::contactNormal(localCell, outerCell));
			}
		}
	}
	
	if (facesNormals.empty()) { return RealD::Zeros(); };
	
	return linal::normalize(std::accumulate(
			facesNormals.begin(), facesNormals.end(), RealD::Zeros()));
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
typename SimplexGrid<Dimensionality, TriangulationT>::Cell
SimplexGrid<Dimensionality, TriangulationT>::
findOwnerCell(const Iterator& it, const RealD& shift) const {
	
	/// Struct allows to go along the line cell-by-cell in triangulation
	typedef LineWalker<Triangulation, DIMENSIONALITY>      LINE_WALKER;
	LINE_WALKER lineWalker(triangulation, vertexHandle(it), shift, id);
	RealD query = coordsD(it) + shift;
	
	CellHandle current = lineWalker.currentCell();
	CellHandle previous = current;
	if (current == NULL) { return createCell(); }
	
	while ( belongsToTheGrid(current) &&
	        !Triangulation::contains(current, query) ) {
	// along the line until found or go outside the grid
		previous = current;
		current = lineWalker.next();
	}
	
	
	if ( belongsToTheGrid(current) ) {
	// cell found inside the grid
		assert_true(belongsToTheGrid(previous));
		return createCell(current);
	}
	
	if ( belongsToTheGrid(previous) ) {
	// seems to go outside the grid through the border
		Cell ans = createCell(
				Triangulation::commonVertices(current, previous), previous);
		if ( isInner(it) || !ans.has(it) ) {
		// it is not immediate outgoing from a border node
			return ans;
		}
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
				
			default:
				THROW_BAD_MESH("Unknown border state");
				break;
		}
	}
	
	assert_eq(contactIndices.size() + borderIndices.size() +
			innerIndices.size(), sizeOfAllNodes());
	
	LOG_DEBUG("Number of contact vertices: " << contactIndices.size());
	LOG_DEBUG("Number of border vertices: " << borderIndices.size());
	LOG_DEBUG("Number of inner vertices: " << innerIndices.size());
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

