#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>
#include <libgcm/grid/simplex/VertexInfoAndCellInfo.hpp>
#include <libgcm/util/snapshot/VtkUtils.hpp>

using namespace gcm;


template<int Dimensionality, typename VertexInfo, typename CellInfo>
CgalTriangulation<Dimensionality, VertexInfo, CellInfo>::
CgalTriangulation(const Task& task) : Base(task) {
	static_assert(CELL_POINTS_NUMBER == Base::CELL_SIZE, "");
	static_assert(DIMENSIONALITY == Base::DIMENSIONALITY, "");
	rescale(task.simplexGrid.scale);
	
	int hangedCellsCounter = 1;
	int disconnectedCellsSetsCaseCounter = 1;
	int iterationsCounter = 0;
	while (disconnectedCellsSetsCaseCounter > 0 || hangedCellsCounter > 0) {
		if (++iterationsCounter > 10 ) { break; }
		LOG_INFO("Cleaning the triangulation from bad material cases: "
				<< "iteration " << iterationsCounter);
		
		hangedCellsCounter = correctHangedCells();
		LOG_INFO(hangedCellsCounter << " hanged cells cases was found");
		
		
		disconnectedCellsSetsCaseCounter = 0;
		for (auto vertexIter = verticesBegin();
		          vertexIter != verticesEnd(); ++vertexIter) {
			while (clearFromDisconnectedCellSets(vertexIter)) {
				++disconnectedCellsSetsCaseCounter;
			}
		}
		LOG_INFO(disconnectedCellsSetsCaseCounter << " DCS cases was found");
	}
}


template<int Dimensionality, typename VertexInfo, typename CellInfo>
int
CgalTriangulation<Dimensionality, VertexInfo, CellInfo>::
correctHangedCells() {
	int hangsCounter = 0;
	for (auto cell = this->allCellsBegin();
			cell != this->allCellsEnd(); ++cell) if (!isInfinite(cell)) {
		std::multiset<GridId> neighborsId;
		for (int i = 0; i < CELL_POINTS_NUMBER; i++) {
			neighborsId.insert(cell->neighbor(i)->info().getGridId());
		}
		if (neighborsId.find(cell->info().getGridId()) == neighborsId.end()) {
		/// the cell has no neighbors with the same material id
			++hangsCounter;
			GridId theMostCommonId = *neighborsId.begin();
			for (const GridId id : neighborsId) {
				if (neighborsId.count(id) > neighborsId.count(theMostCommonId)) {
					theMostCommonId = id;
				}
			}
			cell->info().setGridId(theMostCommonId);
		}
	}
	return hangsCounter;
}


template<int Dimensionality, typename VertexInfo, typename CellInfo>
bool
CgalTriangulation<Dimensionality, VertexInfo, CellInfo>::
clearFromDisconnectedCellSets(VertexHandle vh) {
	std::set<GridId> uniqueMaterials = incidentGridsIds(vh);
	if (uniqueMaterials.size() == 1) { return false; }
	CCSMultiset connectedCellsSets = constructCCSMultiset(vh);
	
	struct CcsInfo {
		GridId id;
		size_t numberOfCcsWithId = 0;
		size_t numberOfCellsWithId = 0;
	};
	std::list<CcsInfo> idsSortedByNumberOfCells;
	for (const GridId id : uniqueMaterials) {
		CcsInfo ccsInfo;
		ccsInfo.id = id;
		const auto range = connectedCellsSets.equal_range(ConnectedCellsSet(id));
		for (auto it = range.first; it != range.second; ++it) {
			++ccsInfo.numberOfCcsWithId;
			ccsInfo.numberOfCellsWithId += it->size();
		}
		idsSortedByNumberOfCells.push_back(ccsInfo);
	}
	idsSortedByNumberOfCells.sort([](const CcsInfo& a, const CcsInfo& b) {
		return a.numberOfCellsWithId > b.numberOfCellsWithId;}); // descending
	std::list<CcsInfo> idsSortedByNumberOfCcs(idsSortedByNumberOfCells);
	idsSortedByNumberOfCcs.sort([](const CcsInfo& a, const CcsInfo& b) {
		return a.numberOfCcsWithId < b.numberOfCcsWithId;}); // ascending
	// now, the id with the biggest number of disconnected cells sets (DCS)
	// and the lowest (among equal number of DCS) number of cells is at the end:
	CcsInfo toRemove = idsSortedByNumberOfCcs.back();
	if (toRemove.numberOfCcsWithId == 1) { return false; }
	// replace with id with the biggest number of cells:
//	CcsInfo toInsert = idsSortedByNumberOfCells.front();
//	if (toInsert.id == toRemove.id) {
//		toInsert = *std::next(idsSortedByNumberOfCells.begin());
//	}
	GridId idToInsert =
		Utils::chooseRandomElementExceptSpecified(uniqueMaterials, toRemove.id);
	
	removeDCS(toRemove.id, idToInsert, connectedCellsSets);
//	std::cout << coordsD(vh);
	return true;
}



template class CgalTriangulation<2, VertexInfo, CellInfoT<3>>;
template class CgalTriangulation<3, VertexInfo, CellInfoT<4>>;
