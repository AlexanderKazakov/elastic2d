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
	
//	for (auto it  = this->allCellsBegin(); it != this->allCellsEnd(); ++it) {
//		if (isInfinite(it)) { continue; }
//		GridId id = it->info().getGridId();
//		if (id != 1) {
//			it->info().setGridId(4);
//		}
//	}
	
	int disconnectedCellsSetsCaseCounter = 1;
	int iterationsCounter = 0;
	while (disconnectedCellsSetsCaseCounter > 0 && ++iterationsCounter < 4) {
		disconnectedCellsSetsCaseCounter = 0;
		LOG_INFO("Start " << iterationsCounter << "'th iteration of "
				 << "cleaning from disconnected cells sets");
		for (auto vertexIter = verticesBegin();
				  vertexIter != verticesEnd(); ++vertexIter) {
			disconnectedCellsSetsCaseCounter += int(
					clearFromDisconnectedCellSets(vertexIter));
		}
		LOG_INFO(disconnectedCellsSetsCaseCounter << " bad cases was found");
	}
}


template<int Dimensionality, typename VertexInfo, typename CellInfo>
bool
CgalTriangulation<Dimensionality, VertexInfo, CellInfo>::
clearFromDisconnectedCellSets(VertexHandle vh) {
	
	struct ConnectedCellsSet {
		ConnectedCellsSet(const GridId gridId) : centerVertex(NULL), id(gridId) { }
		
		ConnectedCellsSet(VertexHandle v, CellHandle c) :
				centerVertex(v), id(c->info().getGridId()) {
			tryToInsert(c);
		}
		
		void tryToInsert(CellHandle c) {
			if (c->has_vertex(centerVertex) && c->info().getGridId() == id) {
				if (set.insert(c).second) {
					for (int i = 0; i < CELL_POINTS_NUMBER; i++) {
						tryToInsert(c->neighbor(i));
					}
				}
			}
		}
		
		bool operator<(const ConnectedCellsSet& other) const {
			return this->id < other.id;
		}
		
		size_t size() const { return set.size(); }
		
		const VertexHandle centerVertex;
		const GridId id;
		std::set<CellHandle> set;
	};
	
	
	const std::list<CellHandle> listAllCells = this->allIncidentCells(vh);
	std::set<GridId> uniqueMaterials;
	for (const CellHandle c : listAllCells) {
		if (isInfinite(c)) { return false; }
		uniqueMaterials.insert(c->info().getGridId());
	}
	if (uniqueMaterials.size() != 2) { return false; }
	
	std::set<CellHandle> allCells(listAllCells.begin(), listAllCells.end());
	typedef std::multiset<ConnectedCellsSet> Multiset;
	typedef typename Multiset::const_iterator Iter;
	Multiset connectedCellsSets;
	while (!allCells.empty()) {
		ConnectedCellsSet newSet(vh, *allCells.begin());
		connectedCellsSets.insert(newSet);
		std::set<CellHandle> difference; // == allCells \ newSet
		std::set_difference(allCells.begin(), allCells.end(),
				newSet.set.begin(), newSet.set.end(),
				std::inserter(difference, difference.begin()));
		allCells = difference;
	}
	
	ConnectedCellsSet m1Key(*uniqueMaterials.rbegin());
	const bool disconnectedSetsFoundInM1 = connectedCellsSets.count(m1Key) > 1;
	ConnectedCellsSet m2Key(*uniqueMaterials.begin());
	const bool disconnectedSetsFoundInM2 = connectedCellsSets.count(m2Key) > 1;
	const bool disconnectedSetsFound =
			disconnectedSetsFoundInM1 || disconnectedSetsFoundInM2;
	
	auto removeDSC = [&](const ConnectedCellsSet idToRemove, const GridId idToInsert) {
		const std::pair<Iter, Iter> range = connectedCellsSets.equal_range(idToRemove);
		Iter largestSet = range.first;
		for (auto it = range.first; it != range.second; ++it) {
			if (largestSet->size() < it->size()) {
				largestSet = it;
			}
		}
		for (auto it = range.first; it != range.second; ++it) {
			if (it != largestSet) {
				for (CellHandle c : it->set) if (!isInfinite(c)) {
					c->info().setGridId(idToInsert);
				}
			}
		}
	};
	
	if (disconnectedSetsFoundInM1) {
		removeDSC(m1Key, m2Key.id);
	} else if (disconnectedSetsFoundInM2) {
		removeDSC(m2Key, m1Key.id);
	}
	
	return disconnectedSetsFound;
}



template class CgalTriangulation<2, VertexInfo, CellInfoT<3>>;
template class CgalTriangulation<3, VertexInfo, CellInfoT<4>>;
