#ifndef LIBGCM_SIMPLEXGLOBALSCENE_HPP
#define LIBGCM_SIMPLEXGLOBALSCENE_HPP

#include <lib/mesh/grid/AbstractGlobalScene.hpp>
#include <lib/mesh/grid/SimplexGrid.hpp>
#include <lib/numeric/gcm/ContactCorrector.hpp>
#include <lib/Engine.hpp>


namespace gcm {

/**
 * Global scene of the program -- a conglomeration of several grids --
 * for the case of calculation on simplex grids
 */
template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
class SimplexGlobalScene : public AbstractGlobalScene {
public:
	
	typedef SimplexGrid<Dimensionality, TriangulationT>    Grid;
	typedef typename Grid::Triangulation                   Triangulation;
	typedef typename Grid::Iterator                        Iterator;
	typedef typename Grid::GridId                          GridId;
	typedef typename Grid::RealD                           RealD;
	
	typedef typename Triangulation::VertexHandle           VertexHandle;
	typedef typename Triangulation::CellHandle             CellHandle;
	
	static const GridId EmptySpaceFlag = Grid::EmptySpaceFlag;
	
	typedef AbstractContactCorrector<Grid>                 ContactCorrector;
	typedef typename ContactCorrector::NodesContact        NodesContact;
	
	
	/// pair of grids in contact
	typedef std::pair<GridId, GridId>                      GridsPair;
	
	struct Contact {
		/// list of nodes pairs in contact
		std::list<NodesContact> nodesInContact;
		/// corrector to handle this contacts
		std::shared_ptr<ContactCorrector> contactCorrector;
	};
	
	
	SimplexGlobalScene(const Task& task, Engine* engine_ = nullptr) :
			engine(engine_),
			triangulation(task),
			movable(task.globalSettings.movable) { }
	virtual ~SimplexGlobalScene() { }
	
	
	/**
	 * Actions to perform after all grids would be constructed --
	 * find and remember all contacts.
	 */
	virtual void afterGridsConstruction(const Task& task) override {
		
		createContacts(task);
		
		for (auto vertexIter  = triangulation.verticesBegin();
		          vertexIter != triangulation.verticesEnd(); ++vertexIter) {
			
			std::set<GridId> incidentGrids = incidentGridsIds(vertexIter);
			incidentGrids.erase((GridId)EmptySpaceFlag); //< FIXME - move border corrector here too
			
			const auto gridPairs = Utils::makePairs(incidentGrids);
			if (gridPairs.size() != 1) { continue; }
			addNodesContact(vertexIter, gridPairs.front());
			
			// FIXME - now, only one pair is handled
//			for (const GridsPair gridsPair : Utils::makePairs(incidentGrids)) {
//				addNodesContact(vertexIter, gridsPair);
//			}
		}
		
		LOG_INFO("Found contacts:");
		for (const auto& contact : contacts) {
			LOG_INFO("For bodies " << contact.first.first << " and "
					<< contact.first.second << " number of contact nodes = "
					<< contact.second.nodesInContact.size());
		}
	}
	
	
	/**
	 * Apply contact correctors
	 */
	virtual void correctContacts() override {
		for (const auto& contact : contacts) {
			contact.second.contactCorrector->apply(
					engine->getAbstractMesh(contact.first.first),
					engine->getAbstractMesh(contact.first.second),
					contact.second.nodesInContact);
		}
	}
	
	
	
private:
	Engine * const engine;
	
	/// global triangulation of the whole calculation space
	Triangulation triangulation;
	
	/// all contacts of all grids
	std::map<GridsPair, Contact> contacts;
	
	/// on/off points motion
	bool movable = false;
	
	
	friend class SimplexGrid<Dimensionality, TriangulationT>;
	USE_AND_INIT_LOGGER("gcm.SimplexGlobalScene")
	
	
	/** Create object in contacts for each grid-grid pair */
	void createContacts(const Task& task) {
		assert_true(engine);
		
		std::set<GridId> gridsIds;
		for (const auto& body : engine->bodies) {
			gridsIds.insert(body.first);
		}
		
		for (const GridsPair gridsPair : Utils::makePairs(gridsIds)) {
			
			ContactConditions::T condition = task.contactCondition.defaultCondition;
			auto conditionIter = task.contactCondition.gridToGridConditions.find(gridsPair);
			if (conditionIter != task.contactCondition.gridToGridConditions.end()) {
				condition = conditionIter->second;
			}
			
			AbstractGrid* firstGrid = engine->getAbstractMesh(gridsPair.first);
			AbstractGrid* secondGrid = engine->getAbstractMesh(gridsPair.second);
			Contact contact;
			contact.contactCorrector = ContactCorrectorFactory<Grid>::create(
					condition, firstGrid->modelId, firstGrid->materialId,
							secondGrid->modelId, secondGrid->materialId);
			
			contacts.insert({gridsPair, contact});
		}
		
		LOG_INFO("There are " << gridsIds.size() << " bodies. "
				<< contacts.size() << " body-to-body contacts are possible");
	}
	
	
	std::set<GridId> incidentGridsIds(const VertexHandle vh) const {
		std::list<CellHandle> incidentCells = triangulation.allIncidentCells(vh);
		std::set<GridId> ans;
		for (CellHandle ch : incidentCells) {
			ans.insert(ch->info().getGridId());
		}
		return ans;
	}
	
	
	void addNodesContact(const VertexHandle vh, const GridsPair gridsIds) {
		
		Grid* firstGrid = dynamic_cast<Grid*>(engine->getAbstractMesh(gridsIds.first));
		assert_true(firstGrid);
		Grid* secondGrid = dynamic_cast<Grid*>(engine->getAbstractMesh(gridsIds.second));
		assert_true(secondGrid);
		
		Iterator first = firstGrid->localVertexIndex(vh);
		Iterator second = secondGrid->localVertexIndex(vh);
		RealD normal = firstGrid->contactNormal(first, gridsIds.second);
		if (normal != RealD::Zeros()) {
			contacts.at(gridsIds).nodesInContact.push_back(
					{first, second, normal});
		}
	}
	
	
};


}

#endif // LIBGCM_SIMPLEXGLOBALSCENE_HPP
