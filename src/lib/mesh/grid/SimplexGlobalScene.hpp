#ifndef LIBGCM_SIMPLEXGLOBALSCENE_HPP
#define LIBGCM_SIMPLEXGLOBALSCENE_HPP

#include <lib/mesh/grid/AbstractGlobalScene.hpp>
#include <lib/mesh/grid/SimplexGrid.hpp>
#include <lib/numeric/gcm/ContactCorrector.hpp>
#include <lib/numeric/gcm/BorderCorrector.hpp>
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
	typedef typename Grid::MatrixDD                        MatrixDD;
	
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
		/// corrector to handle these contacts
		std::shared_ptr<ContactCorrector> contactCorrector;
	};
	
	
	typedef AbstractBorderCorrector<Grid>                  BorderCorrector;
	typedef typename BorderCorrector::NodeBorder           NodeBorder;
	
	struct Border {
		/// list of border nodes
		std::list<NodeBorder> borderNodes;
		/// corrector to handle these nodes
		std::shared_ptr<BorderCorrector> borderCorrector;
		/// used only at instantiation step, not in calculations
		std::shared_ptr<Area> correctionArea;
	};
	
	
	SimplexGlobalScene(const Task& task, Engine* engine_ = nullptr) :
			engine(engine_),
			triangulation(task),
			movable(task.globalSettings.movable) { }
	virtual ~SimplexGlobalScene() { }
	
	
	/**
	 * Actions to perform after all grids would be constructed --
	 * find and remember all contacts and borders.
	 */
	virtual void afterGridsConstruction(const Task& task) override {
		
		createContacts(task);
		createBorders(task);
		
		for (auto vertexIter  = triangulation.verticesBegin();
		          vertexIter != triangulation.verticesEnd(); ++vertexIter) {
			
			std::set<GridId> incidentGrids = incidentGridsIds(vertexIter);
			const auto gridPairs = Utils::makePairs(incidentGrids);
			if (gridPairs.size() != 1) { continue; } // TODO
			addNode(vertexIter, gridPairs.front());
			
		}
		
		LOG_INFO("Found contacts:");
		for (const auto& contact : contacts) {
			LOG_INFO("For bodies " << contact.first.first << " and "
					<< contact.first.second << " number of contact nodes = "
					<< contact.second.nodesInContact.size());
		}
		
		LOG_INFO("Found borders (except non-reflection cases):");
		for (const auto& bodyBorders : borders) {
			for (size_t i = 0; i < bodyBorders.second.size(); i++) {
				LOG_INFO("For body " << bodyBorders.first
						<< " and border condition number " << i
						<< " number of border nodes = "
						<< bodyBorders.second[i].borderNodes.size());
			}
		}
		
	}
	
	
	virtual void nextTimeStep() override {
		
		createNewCalculationBasis();
		
		for (int stage = 0; stage < Dimensionality; stage++) {
			
			for (const auto body : engine->bodies) {
				body.second.solver->beforeStage();
			}
			
			for (const auto body : engine->bodies) {
				body.second.solver->contactStage(stage, Clock::TimeStep());
			}
			
			correctContactsAndBorders(stage);
			
			for (const auto body : engine->bodies) {
				body.second.solver->privateStage(stage, Clock::TimeStep());
			}
			
		}
		
	}
	
	
private:
	Engine * const engine;
	
	/// global triangulation of the whole calculation space
	Triangulation triangulation;
	
	/// on/off points motion
	bool movable = false;
	
	
	/// Current basis of calculations
	MatrixDD calculationBasis = linal::randomBasis(MatrixDD());
	
	
	/// all contact conditions of all grids
	std::map<GridsPair, Contact> contacts;
	
	/// it's possible to have several border conditions for one grid
	/// all border conditions of all grids
	std::map<GridId, std::vector<Border>> borders;
	
	
	friend class SimplexGrid<Dimensionality, TriangulationT>;
	USE_AND_INIT_LOGGER("gcm.SimplexGlobalScene")
	
	
	void createNewCalculationBasis() {
		calculationBasis = linal::randomBasis(calculationBasis);
		LOG_INFO("New calculation basis:" << calculationBasis);
		
		for (const auto& body : engine->bodies) {
			Grid* grid = dynamic_cast<Grid*>(engine->getAbstractMesh(body.first));
			assert_true(grid);
			grid->changeCalculationBasis(calculationBasis);
		}
	}
	
	
	void correctContactsAndBorders(const int stage) {
		
		for (const auto& contact : contacts) {
			contact.second.contactCorrector->apply(
					engine->getAbstractMesh(contact.first.first),
					engine->getAbstractMesh(contact.first.second),
					contact.second.nodesInContact,
					calculationBasis.getColumn(stage));
		}
		
		for (const auto& bodyWithBorders : borders) {
			for (const Border& border : bodyWithBorders.second) {
				border.borderCorrector->apply(
						engine->getAbstractMesh(bodyWithBorders.first),
						border.borderNodes,
						calculationBasis.getColumn(stage));
			}
		}
		
	}
	
	
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
			
			contacts.insert({ gridsPair, contact });
		}
	}
	
	
	/** Create object in contacts for each grid border */
	void createBorders(const Task& task) {
		assert_true(engine);
		
		assert_eq(task.statements.size(), 1); // FIXME - remove statements at all
		Statement statement = task.statements.front();
		
		for (const auto& body : engine->bodies) {
			std::vector<Border> bodyBorderConditions;
			for (const Statement::BorderCondition condition :
					statement.borderConditions) {
				AbstractGrid* grid = engine->getAbstractMesh(body.first);
				Border border;
				border.borderCorrector = BorderCorrectorFactory<Grid>::create(
						condition, grid->modelId, grid->materialId);
				border.correctionArea = condition.area;
				
				bodyBorderConditions.push_back(border);
			}
			borders.insert({ body.first, bodyBorderConditions });
		}
	}
	
	
	void addNode(const VertexHandle vh, const GridsPair gridsIds) {
		
		if (gridsIds.first == EmptySpaceFlag) {
			addBorderNode(vh, gridsIds.second);
		} else if (gridsIds.second == EmptySpaceFlag) {
			addBorderNode(vh, gridsIds.first);
		} else {
			addContactNode(vh, gridsIds);
		}
		
	}
	
	
	void addContactNode(const VertexHandle vh, const GridsPair gridsIds) {
		
		Grid* firstGrid = dynamic_cast<Grid*>(engine->getAbstractMesh(gridsIds.first));
		assert_true(firstGrid);
		Grid* secondGrid = dynamic_cast<Grid*>(engine->getAbstractMesh(gridsIds.second));
		assert_true(secondGrid);
		
		Iterator firstIter = firstGrid->localVertexIndex(vh);
		Iterator secondIter = secondGrid->localVertexIndex(vh);
		RealD normal = firstGrid->contactNormal(firstIter, gridsIds.second);
		if (normal != RealD::Zeros()) {
			contacts.at(gridsIds).nodesInContact.push_back(
					{ firstIter, secondIter, normal });
		}
		
	}
	
	
	void addBorderNode(const VertexHandle vh, const GridId gridId) {
		assert_ne(gridId, EmptySpaceFlag);
		Grid* grid = dynamic_cast<Grid*>(engine->getAbstractMesh(gridId));
		assert_true(grid);
		
		Iterator iter = grid->localVertexIndex(vh);
		
		// for a concrete node, not more than one border condition can be applied
		Border* chosenBorder = nullptr;
		for (Border& border : borders.at(gridId)) {
			if (border.correctionArea->contains(grid->coords(iter))) {
				chosenBorder = &border;
			}
		}
		if (chosenBorder == nullptr) { return; }
		
		RealD normal = grid->borderNormal(iter);
		if (normal != RealD::Zeros()) {
			chosenBorder->borderNodes.push_back({ iter, normal });
		}
		
	}
	
	
	std::set<GridId> incidentGridsIds(const VertexHandle vh) const {
		std::list<CellHandle> incidentCells = triangulation.allIncidentCells(vh);
		std::set<GridId> ans;
		for (CellHandle ch : incidentCells) {
			ans.insert(ch->info().getGridId());
		}
		return ans;
	}
	
	
};


}

#endif // LIBGCM_SIMPLEXGLOBALSCENE_HPP
