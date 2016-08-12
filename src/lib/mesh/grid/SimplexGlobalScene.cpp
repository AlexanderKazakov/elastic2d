#include <lib/mesh/grid/SimplexGlobalScene.hpp>
#include <lib/mesh/grid/cgal/CgalTriangulation.hpp>
#include <lib/Engine.hpp>


using namespace gcm;

template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
SimplexGlobalScene<Dimensionality, TriangulationT>::
SimplexGlobalScene(const Task& task, Engine* engine_) :
		engine(engine_),
		triangulation(task),
		movable(task.globalSettings.movable) { }


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void SimplexGlobalScene<Dimensionality, TriangulationT>::
afterGridsConstruction(const Task& task) {
	
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


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void SimplexGlobalScene<Dimensionality, TriangulationT>::
nextTimeStep() {
	
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


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void SimplexGlobalScene<Dimensionality, TriangulationT>::
createNewCalculationBasis() {
	calculationBasis = linal::randomBasis(calculationBasis);
	LOG_INFO("New calculation basis:" << calculationBasis);
	
	for (const auto& body : engine->bodies) {
		Grid* grid = dynamic_cast<Grid*>(engine->getAbstractMesh(body.first));
		assert_true(grid);
		grid->changeCalculationBasis(calculationBasis);
	}
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void SimplexGlobalScene<Dimensionality, TriangulationT>::
correctContactsAndBorders(const int stage) {
	
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


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void SimplexGlobalScene<Dimensionality, TriangulationT>::
createContacts(const Task& task) {
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


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void SimplexGlobalScene<Dimensionality, TriangulationT>::
createBorders(const Task& task) {
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


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void SimplexGlobalScene<Dimensionality, TriangulationT>::
addNode(const VertexHandle vh, const GridsPair gridsIds) {
	
	if (gridsIds.first == EmptySpaceFlag) {
		addBorderNode(vh, gridsIds.second);
	} else if (gridsIds.second == EmptySpaceFlag) {
		addBorderNode(vh, gridsIds.first);
	} else {
		addContactNode(vh, gridsIds);
	}
	
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void SimplexGlobalScene<Dimensionality, TriangulationT>::
addContactNode(const VertexHandle vh, const GridsPair gridsIds) {
	
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


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void SimplexGlobalScene<Dimensionality, TriangulationT>::
addBorderNode(const VertexHandle vh, const GridId gridId) {
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



template class SimplexGlobalScene<2, CgalTriangulation>;
template class SimplexGlobalScene<3, CgalTriangulation>;

