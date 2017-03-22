#include <libgcm/engine/simplex/Engine.hpp>
#include <libgcm/grid/simplex/cgal/CgalTriangulation.hpp>

#include <limits>

using namespace gcm;
using namespace gcm::simplex;


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
Engine<Dimensionality, TriangulationT>::
Engine(const Task& task) :
		AbstractEngine(task),
		triangulation(task),
		movable(task.simplexGrid.movable),
		borderCalcMode(task.simplexGrid.borderCalcMode),
		gcmType(task.globalSettings.gcmType) {
	
	initializeCalculationBasis(task);
	createMeshes(task);
	createContacts(task);
	
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
	for (const Body& body : bodies) {
		for (size_t i = 0; i < body.borders.size(); i++) {
			LOG_INFO("For body " << body.mesh->id
					<< " and border condition number " << i
					<< " number of border nodes = "
					<< body.borders[i].borderNodes.size());
		}
	}
	
	afterConstruction(task);
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void Engine<Dimensionality, TriangulationT>::
createMeshes(const Task& task) {
	for (const auto& taskBody : task.bodies) {
		Body body;
		const GridId gridId = taskBody.first;
		auto factory = createAbstractFactory(taskBody.second);
		
		body.mesh = factory->createMesh(task, gridId, {&triangulation}, Dimensionality);
		body.mesh->setUpPde(task, calculationBasis.basis, borderCalcMode);
		
		body.gcm = factory->createGcm(gcmType);
		
		for (const Snapshotters::T snapType : task.globalSettings.snapshottersId) {
			body.snapshotters.push_back(
					factory->createSnapshotter(task, snapType));
		}
		
		for (const Odes::T odeType : taskBody.second.odes) {
			body.odes.push_back(factory->createOde(odeType));
		}
		
		for (const Task::BorderCondition& condition : task.borderConditions) {
			Border border;
			border.correctionArea = condition.area;	
			border.borderCorrector = BorderCorrectorFactory<Grid>::create(
					gcmType,
					condition,
					task.bodies.at(body.mesh->id).modelId,
					task.bodies.at(body.mesh->id).materialId);
			body.borders.push_back(border);
		}
		
		bodies.push_back(body);
	}
	
	assert_eq(task.bodies.size(), bodies.size());
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void Engine<Dimensionality, TriangulationT>::
nextTimeStep() {
	changeCalculationBasis();
	
	/// simple first order splitting by summ solutions from all directions
	for (int stage = 0; stage < Dimensionality; stage++) {
		gcmStage(stage, Clock::Time(), Clock::TimeStep());
	}
	for (const Body& body : bodies) {
		body.mesh->averageNewPdeLayersToCurrent();
	}
	
	for (const Body& body : bodies) {
		for (typename Body::OdePtr ode : body.odes) {
			ode->apply(*body.mesh, Clock::TimeStep());
		}
	}
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void Engine<Dimensionality, TriangulationT>::
gcmStage(const int stage, const real currentTime, const real timeStep) {
	for (const Body& body : bodies) {
		body.gcm->beforeStage(stage, *body.mesh);
	}
	for (const Body& body : bodies) {
		body.gcm->contactAndBorderStage(stage, timeStep, *body.mesh);
	}
	correctContactsAndBorders(stage, currentTime + timeStep);
	for (const Body& body : bodies) {
		body.gcm->innerStage(stage, timeStep, *body.mesh);
	}
	for (const Body& body : bodies) {
		body.gcm->afterStage(stage, *body.mesh);
	}
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void Engine<Dimensionality, TriangulationT>::
correctContactsAndBorders(const int stage, const real timeAtNextLayer) {
	switch (borderCalcMode) {
		case BorderCalcMode::GLOBAL_BASIS:
			for (const auto& contact : contacts) {
				contact.second.contactCorrector->applyInGlobalBasis(
						stage,
						getBody(contact.first.first).mesh,
						getBody(contact.first.second).mesh,
						contact.second.nodesInContact);
			}
			for (const Body& body : bodies) {
				for (const Border& border : body.borders) {
					border.borderCorrector->applyInGlobalBasis(
							stage,
							body.mesh,
							border.borderNodes,
							timeAtNextLayer);
				}
			}
			break;
			
		case BorderCalcMode::LOCAL_BASIS:
			if (stage != 0) { return; }
			for (const auto& contact : contacts) {
				contact.second.contactCorrector->applyInLocalBasis(
						getBody(contact.first.first).mesh,
						getBody(contact.first.second).mesh,
						contact.second.nodesInContact);
			}
			for (const Body& body : bodies) {
				for (const Border& border : body.borders) {
					border.borderCorrector->applyInLocalBasis(
							body.mesh,
							border.borderNodes,
							timeAtNextLayer);
				}
			}
			break;
			
		default:
			THROW_BAD_CONFIG("Unknown border calculation mode");
	}
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void Engine<Dimensionality, TriangulationT>::
createContacts(const Task& task) {
	
	std::set<GridId> gridsIds;
	for (const Body& body : bodies) {
		gridsIds.insert(body.mesh->id);
	}
	
	for (const GridsPair gridsPair : Utils::makePairs(gridsIds)) {
		
		ContactConditions::T condition = task.contactCondition.defaultCondition;
		auto conditionIter = task.contactCondition.gridToGridConditions.find(gridsPair);
		if (conditionIter != task.contactCondition.gridToGridConditions.end()) {
			condition = conditionIter->second;
		}
		
		Contact contact;
		contact.contactCorrector = ContactCorrectorFactory<Grid>::create(
				condition,
				task.bodies.at(gridsPair.first).modelId,
				task.bodies.at(gridsPair.first).materialId,
				task.bodies.at(gridsPair.second).modelId,
				task.bodies.at(gridsPair.second).materialId);
		
		contacts.insert({ gridsPair, contact });
	}
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void Engine<Dimensionality, TriangulationT>::
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
void Engine<Dimensionality, TriangulationT>::
addContactNode(const VertexHandle vh, const GridsPair gridsIds) {
	
	Iterator firstIter = getBody(gridsIds.first).mesh->localVertexIndex(vh);
	Iterator secondIter = getBody(gridsIds.second).mesh->localVertexIndex(vh);
	RealD normal = getBody(gridsIds.first).mesh->contactNormal(firstIter, gridsIds.second);
	if (normal != RealD::Zeros()) {
		contacts.at(gridsIds).nodesInContact.push_back(
				{ firstIter, secondIter, normal });
	}
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void Engine<Dimensionality, TriangulationT>::
addBorderNode(const VertexHandle vh, const GridId gridId) {
	assert_ne(gridId, EmptySpaceFlag);
	
	std::shared_ptr<Mesh> mesh = getBody(gridId).mesh;
	Iterator iter = mesh->localVertexIndex(vh);
	
	// for a concrete node, not more than one border condition can be applied
	Border* chosenBorder = nullptr;
	for (Border& border : getBody(gridId).borders) {
		if (border.correctionArea->contains(mesh->coords(iter))) {
			chosenBorder = &border;
		}
	}
	if (chosenBorder == nullptr) { return; }
	
	RealD normal = mesh->borderNormal(iter);
	if (normal != RealD::Zeros()) {
		chosenBorder->borderNodes.push_back({ iter, normal });
	}
	
}


template<int Dimensionality,
         template<int, typename, typename> class TriangulationT>
void Engine<Dimensionality, TriangulationT>::
writeSnapshots(const int step) {
	for (Body& body : bodies) {
		for (typename Body::SnapPtr snapshotter : body.snapshotters) {
			snapshotter->snapshot(body.mesh.get(), step);
		}
	}
}



template class Engine<2, CgalTriangulation>;
template class Engine<3, CgalTriangulation>;

