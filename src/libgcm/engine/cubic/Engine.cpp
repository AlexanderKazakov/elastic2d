#include <libgcm/engine/cubic/Engine.hpp>

#include <libgcm/engine/cubic/AbstractFactory.hpp>
#include <libgcm/engine/mesh/DefaultMesh.hpp>
#include <libgcm/rheology/models/AcousticModel.hpp>
#include <libgcm/rheology/models/ElasticModel.hpp>


using namespace gcm;
using namespace gcm::cubic;


template<int Dimensionality>
Engine<Dimensionality>::Engine(const Task& task) :
		AbstractEngine(task) {
	
	createGridsAndContacts(task);
	
	for (const auto& taskBody : task.bodies) {
		Body& body = getBody(taskBody.first);
		body.grid->setUpPde(task);
		body.gcm = body.factory->createGcm(task);
		
		for (const Snapshotters::T snapType : task.globalSettings.snapshottersId) {
			body.snapshotters.push_back(
					body.factory->createSnapshotter(task, snapType));
		}
		
		for (const Odes::T odeType : taskBody.second.odes) {
			body.odes.push_back(body.factory->createOde(odeType));
		}
	}
	
	afterConstruction(task);
}


template<int Dimensionality>
void Engine<Dimensionality>::
createGridsAndContacts(const Task& task) {
	assert_eq(task.bodies.size(), task.cubicGrid.cubics.size());
	for (const auto& taskBody: task.bodies) {
		Body body;
		const GridId gridId = taskBody.first;
		body.factory = createAbstractFactory(taskBody.second);
		body.grid = body.factory->createMesh(task, gridId,
				createGridConstructionPack(task.cubicGrid, gridId));
		bodies.push_back(body);
	}
	
	for (Body& body : bodies) {
		for (const Body& other : bodies) if (other != body) {
			AABB intersection = AABB::intersection(
					body.grid->aabb(), other.grid->aabb());
			if ( !intersection.valid() ) { continue; } // no intersection
			
			typename Body::Contact contact;
			contact.neighborId = other.grid->id;
			contact.direction = intersection.sliceDirection();
			
			// create box for contact copying
			AABB buffer = intersection;
			if (body.grid->start(contact.direction) >
			   other.grid->start(contact.direction)) {
			// copy from the bottom
				buffer.min(contact.direction) -= body.grid->borderSize;
				buffer.max(contact.direction) -= 1;
			} else {
			// copy from the top
				buffer.min(contact.direction) += 1;
				buffer.max(contact.direction) += body.grid->borderSize;
			}
			contact.copier = body.factory->createContact(
					 body.grid->box( body.grid->globalToLocal(buffer)),
					other.grid->box(other.grid->globalToLocal(buffer)),
					ContactConditions::T::ADHESION,
					other.grid->getModelType(),
					other.grid->getMaterialType());
			
			body.contacts.push_back(contact);
		}
	}
}


template<int Dimensionality>
void Engine<Dimensionality>::
nextTimeStep() {
	for (int stage = 0; stage < Dimensionality; stage++) {
		
		for (Body& body : bodies) {
			for (typename Body::Contact& contact : body.contacts) {
				if (contact.direction == stage) {
					contact.copier->apply(
							*body.grid,
							*getBody(contact.neighborId).grid);
				}
			}
		}
		
		for (Body& body : bodies) {
			body.gcm->stage(stage, Clock::TimeStep(), *body.grid);
			body.grid->swapPdeTimeLayers();
		}
		
	}
	
	for (Body& body : bodies) {
		for (typename Body::OdePtr ode : body.odes) {
			ode->apply(*body.grid, Clock::TimeStep());
		}
	}
}


template<int Dimensionality>
real Engine<Dimensionality>::
estimateTimeStep() {
	/// minimal among all bodies
	real minimalTimeStep = std::numeric_limits<real>::max();
	for (const Body& body : bodies) {
		real bodyTimeStep = body.gcm->calculateTimeStep(*body.grid, CourantNumber);
		if (bodyTimeStep < minimalTimeStep) {
			minimalTimeStep = bodyTimeStep;
		}
	}
	return minimalTimeStep;
}


template<int Dimensionality>
void Engine<Dimensionality>::
writeSnapshots(const int step) {
	for (Body& body : bodies) {
		for (typename Body::SnapPtr snapshotter : body.snapshotters) {
			snapshotter->snapshot(body.grid.get(), step);
		}
	}
}


template<int Dimensionality>
std::shared_ptr<AbstractFactoryBase<typename Engine<Dimensionality>::Grid>>
Engine<Dimensionality>::
createAbstractFactory(const Task::Body& body) {
	switch (body.materialId) {
	
	case Materials::T::ISOTROPIC:
	switch (body.modelId) {
		case (Models::T::ACOUSTIC):
			return std::make_shared<AbstractFactory<
					AcousticModel<Dimensionality>,
					Grid, IsotropicMaterial, DefaultMesh>>();
		case (Models::T::ELASTIC):
			return std::make_shared<AbstractFactory<
					ElasticModel<Dimensionality>,
					Grid, IsotropicMaterial, DefaultMesh>>();
		default:
			THROW_UNSUPPORTED("Unknown model type");
	}
	
	case Materials::T::ORTHOTROPIC:
	switch (body.modelId) {
		case (Models::T::ELASTIC):
			return std::make_shared<AbstractFactory<
					ElasticModel<Dimensionality>,
					Grid, OrthotropicMaterial, DefaultMesh>>();
		default:
			THROW_UNSUPPORTED("Unknown or inappropriate model type");
	}
	
	default:
	THROW_UNSUPPORTED("Unknown material type");
	
	}
}


template<int Dimensionality>
typename Engine<Dimensionality>::GridConstructionPack
Engine<Dimensionality>::
createGridConstructionPack(const Task::CubicGrid& task, const GridId gridId) {
	
	GridConstructionPack pack;
	pack.borderSize = task.borderSize;
	assert_eq(task.h.size(), DIMENSIONALITY);
	pack.h.copyFrom(task.h);
	
	Task::CubicGrid::Cube cube = task.cubics.at(gridId);
	assert_eq(cube.sizes.size(), DIMENSIONALITY);
	pack.sizes.copyFrom(cube.sizes);
	assert_eq(cube.start.size(), DIMENSIONALITY);
	pack.start.copyFrom(cube.start);
	
	return pack;
}



template class Engine<1>;
template class Engine<2>;
template class Engine<3>;

