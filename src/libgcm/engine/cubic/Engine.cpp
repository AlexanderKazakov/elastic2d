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
		body.mesh->setUpPde(task);
		body.gcm = body.factory->createGcm(task);
		body.border = body.factory->createBorder(task, body.mesh);
		
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
		body.mesh = body.factory->createMesh(task, gridId,
				createGridConstructionPack(task.cubicGrid, gridId), 1);
		bodies.push_back(body);
	}
	
	for (Body& body : bodies) {
		for (const Body& other : bodies) if (other != body) {
			
			AABB intersection = AABB::intersection(
					body.mesh->aabb(), other.mesh->aabb());
			if (intersection.valid()) {
				THROW_BAD_MESH("Bodies must not intersect");
			}
			if (intersection.minimalWidth().first != -1) {
				continue; // no contact
			}
			
			typename Body::Contact contact;
			contact.neighborId = other.mesh->id;
			contact.direction = intersection.minimalWidth().second;
			
			// create box for contact copying
			AABB buffer = intersection;
			if (body.mesh->start(contact.direction) >
			   other.mesh->start(contact.direction)) {
			// copy from the bottom
				buffer.min(contact.direction) -= body.mesh->borderSize;
			} else {
			// copy from the top
				buffer.max(contact.direction) += body.mesh->borderSize;
			}
			contact.copier = body.factory->createContact(
					 body.mesh->box( body.mesh->globalToLocal(buffer)),
					other.mesh->box(other.mesh->globalToLocal(buffer)),
					ContactConditions::T::ADHESION,
					other.mesh->getModelType(),
					other.mesh->getMaterialType());
			
			body.contacts.push_back(contact);
		}
	}
}


template<int Dimensionality>
void Engine<Dimensionality>::
nextTimeStep() {
	for (int stage = 0; stage < Dimensionality; stage++) {
		
		for (Body& body : bodies) {
			body.border->apply(*body.mesh, stage);
		}
		
		for (Body& body : bodies) {
			for (typename Body::Contact& contact : body.contacts) {
				if (contact.direction == stage) {
					contact.copier->apply(
							*body.mesh,
							*getBody(contact.neighborId).mesh);
				}
			}
		}
		
		for (Body& body : bodies) {
			body.gcm->stage(stage, Clock::TimeStep(), *body.mesh);
			body.mesh->swapCurrAndNextPdeTimeLayer(0);
		}
		
	}
	
	for (Body& body : bodies) {
		for (typename Body::OdePtr ode : body.odes) {
			ode->apply(*body.mesh, Clock::TimeStep());
		}
	}
}


template<int Dimensionality>
real Engine<Dimensionality>::
estimateTimeStep() {
	real maxEigenvalue = 0;
	RealD h = bodies.front().mesh->h;
	for (const Body& body : bodies) {
		assert_true(h == body.mesh->h);
		real eigenvalue = body.mesh->getMaximalEigenvalue();
		if (eigenvalue > maxEigenvalue) {
			maxEigenvalue = eigenvalue;
		}
	}
	
	return CourantNumber * 
			bodies.front().mesh->getMinimalSpatialStep() /
			maxEigenvalue;
}


template<int Dimensionality>
void Engine<Dimensionality>::
writeSnapshots(const int step) {
	for (Body& body : bodies) {
		for (typename Body::SnapPtr snapshotter : body.snapshotters) {
			snapshotter->snapshot(body.mesh.get(), step);
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

