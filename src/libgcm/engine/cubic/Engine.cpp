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
	
	assert_eq(1, task.bodies.size()); // TODO
	
	for (const auto& taskBody : task.bodies) {
		Body body;
		const GridId gridId = taskBody.first;
		auto factory = createAbstractFactory(taskBody.second);
		
		body.grid = factory->createMesh(task, gridId,
				createGridConstructionPack(task, gridId));
		
		body.gcm = factory->createGcm(task);
		
		for (const Snapshotters::T snapType : task.globalSettings.snapshottersId) {
			body.snapshotters.push_back(
					factory->createSnapshotter(task, snapType));
		}
		
		for (const Odes::T odeType : taskBody.second.odes) {
			body.odes.push_back(factory->createOde(odeType));
		}
		
		bodies.push_back(body);
	}
	
	afterGridsConstruction(task);
}


template<int Dimensionality>
void Engine<Dimensionality>::
nextTimeStep() {
	
	for (int stage = 0; stage < Dimensionality; stage++) {
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
createGridConstructionPack(const Task& task, const GridId) {
	assert_eq(1, task.bodies.size()); // TODO
	
	GridConstructionPack pack;
	pack.borderSize = task.cubicGrid.borderSize;
	pack.sizes = calculateSizes(task.cubicGrid);
	pack.h = calculateH(task.cubicGrid);
	pack.startR = calculateStartR(task.cubicGrid);
	
	return pack;
}


template<int Dimensionality>
typename Engine<Dimensionality>::IntD
Engine<Dimensionality>::
calculateSizes(const Task::CubicGrid& task) {
	
	IntD _sizes;
	if (!task.h.empty() && !task.lengths.empty()) {
		assert_true(task.sizes.empty());
		assert_eq(task.h.size(), DIMENSIONALITY);
		assert_eq(task.lengths.size(), DIMENSIONALITY);
		
		RealD taskLength;
		taskLength.copyFrom(task.lengths);
		RealD taskH;
		taskH.copyFrom(task.h);
	
		_sizes = linal::plainDivision(taskLength, taskH) + IntD::Ones();
	
	} else {
		assert_eq(task.sizes.size(), DIMENSIONALITY);
		_sizes.copyFrom(task.sizes);
	}
	
	if (Mpi::ForceSequence()) {
		return _sizes;
	}
	
	// MPI - we divide the grid among processes equally along x-axis
	_sizes(0) = numberOfNodesAlongXPerOneCore(task);
	if (Mpi::Rank() == Mpi::Size() - 1) {
	// in order to keep specified in task number of nodes
		_sizes(0) = task.sizes.at(0) -
		            numberOfNodesAlongXPerOneCore(task) * (Mpi::Size() - 1);
	}
	
	return _sizes;
}


template<int Dimensionality>
typename Engine<Dimensionality>::RealD
Engine<Dimensionality>::
calculateH(const Task::CubicGrid& task) {
	
	RealD _h;
	if (!task.lengths.empty() && !task.sizes.empty()) {
		assert_true(task.h.empty());
		assert_eq(task.sizes.size(), DIMENSIONALITY);
		assert_eq(task.lengths.size(), DIMENSIONALITY);
		
		RealD taskLength; 
		taskLength.copyFrom(task.lengths);
		IntD taskSizes; 
		taskSizes.copyFrom(task.sizes);
		
		_h = linal::plainDivision(taskLength, taskSizes - IntD::Ones());
	
	} else {
		assert_eq(task.h.size(), DIMENSIONALITY);
		_h.copyFrom(task.h);
	}
	
	return _h;
}


template<int Dimensionality>
typename Engine<Dimensionality>::RealD
Engine<Dimensionality>::
calculateStartR(const Task::CubicGrid &task) {
	RealD _startR = RealD::Zeros();
	
	if (!task.startR.empty()) {
		assert_eq(task.startR.size(), DIMENSIONALITY);
		_startR.copyFrom(task.startR);
	}
	
	if (Mpi::ForceSequence()) {
		return _startR;
	}
	
	// MPI - divide the grid among processes equally along x-axis
	_startR(0) += Mpi::Rank() * numberOfNodesAlongXPerOneCore(task) *
			calculateH(task)(0);
	
	return _startR;
}



template class Engine<1>;
template class Engine<2>;
template class Engine<3>;

