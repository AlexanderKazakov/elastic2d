#include <lib/Engine.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TGrid>
void Engine<TGrid>::initialize(const Task &task) {
	LOG_INFO("Start initialization");
	solver = new DefaultSolver<TGrid>();
	solver->initialize(task);
	snapshotter = new VtkTextStructuredSnapshotter<TGrid>();
	snapshotter->initialize(task);

	real tau = solver->calculateTau();
	requiredTime = task.numberOfSnaps * tau;
	if (task.numberOfSnaps == 0) requiredTime = task.T;
	assert_gt(requiredTime, 0);
}

template<class TGrid>
Engine<TGrid>::~Engine() {
	delete solver;
	delete snapshotter;
}

template<class TGrid>
void Engine<TGrid>::run() {
	LOG_INFO("Start calculations");
	int step = 0;
	snapshotter->snapshot(solver->mesh, step);

	while (solver->currentTime < requiredTime) {
		solver->nextTimeStep();
		step++;
		snapshotter->snapshot(solver->mesh, step);
	}
}



template class Engine<StructuredGrid<Elastic1DModel>>;
template class Engine<StructuredGrid<Elastic2DModel>>;
template class Engine<StructuredGrid<Elastic3DModel>>;

template class Engine<StructuredGrid<PlasticFlow1DModel>>;
template class Engine<StructuredGrid<PlasticFlow2DModel>>;
template class Engine<StructuredGrid<PlasticFlow3DModel>>;

template class Engine<StructuredGrid<OrthotropicElastic3DModel>>;
template class Engine<StructuredGrid<OrthotropicPlasticFlow3DModel>>;
