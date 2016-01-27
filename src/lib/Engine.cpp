#include <lib/Engine.hpp>

using namespace gcm;

template<class TNode>
void Engine<TNode>::initialize(const Task &task) {
	LOG_INFO("Start initialization");
	solver = new DefaultSolver<TNode>();
	solver->initialize(task);
	snapshotter = new VtkTextStructuredSnapshotter<TNode>();
	snapshotter->initialize(task);

	real tau = solver->calculateTau();
	requiredTime = task.numberOfSnaps * tau;
	if (task.numberOfSnaps == 0) requiredTime = task.T;
	assert_gt(requiredTime, 0);
}

template<class TNode>
Engine<TNode>::~Engine() {
	delete solver;
	delete snapshotter;
}

template<class TNode>
void Engine<TNode>::run() {
	LOG_INFO("Start calculations");
	int step = 0;
	snapshotter->snapshot(solver->mesh, step);

	while (solver->currentTime < requiredTime) {
		solver->nextTimeStep();
		step++;
		snapshotter->snapshot(solver->mesh, step);
	}
}


template class Engine<IdealElastic1DNode>;
template class Engine<IdealElastic2DNode>;
template class Engine<IdealElastic3DNode>;