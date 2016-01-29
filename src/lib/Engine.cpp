#include <lib/Engine.hpp>
#include <lib/numeric/solvers/DefaultSolver.hpp>
#include <lib/util/snapshot/VtkTextStructuredSnapshotter.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

void Engine::initialize(const Task &task) {
	LOG_INFO("Start initialization");
	solver = new DefaultSolver<StructuredGrid<Elastic2DModel>>();
	solver->initialize(task);
	snapshotter = new VtkTextStructuredSnapshotter<StructuredGrid<Elastic2DModel>>();
	snapshotter->initialize(task);

	real tau = solver->calculateTau();
	requiredTime = task.numberOfSnaps * tau;
	if (task.numberOfSnaps == 0) requiredTime = task.T;
	assert_gt(requiredTime, 0);
}

Engine::~Engine() {
	delete solver;
	delete snapshotter;
}

void Engine::run() {
	LOG_INFO("Start calculations");
	int step = 0;
	snapshotter->snapshot(solver->getGrid(), step);

	while (solver->getCurrentTime() < requiredTime) {
		solver->nextTimeStep();
		step++;
		snapshotter->snapshot(solver->getGrid(), step);
	}
}
