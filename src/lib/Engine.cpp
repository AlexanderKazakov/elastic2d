#include <lib/Engine.hpp>

using namespace gcm;

void Engine::initialize(const Task &task) {
	assert_true(solver != nullptr);
	LOG_INFO("Start initialization");
	solver->initialize(task);
	for (auto snap : snapshotters) {
		snap->initialize(task);
	}
	real tau = solver->calculateTau();
	requiredTime = task.numberOfSnaps * task.stepsPerSnap * tau;
	if (task.numberOfSnaps == 0) requiredTime = task.T;
	assert_gt(requiredTime, 0);
}

void Engine::run() {
	LOG_INFO("Start calculations");
	int step = 0;
	for (auto snap : snapshotters) {
		snap->beforeCalculation(solver);
		snap->snapshot(solver->getGrid(), step);
	}

	while (solver->getCurrentTime() < requiredTime) {
		solver->nextTimeStep();
		step++;
		for (auto snap : snapshotters) {
			snap->snapshot(solver->getGrid(), step);
		}
	}
	
	for (auto snap : snapshotters) {
		snap->afterCalculation();
	}
}
