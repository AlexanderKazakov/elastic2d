#include <lib/Engine.hpp>

using namespace gcm;

void Engine::initialize(const Task &task) {
	assert_true(solver != nullptr);
	LOG_INFO("Start initialization");
	solver->initialize(task);
	if (snapshotter) snapshotter->initialize(task);

	real tau = solver->calculateTau();
	requiredTime = task.numberOfSnaps * task.stepsPerSnap * tau;
	if (task.numberOfSnaps == 0) requiredTime = task.T;
	assert_gt(requiredTime, 0);
}

void Engine::run() {
	LOG_INFO("Start calculations");
	int step = 0;
	if (snapshotter) snapshotter->snapshot(solver->getGrid(), step);

	while (solver->getCurrentTime() < requiredTime) {
		solver->nextTimeStep();
		step++;
		if (snapshotter) snapshotter->snapshot(solver->getGrid(), step);
	}
}
