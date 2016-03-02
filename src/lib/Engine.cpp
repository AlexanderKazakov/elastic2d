#include <lib/Engine.hpp>

using namespace gcm;

void Engine::initialize(const Task &task_) {
	LOG_INFO("Start initialization");
	task = task_;
	assert_gt(task.statements.size(), 0);
	assert_true(solver != nullptr);
	solver->initialize(task);
	for (auto snap : snapshotters) {
		snap->initialize(task);
	}
}

void Engine::run() {
	for (const auto& statement : task.statements) {
		std::cout << "Start statement " << statement.id << std::endl;
		beforeStatement(statement);
		runStatement();
	}
}

void Engine::beforeStatement(const Statement& statement) {
	solver->beforeStatement(statement);
	real tau = solver->calculateTau();
	std::cout << "Tau = " << tau << std::endl;
	requiredTime = statement.numberOfSnaps * statement.stepsPerSnap * tau;
	if (statement.numberOfSnaps == 0) requiredTime = statement.T;
	assert_gt(requiredTime, 0);
	
	for (auto snap : snapshotters) {
		snap->beforeStatement(statement);
		snap->snapshot(solver->getGrid(), 0);
	}
}

void Engine::runStatement() {
	LOG_INFO("Start calculations");

	int step = 1;
	while (solver->getCurrentTime() < requiredTime) {
		solver->nextTimeStep();
		step++;
		for (auto snap : snapshotters) {
			snap->snapshot(solver->getGrid(), step);
		}
	}
	
	for (auto snap : snapshotters) {
		snap->afterStatement();
	}
	solver->afterStatement();
}
