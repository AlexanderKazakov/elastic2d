#include <lib/Engine.hpp>
#include <lib/AbstractFactory.hpp>
#include <lib/numeric/solvers/Solver.hpp>
#include <lib/util/snapshot/Snapshotter.hpp>

using namespace gcm;

void Engine::initialize(const Task &task_) {
	LOG_INFO("Start initialization");
	task = task_;
	assert_gt(task.statements.size(), 0);
	
	forceSequence = task.globalSettings.forceSequence;
	CubicGrid::preprocessTask(task.cubicGrid);
	
	auto factory = Factory::create(task);
	solver = factory->createSolver(task);
	
	for (auto snapId : task.snapshottersId) {
		auto snapshotter = factory->createSnapshotter(snapId);
		snapshotter->initialize(task);
		snapshotters.push_back(snapshotter);
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
	currentTime = 0;
	solver->beforeStatement(statement);

	requiredTime = statement.globalSettings.numberOfSnaps * 
	               statement.globalSettings.stepsPerSnap * estimateTimeStep();
	if (statement.globalSettings.numberOfSnaps <= 0) {
		requiredTime = statement.globalSettings.requiredTime;
	}
	assert_gt(requiredTime, 0);
	
	for (auto snap : snapshotters) {
		snap->beforeStatement(statement);
		snap->snapshot(solver->getActualGrid(), 0);
	}
}

void Engine::runStatement() {
	LOG_INFO("Start calculations");

	int step = 1; currentTime = 0;
	while (getCurrentTime() < requiredTime) {
		real timeStep = estimateTimeStep();
		solver->nextTimeStep(timeStep);
		for (auto snap : snapshotters) {
			snap->snapshot(solver->getActualGrid(), step);
		}
		step++; currentTime += timeStep;
	}
	
	for (auto snap : snapshotters) {
		snap->afterStatement();
	}
	solver->afterStatement();
}

real Engine::estimateTimeStep() {
	// when more than one solver will be ...
	return solver->calculateTimeStep();
}
