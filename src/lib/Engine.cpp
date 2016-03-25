#include <lib/Engine.hpp>
#include <lib/AbstractFactory.hpp>
#include <lib/numeric/solvers/Solver.hpp>
#include <lib/util/snapshot/Snapshotter.hpp>

using namespace gcm;

real Clock::time = 0;
int Mpi::rank = 0;
int Mpi::size = 1;
bool Mpi::forceSequence = false;

Engine::Engine(const Task &task_) :
		task(task_) {
	LOG_INFO("Start Engine");
	Clock::time = 0;
	Mpi::initialize(task.globalSettings.forceSequence);
	
	assert_gt(task.statements.size(), 0);
	if (task.gridId == Grids::T::CUBIC) {
		CubicGrid::preprocessTask(task.cubicGrid);
	}
	
	auto factory = Factory::create(task);
	solver = factory->createSolver(task);
	
	for (auto snapId : task.snapshottersId) {
		auto snapshotter = factory->createSnapshotter(snapId);
		snapshotter->initialize(task);
		snapshotters.push_back(snapshotter);
	}
}

Engine::~Engine() {
	for (auto snap : snapshotters) {
		delete snap;
	}
	delete solver;
}

void Engine::run() {
	for (const auto& statement : task.statements) {
		LOG_INFO("Start statement " << statement.id);
		beforeStatement(statement);
		runStatement();
	}
}

void Engine::beforeStatement(const Statement& statement) {
	Clock::time = 0;
	solver->beforeStatement(statement);

	requiredTime = statement.globalSettings.numberOfSnaps * 
	               statement.globalSettings.stepsPerSnap * estimateTimeStep();
	if (statement.globalSettings.numberOfSnaps <= 0) {
		requiredTime = statement.globalSettings.requiredTime;
	}
	assert_gt(requiredTime, 0);
	
	for (auto snap : snapshotters) {
		snap->beforeStatement(statement);
		snap->snapshot(solver->getActualMesh(), 0);
	}
}

void Engine::runStatement() {
	int step = 0;
	while (Clock::Time() < requiredTime) {
		real timeStep = estimateTimeStep();
		solver->nextTimeStep(timeStep);
		step++; Clock::time += timeStep;
		for (auto snap : snapshotters) {
			snap->snapshot(solver->getActualMesh(), step);
		}
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
