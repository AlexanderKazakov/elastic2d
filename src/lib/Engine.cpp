#include <lib/Engine.hpp>
#include <lib/AbstractFactory.hpp>
#include <lib/numeric/solvers/Solver.hpp>
#include <lib/util/snapshot/Snapshotter.hpp>

using namespace gcm;

real Clock::time = 0;
real Clock::timeStep = 0;
int Mpi::rank = 0;
int Mpi::size = 1;
bool Mpi::forceSequence = false;

Engine::
Engine(const Task& task_) :
	task(task_) {
	LOG_INFO("Start Engine");
	Clock::setZero();
	Mpi::initialize(task.globalSettings.forceSequence);

	assert_gt(task.statements.size(), 0);

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


void Engine::
run() {
	for (const auto& statement : task.statements) {
		LOG_INFO("Start statement " << statement.id);
		beforeStatement(statement);
		runStatement();
	}
}


void Engine::
beforeStatement(const Statement& statement) {
	Clock::setZero();
	solver->beforeStatement(statement);
	estimateTimeStep();

	requiredTime = statement.globalSettings.numberOfSnaps *
	               statement.globalSettings.stepsPerSnap * Clock::TimeStep();
	if (statement.globalSettings.numberOfSnaps <= 0) {
		requiredTime = statement.globalSettings.requiredTime;
	}
	assert_gt(requiredTime, 0);

	for (auto snap : snapshotters) {
		snap->beforeStatement(statement);
		snap->snapshot(solver->getActualMesh(), 0);
	}
}


void Engine::
runStatement() {
	int step = 0;
	while (Clock::Time() < requiredTime) {
		estimateTimeStep();
		
		LOG_DEBUG("Start next time step. Time = " << Clock::Time());
		solver->nextTimeStep();
		step++; Clock::tickTack();
		
		for (auto snap : snapshotters) {
			snap->snapshot(solver->getActualMesh(), step);
		}
	}

	for (auto snap : snapshotters) {
		snap->afterStatement();
	}
	solver->afterStatement();
}


void Engine::
estimateTimeStep() {
	// when more solvers will be ..
	Clock::timeStep = solver->calculateTimeStep();
}


