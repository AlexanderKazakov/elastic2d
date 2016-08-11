#include <lib/Engine.hpp>
#include <lib/GlobalVariables.hpp>
#include <lib/AbstractFactory.hpp>
#include <lib/numeric/solvers/Solver.hpp>
#include <lib/util/snapshot/Snapshotter.hpp>

#include <limits>

using namespace gcm;


Engine::
Engine(const Task& task_) : task(task_) {
	
	LOG_INFO("Start Engine");
	Clock::setZero();
	Mpi::initialize(task.globalSettings.forceSequence);
	
	assert_gt(task.statements.size(), 0);
	
	globalScene = Factory::createGlobalScene(task, this);
	
	for (const auto& body : task.bodies) {
		auto factory = Factory::create(task, body.second);
		
		Solver* solver = factory->createSolver(task, body.first, globalScene);
		
		std::vector<Snapshotter*> snapshotters;
		for (const auto snapId : task.globalSettings.snapshottersId) {
			Snapshotter* snapshotter = factory->createSnapshotter(snapId);
			snapshotter->initialize(task);
			snapshotters.push_back(snapshotter);
		}
		
		bodies.insert({body.first, {solver, snapshotters}});
	}
	
	globalScene->afterGridsConstruction(task);
}


Engine::~Engine() {
	for (const auto body : bodies) {
		for (Snapshotter* snapshotter : body.second.snapshotters) {
			delete snapshotter;
		}
		delete body.second.solver;
	}
	delete globalScene;
}


void Engine::
run() {
	for (const Statement& statement : task.statements) {
		LOG_INFO("Start statement " << statement.id);
		beforeStatement(statement);
		runStatement();
	}
}


void Engine::
beforeStatement(const Statement& statement) {
	Clock::setZero();
	for (const auto body : bodies) {
		body.second.solver->beforeStatement(statement);
		for (Snapshotter* snapshotter : body.second.snapshotters) {
			snapshotter->beforeStatement(statement);
			snapshotter->snapshot(body.second.solver->getAbstractMesh(), 0);
		}
	}
	
	estimateTimeStep();
	
	requiredTime = statement.globalSettings.numberOfSnaps *
	               statement.globalSettings.stepsPerSnap * Clock::TimeStep();
	if (statement.globalSettings.numberOfSnaps <= 0) {
		requiredTime = statement.globalSettings.requiredTime;
	}
	assert_gt(requiredTime, 0);
}


void Engine::
runStatement() {
	
	int step = 0;
	Utils::seedRand();
	
	while (Clock::Time() < requiredTime) {
		estimateTimeStep();
		
		LOG_INFO("Start next time step. Time = " << Clock::Time()
				<< ", TimeStep = " << Clock::TimeStep());
		nextTimeStep();
		step++; Clock::tickTack();
		
		for (const auto body : bodies) {
			for (Snapshotter* snapshotter : body.second.snapshotters) {
				snapshotter->snapshot(body.second.solver->getAbstractMesh(), step);
			}
		}
	}
	
	for (const auto body : bodies) {
		body.second.solver->afterStatement();
		for (Snapshotter* snapshotter : body.second.snapshotters) {
			snapshotter->afterStatement();
		}
	}
}


void Engine::
nextTimeStep() {
	globalScene->nextTimeStep();
	for (const auto body : bodies) {
		body.second.solver->afterStages(Clock::TimeStep());
	}
}


void Engine::
estimateTimeStep() {
	/// minimal among all solvers
	real minimalTimeStep = std::numeric_limits<real>::max();
	for (const auto body : bodies) {
		real solverTimeStep = body.second.solver->calculateTimeStep();
		if (solverTimeStep < minimalTimeStep) {
			minimalTimeStep = solverTimeStep;
		}
	}
	Clock::timeStep = minimalTimeStep;
}


