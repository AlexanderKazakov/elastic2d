#include <lib/Engine.hpp>
#include <lib/AbstractFactory.hpp>
#include <lib/numeric/solvers/Solver.hpp>
#include <lib/util/snapshot/Snapshotter.hpp>

#include <limits>

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
	
	globalScene = Factory::createGlobalScene(task);
	
	for (const Task::Body& body : task.bodies) {
		auto factory = Factory::create(task, body);
		
		Solver* solver = factory->createSolver(task, body, globalScene);
		
		std::vector<Snapshotter*> snapshotters;
		for (const auto snapId : task.globalSettings.snapshottersId) {
			Snapshotter* snapshotter = factory->createSnapshotter(snapId);
			snapshotter->initialize(task);
			snapshotters.push_back(snapshotter);
		}
		
		bodies.push_back({solver, snapshotters});
	}
	
}


Engine::~Engine() {
	for (Body& body : bodies) {
		for (Snapshotter* snapshotter : body.snapshotters) {
			delete snapshotter;
		}
		delete body.solver;
	}
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
	for (Body& body : bodies) {
		body.solver->beforeStatement(statement);
		for (Snapshotter* snapshotter : body.snapshotters) {
			snapshotter->beforeStatement(statement);
			snapshotter->snapshot(body.solver->getActualMesh(), 0);
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
	
	while (Clock::Time() < requiredTime) {
		estimateTimeStep();
		
		LOG_INFO("Start next time step. Time = " << Clock::Time()
				<< ", TimeStep = " << Clock::TimeStep());
		for (Body& body : bodies) {
			body.solver->nextTimeStep();
		}
		step++; Clock::tickTack();
		
		for (Body& body : bodies) {
			for (Snapshotter* snapshotter : body.snapshotters) {
				snapshotter->snapshot(body.solver->getActualMesh(), step);
			}
		}
	}
	
	for (Body& body : bodies) {
		body.solver->afterStatement();
		for (Snapshotter* snapshotter : body.snapshotters) {
			snapshotter->afterStatement();
		}
	}
}


void Engine::
estimateTimeStep() {
	/// minimal among all solvers
	real minimalTimeStep = std::numeric_limits<real>::max();
	for (Body& body : bodies) {
		real solverTimeStep = body.solver->calculateTimeStep();
		if (solverTimeStep < minimalTimeStep) {
			minimalTimeStep = solverTimeStep;
		}
	}
	Clock::timeStep = minimalTimeStep;
}


