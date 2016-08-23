#include <lib/Engine.hpp>
#include <lib/GlobalVariables.hpp>
#include <lib/AbstractFactory.hpp>
#include <lib/numeric/solvers/Solver.hpp>
#include <lib/util/snapshot/Snapshotter.hpp>

#include <limits>

using namespace gcm;


Engine::
Engine(const Task& task) {
	
	LOG_INFO("Start Engine");
	Clock::setZero();
	Mpi::initialize(task.globalSettings.forceSequence);
	
	globalScene = Factory::createGlobalScene(task, this);
	
	for (const auto& body : task.bodies) {
		auto factory = Factory::create(task, body.second);
		
		Solver* solver = factory->createSolver(task, body.first, globalScene);
		
		std::vector<Snapshotter*> snapshotters;
		for (const auto snapId : task.globalSettings.snapshottersId) {
			Snapshotter* snapshotter = factory->createSnapshotter(task, snapId);
			snapshotters.push_back(snapshotter);
		}
		
		bodies.insert({body.first, {solver, snapshotters}});
	}
	
	globalScene->afterGridsConstruction(task);
	
	Clock::setZero();
	for (const auto body : bodies) {
		for (Snapshotter* snapshotter : body.second.snapshotters) {
			snapshotter->snapshot(body.second.solver->getAbstractMesh(), 0);
		}
	}
	
	estimateTimeStep();
	
	requiredTime = task.globalSettings.numberOfSnaps *
	               task.globalSettings.stepsPerSnap * Clock::TimeStep();
	if (task.globalSettings.numberOfSnaps <= 0) {
		requiredTime = task.globalSettings.requiredTime;
	}
	assert_gt(requiredTime, 0);
	
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


void Engine::run() {
	
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


