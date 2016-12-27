#include <libgcm/engine/AbstractEngine.hpp>
#include <libgcm/engine/GlobalVariables.hpp>
#include <libgcm/util/snapshot/Snapshotter.hpp>


using namespace gcm;


AbstractEngine::AbstractEngine(const Task& task) :
		CourantNumber(task.globalSettings.CourantNumber),
		verboseTimeSteps(task.globalSettings.verboseTimeSteps) {
	LOG_INFO("Start initialization ...");
	Clock::setZero();
	Mpi::initialize(task.globalSettings.forceSequence);
}


void AbstractEngine::afterConstruction(const Task& task) {
	Clock::timeStep = estimateTimeStep();
	requiredTime = Clock::TimeStep() *
			task.globalSettings.numberOfSnaps *
			task.globalSettings.stepsPerSnap;
	if (task.globalSettings.numberOfSnaps <= 0) {
		requiredTime = task.globalSettings.requiredTime;
	}
	assert_gt(requiredTime, 0);
}


void AbstractEngine::run() {
	int step = 0;
	writeSnapshots(step);
	Utils::seedRand();
	
	while (Clock::Time() < requiredTime) {
		Clock::timeStep = estimateTimeStep();
		if (verboseTimeSteps) {
			LOG_INFO("Start step " << step
					<< ". Time = " << Clock::Time()
					<< ". TimeStep = " << Clock::TimeStep());
		}
		nextTimeStep();
		step++; Clock::tickTack();
		writeSnapshots(step);
	}
}


