#include <libgcm/engine/AbstractEngine.hpp>
#include <libgcm/engine/GlobalVariables.hpp>
#include <libgcm/util/snapshot/Snapshotter.hpp>


using namespace gcm;


AbstractEngine::AbstractEngine(const Task& task) :
		CourantNumber(task.globalSettings.CourantNumber) {
	LOG_INFO("Start initialization ...");
	Clock::setZero();
	Mpi::initialize(task.globalSettings.forceSequence);
}


void AbstractEngine::afterGridsConstruction(const Task& task) {
	Clock::timeStep = estimateTimeStep();
	requiredTime = Clock::TimeStep() *
			task.globalSettings.numberOfSnaps *
			task.globalSettings.stepsPerSnap;
	if (task.globalSettings.numberOfSnaps <= 0) {
		requiredTime = task.globalSettings.requiredTime;
	}
	assert_gt(requiredTime, 0);
	
	writeSnapshots(0);
}


void AbstractEngine::run() {
	
	int step = 0;
	Utils::seedRand();
	
	while (Clock::Time() < requiredTime) {
		Clock::timeStep = estimateTimeStep();
		
		LOG_INFO("Start step " << step
				<< ". Time = " << Clock::Time()
				<< ". TimeStep = " << Clock::TimeStep());
		nextTimeStep();
		step++; Clock::tickTack();
		
		writeSnapshots(step);
	}
}


