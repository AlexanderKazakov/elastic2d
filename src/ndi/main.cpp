#include <chrono>

#include <libgcm/engine/EngineFactory.hpp>

#include <ndi/ndi.hpp>


using namespace gcm;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	USE_AND_INIT_LOGGER("gcm.main");
	
//	Task task = test(); 
	Task task = ndi();
	try {
		auto t1 = std::chrono::high_resolution_clock::now();
		createEngine(task)->run();
		auto t2 = std::chrono::high_resolution_clock::now();
		SUPPRESS_WUNUSED(t1); SUPPRESS_WUNUSED(t2);

		LOG_INFO("Time of calculation, microseconds = " << 
			std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count());
	} catch (Exception e) {
		LOG_FATAL(e.what());
	}

	MPI_Finalize();
	return 0;
}

