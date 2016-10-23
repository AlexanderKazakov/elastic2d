#include <csignal>
#include <chrono>

#include <libgcm/engine/EngineFactory.hpp>

#include <ndi/ndi.hpp>


using namespace gcm;


inline void signalHandler(int a) {
	std::cout << "Receive signal " << a << ". Terminate." << std::endl;
	abort();
}


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	signal(SIGFPE, signalHandler);
	USE_AND_INIT_LOGGER("gcm.main");
	
	if (argc != 3) {
		LOG_FATAL("Specify composite width in mm and number of nodes per 1mm along Y-axis. "
				<< "Example:\n"
				<< "\"" << argv[0] << " 6 10\" - six millimeters and ten nodes per millimeter\n"
				<< "\"" << argv[0] << " 0 10\" - empty (without composite)\n"
				<< "\"" << argv[0] << " 100 10\" - non-reflection from the bottom");
		return -1;
	}
	
	try {
		const real thickness = std::atof(argv[1]);
		assert_ge(thickness, 0);
		const int nodesPerMm = std::atoi(argv[2]);
		assert_ge(nodesPerMm, 1);
		Task task = ndi(thickness, nodesPerMm);
		
		auto t1 = std::chrono::high_resolution_clock::now();
		createEngine(task)->run();
		auto t2 = std::chrono::high_resolution_clock::now();
		SUPPRESS_WUNUSED(t1); SUPPRESS_WUNUSED(t2);
		
		LOG_INFO("Time of calculation, microseconds = " << 
			std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count());
		
	} catch (Exception& e) {
		LOG_FATAL(e.what());
		MPI_Finalize();
		return -1;
		
	} catch (std::exception& e) {
		LOG_FATAL("std::exception was thrown: " << e.what() << "Exit");
		MPI_Finalize();
		return -2;
		
	} catch ( ... ) {
		LOG_FATAL("Unknown exception was thrown. Exit");
		MPI_Finalize();
		return -3;
	}
	
	MPI_Finalize();
	return 0;
}

