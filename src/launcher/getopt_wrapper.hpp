#include <getopt.h>
#include <string>


/**
 * @see getopt_long_only documentation
 * @param argc program arguments number
 * @param argv program arguments array
 */
inline void getopt_wrapper(int argc, char** argv,
		std::string& taskId,
		std::string& outputDirectory) {
	taskId = "";
	outputDirectory = "";
	
	while (true) {
	// through all options
		static struct option long_options[] = {
			{"task",    required_argument, 0, 't'},
			{"out",     required_argument, 0, 'o'},
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
		int c = getopt_long_only(argc, argv, "t:o:", long_options, &option_index);
		if (c == -1) { break; } // end of the options
		
		switch (c) {
			case 't':
				taskId = optarg;
				break;
			
			case 'o':
				outputDirectory = optarg;
				break;
			
			default:
				break;
		}
	}
}
