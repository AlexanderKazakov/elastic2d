#ifndef LIBGCM_TASK_HPP
#define LIBGCM_TASK_HPP

#include <map>

#include "lib/util/Types.hpp"
#include "lib/config.hpp"

namespace gcm {
	class Task {
	public:
		/* ------------------ Properties and conditions ------------------ */

		real xLength = 0.0;
		real yLength = 0.0;

		int accuracyOrder = 0; // order of accuracy of spatial interpolation

		int X = 0; // number of nodes along x direction
		int Y = 0; // number of nodes along y direction

		real rho0 = 0.0; // default density
		real lambda0 = 0.0; // default Lame parameter
		real mu0 = 0.0; // default Lame parameter

		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition

		int numberOfSnaps = 0; // how many snaps to calculate
		real T = 0.0; // optional, required time if (numberOfSnaps == 0)

		InitialConditions initialConditions = InitialConditions::Zero;

		std::map <std::string, BorderConditions> borderConditions = {
				{"left",   BorderConditions::NonReflection},
				{"right",  BorderConditions::NonReflection},
		        {"bottom", BorderConditions::NonReflection},
		        {"up",     BorderConditions::NonReflection}};

		/* ------------------ Properties and conditions (end) ------------------ */

		Task();
	};
}

#endif //LIBGCM_TASK_HPP
