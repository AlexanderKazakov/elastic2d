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
		real zLength = 0.0;

		int accuracyOrder = 0; // order of accuracy of spatial interpolation

		int X = 0; // number of nodes along x direction
		int Y = 0; // number of nodes along y direction
		int Z = 0; // number of nodes along z direction

		real rho0 = 0.0; // default density
		real lambda0 = 0.0; // default Lame parameter
		real mu0 = 0.0; // default Lame parameter

		real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition

		int numberOfSnaps = 0; // how many snaps to calculate
		real T = 0.0; // optional, required time if (numberOfSnaps == 0)

		bool splittingSecondOrder = false;
		bool enableSnapshotting = false;
		bool forceSequence = false; // if true make the grid thinking that the number of
		// processes is one, even if it isn's so actually (for testing purposes)

		InitialConditions initialConditions = InitialConditions::Zero;

		std::map <Border, BorderConditions> borderConditions = {
				{Border::X_LEFT,  BorderConditions::NonReflection},
				{Border::X_RIGHT, BorderConditions::NonReflection},
				{Border::Y_LEFT,  BorderConditions::NonReflection},
				{Border::Y_RIGHT, BorderConditions::NonReflection},
				{Border::Z_LEFT,  BorderConditions::NonReflection},
				{Border::Z_RIGHT, BorderConditions::NonReflection}};

		/* ------------------ Properties and conditions (end) ------------------ */

		Task();
	};
}

#endif //LIBGCM_TASK_HPP
