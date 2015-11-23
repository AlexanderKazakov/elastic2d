#ifndef ELASTIC2D_TASK_HPP
#define ELASTIC2D_TASK_HPP

#include "lib/config.hpp"


class Task {
public:
	/* ------------------ Properties and conditions ------------------ */

	real xLength = 0.0;
	real yLength = 0.0;

	uint accuracyOrder = 0; // order of accuracy of spatial interpolation

	uint X = 0; // number of nodes along x direction
	uint Y = 0; // number of nodes along y direction

	real rho0 = 0.0; // default density
	real lambda0 = 0.0; // default Lame parameter
	real mu0 = 0.0; // default Lame parameter

	real CourantNumber = 0.0; // number from Courant–Friedrichs–Lewy condition

	uint numberOfSnaps = 0; // how many snaps to calculate
	real T = -1.0; // optional, required time if (numberOfSnaps == 0)

	InitialConditions initialConditions = InitialConditions::Zero;

	/* ------------------ Properties and conditions (end) ------------------ */

	Task();
};


#endif //ELASTIC2D_TASK_HPP
