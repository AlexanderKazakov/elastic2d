#include "Task.hpp"


Task::Task() {
	xLength = 1.0;
	yLength = 1.0;

	accuracyOrder = 2; // order of accuracy of spatial interpolation

	X = 51; // number of nodes along x direction
	Y = 51; // number of nodes along y direction

	rho0 = 8.0; // default density
	lambda0 = 12e+4; // default Lame parameter
	mu0 = 77e+3; // default Lame parameter

	CourantNumber = 1.2; // number from Courant–Friedrichs–Lewy condition
	numberOfSnaps = 20;

	initialConditions = InitialConditions::PWaveX;
}