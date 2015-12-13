#include "Task.hpp"


Task::Task() {
	xLength = 2.0;
	yLength = 3.0;

	accuracyOrder = 2; // order of accuracy of spatial interpolation

	X = 4001; // number of nodes along x direction
	Y = 8001; // number of nodes along y direction

	rho0 = 8.0; // default density
	lambda0 = 12e+4; // default Lame parameter
	mu0 = 77e+3; // default Lame parameter

	CourantNumber = 1.2; // number from Courant–Friedrichs–Lewy condition
	numberOfSnaps = 40;

	initialConditions = InitialConditions::PWaveX;

	borderConditions.at("left") = BorderConditions::FreeBorder;
	/*borderConditions.at("right") = BorderConditions::FreeBorder;
	borderConditions.at("bottom") = BorderConditions::FreeBorder;
	borderConditions.at("up") = BorderConditions::FreeBorder;*/
}
