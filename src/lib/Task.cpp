#include "Task.hpp"

using namespace gcm;

Task::Task() {
	xLength = 2.0;
	yLength = 3.0;
	zLength = 0.0;

	accuracyOrder = 2; // order of accuracy of spatial interpolation

	X = 51; // number of nodes along x direction
	Y = 51; // number of nodes along y direction
	Z = 1; // number of nodes along z direction

	rho0 = 8.0; // default density
	lambda0 = 12e+4; // default Lame parameter
	mu0 = 77e+3; // default Lame parameter

	CourantNumber = 1.2; // number from Courant–Friedrichs–Lewy condition
	numberOfSnaps = 50;

	initialConditions = InitialConditions::ExplosionAtTheLeft;

	borderConditions.at("left") = BorderConditions::FreeBorder;
	/*borderConditions.at("right") = BorderConditions::FreeBorder;
	borderConditions.at("bottom") = BorderConditions::FreeBorder;
	borderConditions.at("up") = BorderConditions::FreeBorder;*/
}