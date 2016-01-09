#include "Task.hpp"

using namespace gcm;

Task::Task() {
	xLength = 2.0;
	yLength = 3.0;
	zLength = 2.0;

	accuracyOrder = 2; // order of accuracy of spatial interpolation

	X = 21; // number of nodes along x direction
	Y = 21; // number of nodes along y direction
	Z = 1; // number of nodes along z direction

	rho0 = 8.0; // default density
	lambda0 = 12e+4; // default Lame parameter
	mu0 = 77e+3; // default Lame parameter

	CourantNumber = 1.2; // number from Courant–Friedrichs–Lewy condition
	numberOfSnaps = 20;

	initialConditions = InitialConditions::Explosion;

	borderConditions.at(Border::X_LEFT) = BorderConditions::FreeBorder;

}