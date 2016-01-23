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

	real rho0 = 8.0; // default density
	real lambda0 = 12e+4; // default Lame parameter
	real mu0 = 77e+3; // default Lame parameter
	material = IsotropicMaterial(rho0, lambda0, mu0);

	CourantNumber = 1.2; // number from Courant–Friedrichs–Lewy condition
	numberOfSnaps = 20;

}