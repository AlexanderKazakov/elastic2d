#include <fstream>
#include <algorithm>
#include "lib/Mesh.hpp"

using namespace gcm;


void Mesh::initialize(const Task &task, const bool forceSequence) {

	rank = MPI::COMM_WORLD.Get_rank();
	numberOfWorkers = MPI::COMM_WORLD.Get_size();

	if (forceSequence) {rank = 0; numberOfWorkers = 1;}

	/* ------------------ Properties and conditions ------------------ */

	accuracyOrder = task.accuracyOrder; // order of accuracy of spatial interpolation

	X = task.X; // number of nodes along x direction

	// we divide the mesh among processes equally along y-axis
	int numberOfNodesAlongYPerOneCore = (int)std::round((real)task.Y / numberOfWorkers);
	Y = numberOfNodesAlongYPerOneCore; // number of nodes along y direction
	if (rank == numberOfWorkers - 1) Y = task.Y - numberOfNodesAlongYPerOneCore * (numberOfWorkers - 1);

	globalY = task.Y;

	h[0] = task.xLength / (X - 1); /* x spatial step */
	h[1] = task.yLength / (task.Y - 1); /* y spatial step */

	real c0 = sqrt((task.lambda0 + 2 * task.mu0) / task.rho0); // default acoustic velocity

	tau = task.CourantNumber * fmin(h[0], h[1]) / c0; // time step

	T = task.numberOfSnaps * tau; // required time
	if (task.numberOfSnaps == 0) T = task.T;

	initialConditions = task.initialConditions;
	borderConditions = task.borderConditions;

	/* ------------------ Properties and conditions (end) ------------------ */

	startY = rank * numberOfNodesAlongYPerOneCore;

	nodes = new Node[(Y + 2 * accuracyOrder) * (X + 2 * accuracyOrder)];

	for (int i = 0; i < (Y + 2 * accuracyOrder) * (X + 2 * accuracyOrder); ++i) {
		for (int k = 0; k < N; ++k) {
			nodes[i].u(k) = 0;
		}
		nodes[i].matrix = nullptr;
	}

	defaultMatrix = std::make_shared<PDEMatrices>(task.rho0, task.lambda0, task.mu0);
	for (int y = 0; y < Y; ++y) {
		for (int x = 0; x < X; ++x) {
			(*this)(y, x).matrix = defaultMatrix;
		}
	}

	applyInitialConditions();
}


Matrix Mesh::interpolateValuesAround(const int stage, const int y, const int x,
                                     const Vector &dx) const {

	Matrix ans;
	std::vector<Vector> src(accuracyOrder + 1);
	Vector res;
	for (int k = 0; k < N; k++) {
		findSourcesForInterpolation(stage, y, x, dx.get(k), src);
		interpolator.minMaxInterpolate(res, src, fabs(dx.get(k)) / h[stage]);
		ans.setColumn(k, res);
	}

	return ans;
}


void Mesh::findSourcesForInterpolation(const int stage, const int y, const int x,
                                       const real &dx, std::vector<Vector>& src) const {

	const int alongX = (stage == 0) * ( (dx > 0) ? 1 : -1 );
	const int alongY = (stage == 1) * ( (dx > 0) ? 1 : -1 );
	for (int k = 0; k < src.size(); k++) {
		src[k] = get(y + alongY * k, x + alongX * k).u;
	}

}


void Mesh::snapshot(int step) const {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d.vtk", "snaps/core", rank, "_snapshot", step);
	std::fstream f(buffer, std::ios::out);
	if (!f) {
		std::cerr << "Unable to open file " << buffer << std::endl;
		return;
	}
	f << "# vtk DataFile Version 3.0" << std::endl;
	f << "U data" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET STRUCTURED_POINTS" << std::endl;
	f << "DIMENSIONS " << X << " " << Y << " 1" << std::endl;
	f << "SPACING " << h[0] << " " << h[1] << " 1" << std::endl;
	f << "ORIGIN " << "0 " << startY * h[1] << " 0" << std::endl;
	f << "POINT_DATA " << X * Y << std::endl;

	// TODO - will it works with floats?
	f << "VECTORS V double" << std::endl;
	for (int y = 0; y < Y; y++)
		for (int x = 0; x < X; x++)
			f << get(y, x).get(NodeMap::Vx) << " " << get(y, x).get(NodeMap::Vy) << " 0" << std::endl;

	f << "SCALARS Sxx double" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int y = 0; y < Y; y++)
		for (int x = 0; x < X; x++)
			f << get(y, x).get(NodeMap::Sxx) << std::endl;

	f << "SCALARS Sxy double" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int y = 0; y < Y; y++)
		for (int x = 0; x < X; x++)
			f << get(y, x).get(NodeMap::Sxy) << std::endl;

	f << "SCALARS Syy double" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int y = 0; y < Y; y++)
		for (int x = 0; x < X; x++)
			f << get(y, x).get(NodeMap::Syy) << std::endl;

	f << "SCALARS pressure double" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int y = 0; y < Y; y++)
		for (int x = 0; x < X; x++)
			f << - (get(y, x).get(NodeMap::Sxx) + get(y, x).get(NodeMap::Syy)) / 2 << std::endl;

	f.close();
}


template<typename T>
static void put(std::fstream &f, const T value) {
	union {
		char buf[sizeof(T)];
		T val;
	} helper;
	helper.val = value;
	std::reverse(helper.buf, helper.buf + sizeof(T));
	f.write(helper.buf, sizeof(T));
}

void Mesh::_snapshot(int step) const {
	char buffer[50];
	sprintf(buffer, "%s%02d%s%05d.vtk", "snaps/core", rank, "_snapshot", step);
	std::fstream f(buffer, std::ios::out);
	if (!f) {
		std::cerr << "Unable to open file " << buffer << std::endl;
		return;
	}
	f << "# vtk DataFile Version 3.0" << std::endl;
	f << "U data" << std::endl;
	f << "BINARY" << std::endl;
	f << "DATASET STRUCTURED_POINTS" << std::endl;
	f << "DIMENSIONS " << X + 2 * accuracyOrder << " " << Y + 2 * accuracyOrder << " 1" << std::endl;
	f << "SPACING " << h[0] << " " << h[1] << " 1" << std::endl;
	f << "ORIGIN " << - accuracyOrder * h[0] << " "
	<< (rank * (globalY / numberOfWorkers + 2 * accuracyOrder) - accuracyOrder) * h[1]
	<< " 0" << std::endl;
	f << "POINT_DATA " << (X + 2 * accuracyOrder) * (Y + 2 * accuracyOrder) << std::endl;

	f << "SCALARS Vx double" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < (Y + 2 * accuracyOrder) * (X + 2 * accuracyOrder); ++i) {
		put(f, nodes[i].get(NodeMap::Vx));
	}

	f << "SCALARS Vy double" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < (Y + 2 * accuracyOrder) * (X + 2 * accuracyOrder); ++i) {
		put(f, nodes[i].get(NodeMap::Vy));
	}

	f << "SCALARS Sxx double" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < (Y + 2 * accuracyOrder) * (X + 2 * accuracyOrder); ++i) {
		put(f, nodes[i].get(NodeMap::Sxx));
	}

	f << "SCALARS Sxy double" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < (Y + 2 * accuracyOrder) * (X + 2 * accuracyOrder); ++i) {
		put(f, nodes[i].get(NodeMap::Sxy));
	}

	f << "SCALARS Syy double" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < (Y + 2 * accuracyOrder) * (X + 2 * accuracyOrder); ++i) {
		put(f, nodes[i].get(NodeMap::Syy));
	}

	f.close();
}


void Mesh::changeRheology(const real& rho2rho0, const real& lambda2lambda0, const real& mu2mu0) {

	auto newRheologyMatrix = std::make_shared<PDEMatrices>(rho2rho0 * defaultMatrix->rho,
	                                                       lambda2lambda0 * defaultMatrix->lambda,
	                                                       mu2mu0 * defaultMatrix->mu);

	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++) {
			if ((y + startY) * h[1] >= 0.5) {
				(*this)(y, x).matrix = newRheologyMatrix;
			}
		}
	}

	const real defaultAcousticVelocity = sqrt((defaultMatrix->lambda + 2 * defaultMatrix->mu) /
			                                          defaultMatrix->rho);
	const real newAcousticVelocity = sqrt((newRheologyMatrix->lambda + 2 * newRheologyMatrix->mu) / 
			                                      newRheologyMatrix->rho);

	if (newAcousticVelocity > defaultAcousticVelocity) {
		tau /= newAcousticVelocity / defaultAcousticVelocity;
	}

}


void Mesh::changeRheology2(const real& rho2rho0, const real& lambda2lambda0, const real& mu2mu0) {

	auto newRheologyMatrix = std::make_shared<PDEMatrices>(rho2rho0 * defaultMatrix->rho,
	                                                       lambda2lambda0 * defaultMatrix->lambda,
	                                                       mu2mu0 * defaultMatrix->mu);
	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++) {
			if ((y + startY) > 0.2 * globalY && (y + startY) < 0.6 * globalY &&
					x > 0.4 * X && x < 0.6 * X) {
				(*this)(y, x).matrix = newRheologyMatrix;
			}
		}
	}

	const real defaultAcousticVelocity = sqrt((defaultMatrix->lambda + 2 * defaultMatrix->mu) /
	                                          defaultMatrix->rho);
	const real newAcousticVelocity = sqrt((newRheologyMatrix->lambda + 2 * newRheologyMatrix->mu) /
	                                      newRheologyMatrix->rho);

	if (newAcousticVelocity > defaultAcousticVelocity) {
		tau /= newAcousticVelocity / defaultAcousticVelocity;
	}

}


void Mesh::applyBorderConditions() {
	if (borderConditions.at("left") == BorderConditions::FreeBorder) {
		for (int y = 0; y < Y; y++) {
			for (int i = 1; i <= accuracyOrder; i++) {
				(*this)(y, -i)(NodeMap::Vx)  =   get(y, i).get(NodeMap::Vx);
				(*this)(y, -i)(NodeMap::Vy)  =   get(y, i).get(NodeMap::Vy);
				(*this)(y, -i)(NodeMap::Sxx) = - get(y, i).get(NodeMap::Sxx);
				(*this)(y, -i)(NodeMap::Sxy) = - get(y, i).get(NodeMap::Sxy);
				(*this)(y, -i)(NodeMap::Syy) = - get(y, i).get(NodeMap::Syy);
			}
		}
	}
	if (borderConditions.at("right") == BorderConditions::FreeBorder) {
		for (int y = 0; y < Y; y++) {
			for (int i = 1; i <= accuracyOrder; i++) {
				(*this)(y, X - 1 + i)(NodeMap::Vx)  =   get(y, X - 1 - i).get(NodeMap::Vx);
				(*this)(y, X - 1 + i)(NodeMap::Vy)  =   get(y, X - 1 - i).get(NodeMap::Vy);
				(*this)(y, X - 1 + i)(NodeMap::Sxx) = - get(y, X - 1 - i).get(NodeMap::Sxx);
				(*this)(y, X - 1 + i)(NodeMap::Sxy) = - get(y, X - 1 - i).get(NodeMap::Sxy);
				(*this)(y, X - 1 + i)(NodeMap::Syy) = - get(y, X - 1 - i).get(NodeMap::Syy);
			}
		}
	}
	
	if (rank == 0 && borderConditions.at("bottom") == BorderConditions::FreeBorder) {
		for (int x = 0; x < X; x++) {
			for (int i = 1; i <= accuracyOrder; i++) {
				(*this)(-i, x)(NodeMap::Vx)  =   get(i, x).get(NodeMap::Vx);
				(*this)(-i, x)(NodeMap::Vy)  =   get(i, x).get(NodeMap::Vy);
				(*this)(-i, x)(NodeMap::Sxx) = - get(i, x).get(NodeMap::Sxx);
				(*this)(-i, x)(NodeMap::Sxy) = - get(i, x).get(NodeMap::Sxy);
				(*this)(-i, x)(NodeMap::Syy) = - get(i, x).get(NodeMap::Syy);
			}
		}
	}
	if (rank == numberOfWorkers - 1 && borderConditions.at("up") == BorderConditions::FreeBorder) {
		for (int x = 0; x < X; x++) {
			for (int i = 1; i <= accuracyOrder; i++) {
				(*this)(Y - 1 + i, x)(NodeMap::Vx)  =   get(Y - 1 - i, x).get(NodeMap::Vx);
				(*this)(Y - 1 + i, x)(NodeMap::Vy)  =   get(Y - 1 - i, x).get(NodeMap::Vy);
				(*this)(Y - 1 + i, x)(NodeMap::Sxx) = - get(Y - 1 - i, x).get(NodeMap::Sxx);
				(*this)(Y - 1 + i, x)(NodeMap::Sxy) = - get(Y - 1 - i, x).get(NodeMap::Sxy);
				(*this)(Y - 1 + i, x)(NodeMap::Syy) = - get(Y - 1 - i, x).get(NodeMap::Syy);
			}
		}
	}
}


void Mesh::applyInitialConditions() {
	if (initialConditions == InitialConditions::Zero) {
		return;

	} else if (initialConditions == InitialConditions::TestExplosion) {
		if (numberOfWorkers != 1) throw "This condition only for sequence version";
		int xTo = (X % 2 /* X is odd? */) ? X / 2 : X / 2 + 1;
		int xFrom = X / 2;
		int yTo = (Y % 2 /* Y is odd? */) ? Y / 2 : Y / 2 + 1;
		int yFrom = Y / 2;
		for (int x = xFrom; x <= xTo; x++) {
			for (int y = yFrom; y <= yTo; y++) {
				(*this)(y, x)(NodeMap::Sxx) = (*this)(y, x)(NodeMap::Syy) = 1.0;
			}
		}
		return;

	} else if(initialConditions == InitialConditions::Explosion) {
		real R = 0.1 * fmin(X, globalY);
		for (int x = 0; x <= X; x++) {
			for (int y = 0; y <= Y; y++) {
				if ( (x - X / 2)*(x - X / 2) + (y + startY - globalY / 2)*(y + startY - globalY / 2) <= R*R )
					(*this)(y, x)(NodeMap::Sxx) = (*this)(y, x)(NodeMap::Syy)  = - 1.0;
			}
		}
		return;

	} else if (initialConditions == InitialConditions::PWaveX) {
		for (int x = 2; x < 0.15 * X + 2; x++) {
			for (int y = 0; y < Y; y++) {
				(*this)(y, x).u = defaultMatrix->A(0).U1.getColumn(1);
			}
		}

	} else if (initialConditions == InitialConditions::PWaveY) {
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				if (y + startY >= 2 && y + startY < 0.45 * globalY + 2)
					(*this)(y, x).u = defaultMatrix->A(1).U1.getColumn(1);
			}
		}

	} else if (initialConditions == InitialConditions::SWaveX) {
		for (int x = 2; x < 0.15 * X + 2; x++) {
			for (int y = 0; y < Y; y++) {
				(*this)(y, x).u = defaultMatrix->A(0).U1.getColumn(3);
			}
		}

	} else if (initialConditions == InitialConditions::SWaveY) {
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				if (y + startY >= 2 && y + startY < 0.15 * globalY + 2)
					(*this)(y, x).u = defaultMatrix->A(1).U1.getColumn(3);
			}
		}

	} else if (initialConditions == InitialConditions::SWaveXBackward) {
		for (int x = (int) (0.85 * X - 2); x < X - 2; x++) {
			for (int y = 0; y < Y; y++) {
				(*this)(y, x).u = defaultMatrix->A(0).U1.getColumn(2);
			}
		}

	} else if (initialConditions == InitialConditions::SxxOnly) {
		if (numberOfWorkers != 1) throw "This condition only for sequence version";
		(*this)(Y / 2, X / 2)(NodeMap::Sxx) = 5.5;

	} else if (initialConditions == InitialConditions::PWaveXBackward) {
		for (int x = (int)(0.15 * X); x < 0.3 * X; x++) {
			for (int y = 0; y < Y; y++) {
				(*this)(y, x).u = defaultMatrix->A(0).U1.getColumn(0);
			}
		}

	} else if (initialConditions == InitialConditions::PWaveYBackward) {
		for (int y = (int)(0.15 * Y); y < 0.3 * Y; y++) {
			for (int x = 0; x < X; x++) {
				(*this)(y, x).u = defaultMatrix->A(1).U1.getColumn(1);
			}
		}

	} else if (initialConditions == InitialConditions::ExplosionAtTheLeft) {
		real R = 0.05 * fmin(X, globalY);
		for (int y = 0; y < Y; y++) {
			for (int x = 0; x < X; x++) {
				if ( (x - X / 6)*(x - X / 6) + (y + startY - globalY / 2)*(y + startY - globalY / 2) <= R*R )
					(*this)(y, x)(NodeMap::Sxx) = (*this)(y, x)(NodeMap::Syy)  = - 1.0;
			}
		}

	} else {
		throw("Unknown initial condition");
	}
}
