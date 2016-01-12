#include <fstream>
#include <algorithm>

#include "lib/model/IdealElastic1DModel.hpp"
#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/model/IdealElastic3DModel.hpp"
#include "lib/grid/StructuredGrid.hpp"

using namespace gcm;

template<class TModel>
void StructuredGrid<TModel>::initializeImpl(const Task &task) {
	LOG_INFO("Start initialization");

	accuracyOrder = task.accuracyOrder; // order of accuracy of spatial interpolation
	assert_gt(accuracyOrder, 0);

	globalX = task.X; // number of nodes along x direction of all meshes (from all cores)
	assert_gt(globalX, 0);
	Y = task.Y; // number of nodes along y direction
	assert_gt(Y, 0);
	Z = task.Z; // number of nodes along z direction
	assert_gt(Z, 0);

	// we divide the grid among processes equally along x-axis
	int numberOfNodesAlongXPerOneCore = (int) std::round((real) task.X / this->numberOfWorkers);
	assert_gt(numberOfNodesAlongXPerOneCore, 0);
	X = numberOfNodesAlongXPerOneCore; // number of nodes along x direction on this mesh
	// in order to keep specified in task number of nodes
	if (this->rank == this->numberOfWorkers - 1) X = task.X - numberOfNodesAlongXPerOneCore * (this->numberOfWorkers - 1);
	assert_ge(X, 2 * accuracyOrder);

	h[0] = task.xLength / (task.X - 1);
	h[1] = task.yLength / (task.Y - 1);
	h[2] = task.zLength / (task.Z - 1);

	if (task.X == 1) h[0] = std::numeric_limits<real>::max();
	if (task.Y == 1) h[1] = std::numeric_limits<real>::max();
	if (task.Z == 1) h[2] = std::numeric_limits<real>::max();

	for (int j = 0; j < 3; j++) {
		assert_gt(h[j], 0.0);
		assert_eq(h[j], h[j]); // this is supposed to catch NaN
	}

	globalStartXindex = this->rank * numberOfNodesAlongXPerOneCore;

	startX = task.startX + globalStartXindex * h[0];
	startY = task.startY;
	startZ = task.startZ;
	assert_eq(startX, startX); assert_eq(startY, startY); assert_eq(startZ, startZ);

	this->nodes.resize( (unsigned long) (X + 2 * accuracyOrder) * (Y + 2 * accuracyOrder) * (Z + 2 * accuracyOrder) );
}

template<class TModel>
typename StructuredGrid<TModel>::Matrix StructuredGrid<TModel>::interpolateValuesAround
		(const int stage, const int x, const int y, const int z, const Vector& dx) const {

	Matrix ans;
	std::vector<Vector> src(accuracyOrder + 1);
	Vector res;
	for (int k = 0; k < Node::M; k++) {
		findSourcesForInterpolation(stage, x, y, z, dx(k), src);
		interpolator.minMaxInterpolate(res, src, fabs(dx(k)) / h[stage]);
		ans.setColumn(k, res);
	}

	return ans;
}

template<class TModel>
void StructuredGrid<TModel>::findSourcesForInterpolation(const int stage, const int x, const int y, const int z,
                                                         const real &dx, std::vector<Vector>& src) const {

	const int alongX = (stage == 0) * ( (dx > 0) ? 1 : -1 );
	const int alongY = (stage == 1) * ( (dx > 0) ? 1 : -1 );
	const int alongZ = (stage == 2) * ( (dx > 0) ? 1 : -1 );
	for (int k = 0; k < src.size(); k++) {
		src[k] = get(x + alongX * k, y + alongY * k, z + alongZ * k);
	}
}

template<class TModel>
void StructuredGrid<TModel>::changeRheology(const real& rho2rho0, const real& lambda2lambda0, const real& mu2mu0) {

	auto oldMatrix = this->defaultMatrix;
	auto newRheologyMatrix = std::make_shared<typename TModel::GcmMatrices>(rho2rho0 * oldMatrix->rho,
	                                                       lambda2lambda0 * oldMatrix->lambda,
	                                                       mu2mu0 * oldMatrix->mu);
	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++) {
			for (int z = 0; z < Z; z++) {
				if (y * h[1] >= 0.5) {
					(*this)(x, y, z).matrix = newRheologyMatrix;
				}
			}
		}
	}

	this->maximalLambda = fmax(this->maximalLambda, newRheologyMatrix->getMaximalEigenvalue());
}

template<class TModel>
void StructuredGrid<TModel>::applyBorderConditions() {

	if (this->rank == 0 && this->borderConditions.at(Border::X_LEFT) == BorderConditions::FreeBorder) {
		for (int y = 0; y < Y; y++) {
			for (int z = 0; z < Z; z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < TModel::Node::V_SIZE; j++) {
						(*this)( - i, y, z).V[j] = get(i, y, z).V[j];
					}
					for (int j = 0; j < TModel::Node::S_SIZE; j++) {
						(*this)( - i, y, z).S[j] = - get(i, y, z).S[j];
					}
				}
			}
		}
	}
	if (this->rank == this->numberOfWorkers - 1 && this->borderConditions.at(Border::X_RIGHT) == BorderConditions::FreeBorder) {
		for (int y = 0; y < Y; y++) {
			for (int z = 0; z < Z; z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < TModel::Node::V_SIZE; j++) {
						(*this)(X - 1 + i, y, z).V[j] = get(X - 1 - i, y, z).V[j];
					}
					for (int j = 0; j < TModel::Node::S_SIZE; j++) {
						(*this)(X - 1 + i, y, z).S[j] = - get(X - 1 - i, y, z).S[j];
					}
				}
			}
		}
	}
	
	if (this->borderConditions.at(Border::Y_LEFT) == BorderConditions::FreeBorder) {
		for (int x = 0; x < X; x++) {
			for (int z = 0; z < Z; z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < TModel::Node::V_SIZE; j++) {
						(*this)(x, - i, z).V[j] = get(x, i, z).V[j];
					}
					for (int j = 0; j < TModel::Node::S_SIZE; j++) {
						(*this)(x, - i, z).S[j] = - get(x, i, z).S[j];
					}
				}
			}
		}
	}
	if (this->borderConditions.at(Border::Y_RIGHT) == BorderConditions::FreeBorder) {
		for (int x = 0; x < X; x++) {
			for (int z = 0; z < Z; z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < TModel::Node::V_SIZE; j++) {
						(*this)(x, Y - 1 + i, z).V[j] = get(x, Y - 1 - i, z).V[j];
					}
					for (int j = 0; j < TModel::Node::S_SIZE; j++) {
						(*this)(x, Y - 1 + i, z).S[j] = - get(x, Y - 1 - i, z).S[j];
					}
				}
			}
		}
	}

	if (this->borderConditions.at(Border::Z_LEFT) == BorderConditions::FreeBorder) {
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < TModel::Node::V_SIZE; j++) {
						(*this)(x, y, - i).V[j] = get(x, y, i).V[j];
					}
					for (int j = 0; j < TModel::Node::S_SIZE; j++) {
						(*this)(x, y, - i).S[j] = - get(x, y, i).S[j];
					}
				}
			}
		}
	}
	if (this->borderConditions.at(Border::Z_RIGHT) == BorderConditions::FreeBorder) {
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < TModel::Node::V_SIZE; j++) {
						(*this)(x, y, Z - 1 + i).V[j] = get(x, y, Z - 1 - i).V[j];
					}
					for (int j = 0; j < TModel::Node::S_SIZE; j++) {
						(*this)(x, y, Z - 1 + i).S[j] = - get(x, y, Z - 1 - i).S[j];
					}
				}
			}
		}
	}
}

template<class TModel>
void StructuredGrid<TModel>::applyInitialConditions() {
	if (this->initialConditions == InitialConditions::Zero) {
		return;

	} else if (this->initialConditions == InitialConditions::TestExplosion) {
		if (this->numberOfWorkers != 1) THROW_INVALID_ARG("This condition only for sequence version");
		int xTo = (X % 2 /* X is odd? */) ? X / 2 : X / 2 + 1;
		int xFrom = X / 2;
		int yTo = (Y % 2 /* Y is odd? */) ? Y / 2 : Y / 2 + 1;
		int yFrom = Y / 2;
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				for (int z = 0; z < Z; z++) {
					if (x >= xFrom && x <= xTo && y >= yFrom && y <= yTo) {
						(*this)(x, y, z).setPressure(-1.0);
					}
				}
			}
		}
		return;

	} else if(this->initialConditions == InitialConditions::Explosion) {
		real R = 0.1 * (globalX + Y + Z) / 3;
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				for (int z = 0; z < Z; z++) {
					if ( (y - Y / 2) * (y - Y / 2) + (z - Z / 2) * (z - Z / 2) + (x + globalStartXindex - globalX / 2) * (x +
					                                                                                                      globalStartXindex - globalX / 2) <= R * R )
						(*this)(x, y, z).setPressure(1.0);
				}
			}
		}
		return;

	} else if (this->initialConditions == InitialConditions::PWaveX) {
		for (int x = 2; x < 0.15 * X + 2; x++) {
			for (int y = 0; y < Y; y++) {
				for (int z = 0; z < Z; z++) {
					(*this)(x, y, z) = this->defaultMatrix->A(0).U1.getColumn(1);
				}
			}
		}

	} else if (this->initialConditions == InitialConditions::PWaveY) {
		for (int x = 0; x < X; x++) {
			for (int y = 2; y < 0.45 * Y + 2; y++) {
				for (int z = 0; z < Z; z++) {
						(*this)(x, y, z) = this->defaultMatrix->A(1).U1.getColumn(1);
				}
			}
		}

	} else if (this->initialConditions == InitialConditions::SWaveX) {
		for (int x = 2; x < 0.15 * X + 2; x++) {
			for (int y = 0; y < Y; y++) {
				for (int z = 0; z < Z; z++) {
					(*this)(x, y, z) = this->defaultMatrix->A(0).U1.getColumn(3);
				}
			}
		}

	} else if (this->initialConditions == InitialConditions::SWaveY) {
		for (int x = 0; x < X; x++) {
			for (int y = 2; y < 0.15 * Y + 2; y++) {
				for (int z = 0; z < Z; z++) {
						(*this)(x, y, z) = this->defaultMatrix->A(1).U1.getColumn(3);
				}
			}
		}

	} else if (this->initialConditions == InitialConditions::SWaveXBackward) {
		for (int x = (int) (0.85 * X - 2); x < X - 2; x++) {
			for (int y = 0; y < Y; y++) {
				for (int z = 0; z < Z; z++) {
					(*this)(x, y, z) = this->defaultMatrix->A(0).U1.getColumn(2);
				}
			}
		}

	} else if (this->initialConditions == InitialConditions::SxxOnly) {
		if (this->numberOfWorkers != 1) THROW_INVALID_ARG("This condition only for sequence version");
		(*this)(X / 2, Y / 2, Z / 2).Sxx = 5.5;

	} else if (this->initialConditions == InitialConditions::PWaveXBackward) {
		for (int x = (int)(0.15 * X); x < 0.3 * X; x++) {
			for (int y = 0; y < Y; y++) {
				for (int z = 0; z < Z; z++) {
					(*this)(x, y, z) = this->defaultMatrix->A(0).U1.getColumn(0);
				}
			}
		}

	} else if (this->initialConditions == InitialConditions::PWaveYBackward) {
		for (int x = 0; x < X; x++) {
			for (int y = (int)(0.15 * Y); y < 0.3 * Y; y++) {
				for (int z = 0; z < Z; z++) {
					(*this)(x, y, z) = this->defaultMatrix->A(1).U1.getColumn(1);
				}
			}
		}

	} else {
		THROW_INVALID_ARG("Unknown initial condition");
	}
}

template<class TModel>
real StructuredGrid<TModel>::getMinimalSpatialStep() const {
	return fmin(h[0], fmin(h[1], h[2]));
}


template class StructuredGrid<IdealElastic1DModel>;
template class StructuredGrid<IdealElastic2DModel>;
template class StructuredGrid<IdealElastic3DModel>;
