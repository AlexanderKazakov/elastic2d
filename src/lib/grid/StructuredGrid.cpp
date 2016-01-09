#include <fstream>
#include <algorithm>
#include <limits>

#include "lib/model/IdealElastic1DModel.hpp"
#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/model/IdealElastic3DModel.hpp"
#include "lib/grid/StructuredGrid.hpp"

using namespace gcm;

template<class TModel>
void StructuredGrid<TModel>::initializeImpl(const Task &task) {
	LOG_INFO("Start initialization");

	/* ------------------ Properties and conditions ------------------ */

	accuracyOrder = task.accuracyOrder; // order of accuracy of spatial interpolation
	assert_gt(accuracyOrder, 0);

	X = task.X; // number of nodes along x direction
	assert_gt(X, 0);
	Z = task.Z; // number of nodes along z direction
	assert_gt(Z, 0);

	// we divide the grid among processes equally along y-axis
	int numberOfNodesAlongYPerOneCore = (int) std::round((real) task.Y / numberOfWorkers);
	assert_gt(numberOfNodesAlongYPerOneCore, 0);
	Y = numberOfNodesAlongYPerOneCore; // number of nodes along y direction
	if (rank == numberOfWorkers - 1) Y = task.Y - numberOfNodesAlongYPerOneCore * (numberOfWorkers - 1);
	assert_gt(Y, 0);
	assert_ge(Y, 2 * accuracyOrder);

	globalY = task.Y;
	assert_gt(globalY, 0);

	h[0] = task.xLength / (X - 1);
	h[1] = task.yLength / (task.Y - 1);
	h[2] = task.zLength / (Z - 1);

	if (X == 1) h[0] = std::numeric_limits<real>::max();
	if (Y == 1) h[1] = std::numeric_limits<real>::max();
	if (Z == 1) h[2] = std::numeric_limits<real>::max();

	for (int j = 0; j < 3; j++) {
		assert_gt(h[j], 0.0);
		assert_eq(h[j], h[j]);
	}

	/* ------------------ Properties and conditions (end) ------------------ */

	startY = rank * numberOfNodesAlongYPerOneCore;

	nodes = new Node[(Y + 2 * accuracyOrder) * (X + 2 * accuracyOrder)];
	assert_true(nodes);

	for (int i = 0; i < (Y + 2 * accuracyOrder) * (X + 2 * accuracyOrder); ++i) {
		linal::clear(nodes[i]);
		nodes[i].matrix = nullptr;
	}

	defaultMatrix = std::make_shared<typename TModel::GcmMatrices>(task.rho0, task.lambda0, task.mu0);
	assert_true(defaultMatrix);
	maximalLambda = defaultMatrix->getMaximalEigenvalue();
	for (int y = 0; y < Y; ++y) {
		for (int x = 0; x < X; ++x) {
			(*this)(y, x).matrix = defaultMatrix;
		}
	}

	applyInitialConditions();
}

template<class TModel>
typename StructuredGrid<TModel>::Matrix StructuredGrid<TModel>::interpolateValuesAround
		(const int stage, const int y, const int x, const Vector& dx) const {

	Matrix ans;
	std::vector<Vector> src(accuracyOrder + 1);
	Vector res;
	for (int k = 0; k < Node::M; k++) {
		findSourcesForInterpolation(stage, y, x, dx(k), src);
		interpolator.minMaxInterpolate(res, src, fabs(dx(k)) / h[stage]);
		ans.setColumn(k, res);
	}

	return ans;
}

template<class TModel>
void StructuredGrid<TModel>::findSourcesForInterpolation(const int stage, const int y, const int x,
                                                         const real &dx, std::vector<Vector>& src) const {
	assert_ge(x, 0);
	assert_lt(x, X);
	assert_ge(y, 0);
	assert_lt(y, Y);

	const int alongX = (stage == 0) * ( (dx > 0) ? 1 : -1 );
	const int alongY = (stage == 1) * ( (dx > 0) ? 1 : -1 );
	for (int k = 0; k < src.size(); k++) {
		src[k] = get(y + alongY * k, x + alongX * k);
	}
}

template<class TModel>
void StructuredGrid<TModel>::changeRheology(const real& rho2rho0, const real& lambda2lambda0, const real& mu2mu0) {

	auto oldMatrix = defaultMatrix;
	auto newRheologyMatrix = std::make_shared<typename TModel::GcmMatrices>(rho2rho0 * oldMatrix->rho,
	                                                       lambda2lambda0 * oldMatrix->lambda,
	                                                       mu2mu0 * oldMatrix->mu);
	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++) {
			if ((y + startY) * h[1] >= 0.5) {
				(*this)(y, x).matrix = newRheologyMatrix;
			}
		}
	}

	maximalLambda = fmax(maximalLambda, newRheologyMatrix->getMaximalEigenvalue());
}

template<class TModel>
void StructuredGrid<TModel>::applyBorderConditions() {
	if (borderConditions.at("left") == BorderConditions::FreeBorder) {
		for (int y = 0; y < Y; y++) {
			for (int i = 1; i <= accuracyOrder; i++) {
				for (int j = 0; j < TModel::Node::V_SIZE; j++) {
					(*this)(y, - i).V[j] = get(y, i).V[j];
				}
				for (int j = 0; j < TModel::Node::S_SIZE; j++) {
					(*this)(y, - i).S[j] = - get(y, i).S[j];
				}
			}
		}
	}
	if (borderConditions.at("right") == BorderConditions::FreeBorder) {
		for (int y = 0; y < Y; y++) {
			for (int i = 1; i <= accuracyOrder; i++) {
				for (int j = 0; j < TModel::Node::V_SIZE; j++) {
					(*this)(y, X - 1 + i).V[j] = get(y, X - 1 - i).V[j];
				}
				for (int j = 0; j < TModel::Node::S_SIZE; j++) {
					(*this)(y, X - 1 + i).S[j] = - get(y, X - 1 - i).S[j];
				}
			}
		}
	}
	
	if (rank == 0 && borderConditions.at("bottom") == BorderConditions::FreeBorder) {
		for (int x = 0; x < X; x++) {
			for (int i = 1; i <= accuracyOrder; i++) {
				for (int j = 0; j < TModel::Node::V_SIZE; j++) {
					(*this)( - i, x).V[j] = get(i, x).V[j];
				}
				for (int j = 0; j < TModel::Node::S_SIZE; j++) {
					(*this)( - i, x).S[j] = - get(i, x).S[j];
				}
			}
		}
	}
	if (rank == numberOfWorkers - 1 && borderConditions.at("up") == BorderConditions::FreeBorder) {
		for (int x = 0; x < X; x++) {
			for (int i = 1; i <= accuracyOrder; i++) {
				for (int j = 0; j < TModel::Node::V_SIZE; j++) {
					(*this)(Y - 1 + i, x).V[j] = get(Y - 1 - i, x).V[j];
				}
				for (int j = 0; j < TModel::Node::S_SIZE; j++) {
					(*this)(Y - 1 + i, x).S[j] = - get(Y - 1 - i, x).S[j];
				}
			}
		}
	}
}

template<class TModel>
void StructuredGrid<TModel>::applyInitialConditions() {
	if (initialConditions == InitialConditions::Zero) {
		return;

	} else if (initialConditions == InitialConditions::TestExplosion) {
		if (numberOfWorkers != 1) THROW_INVALID_ARG("This condition only for sequence version");
		int xTo = (X % 2 /* X is odd? */) ? X / 2 : X / 2 + 1;
		int xFrom = X / 2;
		int yTo = (Y % 2 /* Y is odd? */) ? Y / 2 : Y / 2 + 1;
		int yFrom = Y / 2;
		for (int x = xFrom; x <= xTo; x++) {
			for (int y = yFrom; y <= yTo; y++) {
				(*this)(y, x).setPressure(-1.0);
			}
		}
		return;

	} else if(initialConditions == InitialConditions::Explosion) {
		real R = 0.1 * fmin(X, globalY);
		for (int x = 0; x <= X; x++) {
			for (int y = 0; y <= Y; y++) {
				if ( (x - X / 2)*(x - X / 2) + (y + startY - globalY / 2)*(y + startY - globalY / 2) <= R*R )
					(*this)(y, x).setPressure(1.0);
			}
		}
		return;

	} else if (initialConditions == InitialConditions::PWaveX) {
		for (int x = 2; x < 0.15 * X + 2; x++) {
			for (int y = 0; y < Y; y++) {
				(*this)(y, x) = defaultMatrix->A(0).U1.getColumn(1);
			}
		}

	} else if (initialConditions == InitialConditions::PWaveY) {
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				if (y + startY >= 2 && y + startY < 0.45 * globalY + 2)
					(*this)(y, x) = defaultMatrix->A(1).U1.getColumn(1);
			}
		}

	} else if (initialConditions == InitialConditions::SWaveX) {
		for (int x = 2; x < 0.15 * X + 2; x++) {
			for (int y = 0; y < Y; y++) {
				(*this)(y, x) = defaultMatrix->A(0).U1.getColumn(3);
			}
		}

	} else if (initialConditions == InitialConditions::SWaveY) {
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				if (y + startY >= 2 && y + startY < 0.15 * globalY + 2)
					(*this)(y, x) = defaultMatrix->A(1).U1.getColumn(3);
			}
		}

	} else if (initialConditions == InitialConditions::SWaveXBackward) {
		for (int x = (int) (0.85 * X - 2); x < X - 2; x++) {
			for (int y = 0; y < Y; y++) {
				(*this)(y, x) = defaultMatrix->A(0).U1.getColumn(2);
			}
		}

	} else if (initialConditions == InitialConditions::SxxOnly) {
		if (numberOfWorkers != 1) THROW_INVALID_ARG("This condition only for sequence version");
		(*this)(Y / 2, X / 2).Sxx = 5.5;

	} else if (initialConditions == InitialConditions::PWaveXBackward) {
		for (int x = (int)(0.15 * X); x < 0.3 * X; x++) {
			for (int y = 0; y < Y; y++) {
				(*this)(y, x) = defaultMatrix->A(0).U1.getColumn(0);
			}
		}

	} else if (initialConditions == InitialConditions::PWaveYBackward) {
		for (int y = (int)(0.15 * Y); y < 0.3 * Y; y++) {
			for (int x = 0; x < X; x++) {
				(*this)(y, x) = defaultMatrix->A(1).U1.getColumn(1);
			}
		}

	} else if (initialConditions == InitialConditions::ExplosionAtTheLeft) {
		real R = 0.05 * fmin(X, globalY);
		for (int y = 0; y < Y; y++) {
			for (int x = 0; x < X; x++) {
				if ( (x - X / 6)*(x - X / 6) + (y + startY - globalY / 2)*(y + startY - globalY / 2) <= R*R )
					(*this)(y, x).setPressure(1.0);
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
