#include <fstream>
#include <algorithm>

#include "lib/grid/StructuredGrid.hpp"

#include "lib/task/InitialCondition.hpp"
#include "lib/model/IdealElastic1DModel.hpp"
#include "lib/model/IdealElastic2DModel.hpp"
#include "lib/model/IdealElastic3DModel.hpp"

using namespace gcm;

template<class TModel>
void StructuredGrid<TModel>::initializeImpl(const Task &task) {
	LOG_INFO("Start initialization");

	accuracyOrder = task.accuracyOrder; // order of accuracy of spatial interpolation
	assert_ge(accuracyOrder, 1);

	globalX = task.X; // number of nodes along x direction of all meshes (from all cores)
	assert_ge(globalX, 2 * accuracyOrder);
	Y = task.Y; // number of nodes along y direction
	if (Y != 1) assert_ge(Y, 2 * accuracyOrder);
	Z = task.Z; // number of nodes along z direction
	if (Z != 1) assert_ge(Z, 2 * accuracyOrder);

	// we divide the grid among processes equally along x-axis
	int numberOfNodesAlongXPerOneCore = (int) std::round((real) task.X / this->numberOfWorkers);
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
	assert_eq(startX, startX); assert_eq(startY, startY); assert_eq(startZ, startZ); // this is supposed to catch NaN

	borderConditions = task.borderConditions;

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
		src[k] = get(x + alongX * k, y + alongY * k, z + alongZ * k).u;
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

	if (this->rank == 0 && borderConditions.at(CUBIC_BORDERS::X_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int y = 0; y < Y; y++) {
			for (int z = 0; z < Z; z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::V_SIZE; j++) {
						(*this)( - i, y, z).u.V[j] = get(i, y, z).u.V[j];
					}
					for (int j = 0; j < Vector::S_SIZE; j++) {
						(*this)( - i, y, z).u.S[j] = - get(i, y, z).u.S[j];
					}
				}
			}
		}
	}
	if (this->rank == this->numberOfWorkers - 1 && borderConditions.at(CUBIC_BORDERS::X_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int y = 0; y < Y; y++) {
			for (int z = 0; z < Z; z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::V_SIZE; j++) {
						(*this)(X - 1 + i, y, z).u.V[j] = get(X - 1 - i, y, z).u.V[j];
					}
					for (int j = 0; j < Vector::S_SIZE; j++) {
						(*this)(X - 1 + i, y, z).u.S[j] = - get(X - 1 - i, y, z).u.S[j];
					}
				}
			}
		}
	}
	
	if (Y != 1 && borderConditions.at(CUBIC_BORDERS::Y_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < X; x++) {
			for (int z = 0; z < Z; z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::V_SIZE; j++) {
						(*this)(x, - i, z).u.V[j] = get(x, i, z).u.V[j];
					}
					for (int j = 0; j < Vector::S_SIZE; j++) {
						(*this)(x, - i, z).u.S[j] = - get(x, i, z).u.S[j];
					}
				}
			}
		}
	}
	if (Y != 1 && borderConditions.at(CUBIC_BORDERS::Y_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < X; x++) {
			for (int z = 0; z < Z; z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::V_SIZE; j++) {
						(*this)(x, Y - 1 + i, z).u.V[j] = get(x, Y - 1 - i, z).u.V[j];
					}
					for (int j = 0; j < Vector::S_SIZE; j++) {
						(*this)(x, Y - 1 + i, z).u.S[j] = - get(x, Y - 1 - i, z).u.S[j];
					}
				}
			}
		}
	}

	if (Z != 1 && borderConditions.at(CUBIC_BORDERS::Z_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::V_SIZE; j++) {
						(*this)(x, y, - i).u.V[j] = get(x, y, i).u.V[j];
					}
					for (int j = 0; j < Vector::S_SIZE; j++) {
						(*this)(x, y, - i).u.S[j] = - get(x, y, i).u.S[j];
					}
				}
			}
		}
	}
	if (Z != 1 && borderConditions.at(CUBIC_BORDERS::Z_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::V_SIZE; j++) {
						(*this)(x, y, Z - 1 + i).u.V[j] = get(x, y, Z - 1 - i).u.V[j];
					}
					for (int j = 0; j < Vector::S_SIZE; j++) {
						(*this)(x, y, Z - 1 + i).u.S[j] = - get(x, y, Z - 1 - i).u.S[j];
					}
				}
			}
		}
	}
}

template<class TModel>
void StructuredGrid<TModel>::applyInitialConditions(const Task& task) {
	InitialCondition<TModel> initialCondition;
	initialCondition.initialize(task, this->defaultMatrix);

	for (int x = 0; x < X; x++) {
		for (int y = 0; y < Y; y++) {
			for (int z = 0; z < Z; z++) {
				initialCondition.apply((*this)(x, y, z), getCoordinates(x, y, z));
			}
		}
	}
}

template<class TModel>
real StructuredGrid<TModel>::getMinimalSpatialStep() const {
	return fmin(h[0], fmin(h[1], h[2]));
}


template class StructuredGrid<IdealElastic1DModel>;
template class StructuredGrid<IdealElastic2DModel>;
template class StructuredGrid<IdealElastic3DModel>;
