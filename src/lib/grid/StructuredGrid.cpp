#include <fstream>
#include <algorithm>

#include <lib/grid/StructuredGrid.hpp>
#include <lib/util/task/InitialCondition.hpp>

using namespace gcm;

template<class TNode>
void StructuredGrid<TNode>::initializeImpl(const Task &task) {
	LOG_INFO("Start initialization");
	accuracyOrder = task.accuracyOrder; // order of accuracy of spatial interpolation
	assert_ge(accuracyOrder, 1);

	h = linal::plainDivision(task.lengthes, task.sizes - linal::VectorInt<3>({1, 1, 1}));
	sizes = task.sizes;
	startR = task.startR;

	// MPI - we divide the grid among processes equally along x-axis
	globalX = task.sizes(0); // number of nodes along x direction of all meshes (from all cores)
	int numberOfNodesAlongXPerOneCore = (int) std::round((real) task.sizes(0) / this->numberOfWorkers);
	sizes(0) = numberOfNodesAlongXPerOneCore; // number of nodes along x direction on this mesh
	// in order to keep specified in task number of nodes
	if (this->rank == this->numberOfWorkers - 1) sizes(0) = task.sizes(0) - numberOfNodesAlongXPerOneCore * (this->numberOfWorkers - 1);
	startR(0) += (this->rank * numberOfNodesAlongXPerOneCore) * h(0); // mpi parallel along X axis
	globalStartXindex = this->rank * numberOfNodesAlongXPerOneCore; // for testing
	// MPI - (end)

	for (int j = 0; j < 3; j++) {
		if (sizes(j) != 1) assert_ge(sizes(j), 2 * accuracyOrder);
		assert_gt(h(j), 0.0);
		assert_eq(h(j), h(j)); // this is supposed to catch NaN
		assert_eq(startR(j), startR(j));
	}
	
	borderConditions = task.borderConditions;

	this->nodes.resize( (unsigned long) (sizes(0) + 2 * accuracyOrder) * (sizes(1) + 2 * accuracyOrder) * (sizes(2) + 2 * accuracyOrder) );
}

template<class TNode>
typename StructuredGrid<TNode>::Matrix StructuredGrid<TNode>::interpolateValuesAround
		(const int stage, const int x, const int y, const int z, const Vector& dx) const {

	Matrix ans;
	std::vector<Vector> src( (unsigned long) (accuracyOrder + 1) );
	Vector res;
	for (int k = 0; k < Node::M; k++) {
		findSourcesForInterpolation(stage, x, y, z, dx(k), src);
		interpolator.minMaxInterpolate(res, src, fabs(dx(k)) / h(stage));
		ans.setColumn(k, res);
	}

	return ans;
}

template<class TNode>
void StructuredGrid<TNode>::findSourcesForInterpolation(const int stage, const int x, const int y, const int z,
                                                         const real &dx, std::vector<Vector>& src) const {

	const int alongX = (stage == 0) * ( (dx > 0) ? 1 : -1 );
	const int alongY = (stage == 1) * ( (dx > 0) ? 1 : -1 );
	const int alongZ = (stage == 2) * ( (dx > 0) ? 1 : -1 );
	for (int k = 0; k < src.size(); k++) {
		src[(unsigned long)k] = get(x + alongX * k, y + alongY * k, z + alongZ * k).u;
	}
}

template<class TNode>
void StructuredGrid<TNode>::changeRheology(const real& rho2rho0, const real& lambda2lambda0, const real& mu2mu0) {

	IsotropicMaterial oldMaterial = (*this)(0, 0, 0).matrix->getMaterial();
	IsotropicMaterial newMaterial(rho2rho0 * oldMaterial.rho, lambda2lambda0 * oldMaterial.lambda, mu2mu0 * oldMaterial.mu);
	auto newRheologyMatrix = std::make_shared<typename TNode::GcmMatrices>(newMaterial);

	for (int x = 0; x < sizes(0); x++) {
		for (int y = 0; y < sizes(1); y++) {
			for (int z = 0; z < sizes(2); z++) {
				if (y * h(1) >= 0.5) {
					(*this)(x, y, z).matrix = newRheologyMatrix;
				}
			}
		}
	}

	this->maximalLambda = fmax(this->maximalLambda, newRheologyMatrix->getMaximalEigenvalue());
}

template<class TNode>
void StructuredGrid<TNode>::applyBorderConditions() {

	if (this->rank == 0 && borderConditions.at(CUBIC_BORDERS::X_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int y = 0; y < sizes(1); y++) {
			for (int z = 0; z < sizes(2); z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::DIMENSIONALITY; j++) {
						(*this)( - i, y, z).u.V[j] = get(i, y, z).u.V[j];
					}
					for (int j = 0; j < ( Vector::DIMENSIONALITY * (Vector::DIMENSIONALITY + 1) ) / 2; j++) {
						(*this)( - i, y, z).u.S[j] = - get(i, y, z).u.S[j];
					}
				}
			}
		}
	}
	if (this->rank == this->numberOfWorkers - 1 && borderConditions.at(CUBIC_BORDERS::X_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int y = 0; y < sizes(1); y++) {
			for (int z = 0; z < sizes(2); z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::DIMENSIONALITY; j++) {
						(*this)(sizes(0) - 1 + i, y, z).u.V[j] = get(sizes(0) - 1 - i, y, z).u.V[j];
					}
					for (int j = 0; j < ( Vector::DIMENSIONALITY * (Vector::DIMENSIONALITY + 1) ) / 2; j++) {
						(*this)(sizes(0) - 1 + i, y, z).u.S[j] = - get(sizes(0) - 1 - i, y, z).u.S[j];
					}
				}
			}
		}
	}
	
	if (sizes(1) != 1 && borderConditions.at(CUBIC_BORDERS::Y_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < sizes(0); x++) {
			for (int z = 0; z < sizes(2); z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::DIMENSIONALITY; j++) {
						(*this)(x, - i, z).u.V[j] = get(x, i, z).u.V[j];
					}
					for (int j = 0; j < ( Vector::DIMENSIONALITY * (Vector::DIMENSIONALITY + 1) ) / 2; j++) {
						(*this)(x, - i, z).u.S[j] = - get(x, i, z).u.S[j];
					}
				}
			}
		}
	}
	if (sizes(1) != 1 && borderConditions.at(CUBIC_BORDERS::Y_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < sizes(0); x++) {
			for (int z = 0; z < sizes(2); z++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::DIMENSIONALITY; j++) {
						(*this)(x, sizes(1) - 1 + i, z).u.V[j] = get(x, sizes(1) - 1 - i, z).u.V[j];
					}
					for (int j = 0; j < ( Vector::DIMENSIONALITY * (Vector::DIMENSIONALITY + 1) ) / 2; j++) {
						(*this)(x, sizes(1) - 1 + i, z).u.S[j] = - get(x, sizes(1) - 1 - i, z).u.S[j];
					}
				}
			}
		}
	}

	if (sizes(2) != 1 && borderConditions.at(CUBIC_BORDERS::Z_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < sizes(0); x++) {
			for (int y = 0; y < sizes(1); y++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::DIMENSIONALITY; j++) {
						(*this)(x, y, - i).u.V[j] = get(x, y, i).u.V[j];
					}
					for (int j = 0; j < ( Vector::DIMENSIONALITY * (Vector::DIMENSIONALITY + 1) ) / 2; j++) {
						(*this)(x, y, - i).u.S[j] = - get(x, y, i).u.S[j];
					}
				}
			}
		}
	}
	if (sizes(2) != 1 && borderConditions.at(CUBIC_BORDERS::Z_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < sizes(0); x++) {
			for (int y = 0; y < sizes(1); y++) {
				for (int i = 1; i <= accuracyOrder; i++) {
					for (int j = 0; j < Vector::DIMENSIONALITY; j++) {
						(*this)(x, y, sizes(2) - 1 + i).u.V[j] = get(x, y, sizes(2) - 1 - i).u.V[j];
					}
					for (int j = 0; j < ( Vector::DIMENSIONALITY * (Vector::DIMENSIONALITY + 1) ) / 2; j++) {
						(*this)(x, y, sizes(2) - 1 + i).u.S[j] = - get(x, y, sizes(2) - 1 - i).u.S[j];
					}
				}
			}
		}
	}
}

template<class TNode>
void StructuredGrid<TNode>::applyInitialConditions(const Task& task) {
	InitialCondition<TNode> initialCondition;
	initialCondition.initialize(task);

	for (int x = 0; x < sizes(0); x++) {
		for (int y = 0; y < sizes(1); y++) {
			for (int z = 0; z < sizes(2); z++) {
				initialCondition.apply((*this)(x, y, z), getCoordinates(x, y, z));
			}
		}
	}
}

template<class TNode>
real StructuredGrid<TNode>::getMinimalSpatialStepImpl() const {
	return fmin(h(0), fmin(h(1), h(2)));
}


template class StructuredGrid<IdealElastic1DNode>;
template class StructuredGrid<IdealElastic2DNode>;
template class StructuredGrid<IdealElastic3DNode>;
