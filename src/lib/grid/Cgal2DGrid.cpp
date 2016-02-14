#include <fstream>
#include <algorithm>

#include <lib/grid/Cgal2DGrid.hpp>
#include <lib/util/task/InitialCondition.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<typename TModel>
void Cgal2DGrid<TModel>::initializeImpl(const Task &task) {
	LOG_INFO("Start initialization");

	// Generate the mesh
	mesh.initialize(task);

	auto initialSize = mesh.verticesSize();
	this->pdeVectors.zeroInitialize(initialSize);
	this->pdeVectorsNew.zeroInitialize(initialSize);
	this->odeValues.zeroInitialize(initialSize);
	this->gcmMatrices.zeroInitialize(initialSize);

	typename TModel::Material material;
	material.initialize(task);
	auto gcmMatricesPtr = new GCM_MATRICES(material);
	for (auto& matrix : this->gcmMatrices) {
		matrix = gcmMatricesPtr;
	}
	maximalLambda = gcmMatricesPtr->getMaximalEigenvalue();
	minimalSpatialStep = mesh.getMinimalSpatialStep();
}

template<typename TModel>
void Cgal2DGrid<TModel>::applyInitialConditions(const Task& task) {
	InitialCondition<TModel> initialCondition;
	initialCondition.initialize(task);

	for (size_t i = 0; i < mesh.verticesSize(); i++) {
		initialCondition.apply((*this)(i), getCoordinates(i));
	}
}

template<typename TModel>
typename Cgal2DGrid<TModel>::Matrix Cgal2DGrid<TModel>::interpolateValuesAround
		(const int stage, const size_t i, const PdeVector& dx) const {

	Matrix ans;
	std::vector<PdeVector> src(1);
	PdeVector res;
	for (int k = 0; k < PdeVector::M; k++) {
		findSourcesForInterpolation(stage, i, dx(k), src);
		/* TODO interpolator.minMaxInterpolate(res, src, fabs(dx(k)) / h(stage));*/
		ans.setColumn(k, res);
	}

	return ans;
}

template<typename TModel>
void Cgal2DGrid<TModel>::findSourcesForInterpolation
		(const int stage, const size_t i,
		 const real &dx, std::vector<PdeVector>& src) const {

	for (int k = 0; k < src.size(); k++) {
		// TODO
		src[(unsigned long)k] = get(i + stage + (size_t)dx);
	}
}

template<typename TModel>
void Cgal2DGrid<TModel>::beforeStageImpl() {
	exchangeNodesWithNeighbors();
	applyBorderConditions();
}

template<typename TModel>
void Cgal2DGrid<TModel>::exchangeNodesWithNeighbors() {
	// TODO
}

template<typename TModel>
void Cgal2DGrid<TModel>::applyBorderConditions() {
	// TODO
}


template class Cgal2DGrid<Elastic2DModel>;
template class Cgal2DGrid<ContinualDamageElastic2DModel>;
template class Cgal2DGrid<IdealPlastic2DModel>;
