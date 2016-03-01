#include <lib/numeric/border_conditions/BorderConditions.hpp>
#include <lib/rheology/models/Model.hpp>

#include "lib/numeric/gcm/GridCharacteristicMethod.hpp"

using namespace gcm;

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::initialize(const Task& task) {
	sizes = task.sizes;
	startR = task.startR;
	lengths = task.lengthes;
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::beforeStatement(const Statement &statement) {
	for (const auto& bc : statement.borderConditions) {
		for (const auto& q : bc.values) {
			assert_eq(PdeVariables::QUANTITIES.count(q.first), 1);
		}
		conditions.push_back(Condition(bc.area, bc.values));
	}
	for (const auto& fr : statement.fractures) {
		for (const auto& q : fr.values) {
			assert_eq(PdeVariables::QUANTITIES.count(q.first), 1);
		}
		int d = fr.direction;
		int index = (int) (sizes(d) * (fr.coordinate - startR(d)) / lengths(d));
		assert_gt(index, 0); assert_lt(index, sizes(d) - 1);
		fractures.push_back(Fracture(d, index,   - 1, fr.area, fr.values));
		fractures.push_back(Fracture(d, index + 1, 1, fr.area, fr.values));
	}
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::applyBorderBeforeStage
(Mesh* mesh_, const real currentTime_, const real timeStep_, const int stage) {
	// handling borders
	mesh = mesh_; currentTime = currentTime_; timeStep = timeStep_; direction = stage;
	// special for x-axis (because MPI partition along x-axis)
	if (direction == 0) {
		if (mesh->getRank() == 0) {
			onTheRight = false;
			handleSide();
		}
		if (mesh->getRank() == mesh->getNumberOfWorkers() - 1) {
			onTheRight = true;
			handleSide();
		}
		return;
	}
	// for other axes
	onTheRight = false;
	handleSide();
	onTheRight = true;
	handleSide();
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::handleSide() const {
	auto borderIter = mesh->slice(direction, 0);
	if (onTheRight)
		borderIter = mesh->slice(direction, mesh->getSizes()(direction) - 1);
	while (borderIter != borderIter.end()) {
		for (const auto& condition : conditions) {
			if (condition.area->contains(mesh->coords(borderIter))) {
				handleBorderPoint(borderIter, condition.values);
			}
		}
		++borderIter;
	}
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::handleBorderPoint
(const Iterator& borderIter, const Map& values) const {
	
	int innerSign = onTheRight ? -1 : 1;
	for (int a = 1; a <= mesh->getAccuracyOrder(); a++) {
		auto realIter = borderIter; realIter(direction) += innerSign * a;
		auto virtIter = borderIter; virtIter(direction) -= innerSign * a;
	
		mesh->_pde(virtIter) = mesh->pde(realIter);
		for (const auto& q : values) {
			const auto& quantity = q.first;
			const auto& timeDependency = q.second;
			real realValue = PdeVariables::QUANTITIES.at(quantity).Get(mesh->pde(realIter));
			real virtValue = - realValue + 2 * timeDependency(currentTime);
			PdeVariables::QUANTITIES.at(quantity).Set(virtValue, mesh->_pde(virtIter));
		}
	}
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::applyBorderAfterStage
(Mesh* mesh_, const real currentTime_, const real timeStep_, const int stage) {
	// handling inner fractures
	mesh = mesh_; currentTime = currentTime_; timeStep = timeStep_; direction = stage;
	for (const auto& fr : fractures) {
		if (fr.direction == stage) {
			allocateHelpMesh();
			auto sliceIter = mesh->slice(direction, fr.index);
			while (sliceIter != sliceIter.end()) {
				if (fr.area->contains(mesh->coords(sliceIter))) {
					handleFracturePoint(sliceIter, fr.values, fr.normal);
				}
				++sliceIter;
			}
			delete helpMesh;
		}
	}
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::allocateHelpMesh() {
	// warning - this is not complete cubic mesh, just an auxiliary smth
	helpMesh = new Mesh();
	*((CubicGrid*)helpMesh) = *((CubicGrid*)mesh);
	helpMesh->sizes = {1, 1, 1};
	helpMesh->sizes(direction) = mesh->getAccuracyOrder();
	helpMesh->allocate();
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::handleFracturePoint
(const Iterator& iter, const Map& values, const int fracNormal) {
	// copy values to helpMesh
	for (int i = 0; i < 2 * mesh->getAccuracyOrder(); i++) {
		Iterator helpMeshIter = {0, 0, 0}; helpMeshIter(direction) += i;
		Iterator realMeshIter = iter;      realMeshIter(direction) += i * fracNormal;
		helpMesh->_pde(helpMeshIter) = mesh->pde(realMeshIter); // todo - what if matrix depends on ode, etc?
		helpMesh->_matrix(helpMeshIter) = mesh->matrix(realMeshIter);
	}
	// apply border conditions before stage on the helpMesh
	onTheRight = false;
	Mesh* tmpMesh = mesh;
	mesh = helpMesh;
	handleBorderPoint({0, 0, 0}, values);
	mesh = tmpMesh;
	// calculate stage on the helpMesh
	GridCharacteristicMethod<Mesh> gcmMetod; // todo static?
	gcmMetod.stage(direction, timeStep * fracNormal, helpMesh);
	// copy calculated values to real mesh
	for (int i = 0; i < mesh->getAccuracyOrder(); i++) {
		Iterator helpMeshIter = {0, 0, 0}; helpMeshIter(direction) += i;
		Iterator realMeshIter = iter;      realMeshIter(direction) += i * fracNormal;
		mesh->_pdeNew(realMeshIter) = helpMesh->pdeNew(helpMeshIter);
	}
}


template class BorderConditions<Elastic1DModel, CubicGrid>;
template class BorderConditions<Elastic2DModel, CubicGrid>;
template class BorderConditions<Elastic3DModel, CubicGrid>;
template class BorderConditions<OrthotropicElastic3DModel, CubicGrid>;
template class BorderConditions<ContinualDamageElastic2DModel, CubicGrid>;
template class BorderConditions<IdealPlastic2DModel, CubicGrid>;
template class BorderConditions<SuperDuperModel, CubicGrid>;
