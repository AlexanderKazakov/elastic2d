#include <lib/numeric/border_conditions/BorderConditions.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/numeric/gcm/GridCharacteristicMethod.hpp>

using namespace gcm;

template<typename TModel, typename TMaterial>
BorderConditions<TModel, CubicGrid, TMaterial>::BorderConditions(const Task& task) {
	sizes = task.cubicGrid.sizes;
	startR = task.cubicGrid.startR;
	lengths = task.cubicGrid.lengthes;
}

template<typename TModel, typename TMaterial>
void BorderConditions<TModel, CubicGrid, TMaterial>::beforeStatement(const Statement &statement) {
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
		innerSurfaces.push_back(InnerSurface(d, index,   - 1, fr.area, fr.values));
		innerSurfaces.push_back(InnerSurface(d, index + 1, 1, fr.area, fr.values));
	}
}

template<typename TModel, typename TMaterial>
void BorderConditions<TModel, CubicGrid, TMaterial>::applyBorderBeforeStage
(Mesh* mesh_, const real timeStep_, const int stage) {
	// handling borders
	mesh = mesh_; timeStep = timeStep_; direction = stage;
	// special for x-axis (because MPI partitioning along x-axis)
	if (direction == 0) {
		if (Engine::Instance().getForceSequence() ||
		    Engine::Instance().MpiRank == 0) {
			onTheRight = false;
			handleSide();
		}
		if (Engine::Instance().getForceSequence() ||
		    Engine::Instance().MpiRank == Engine::Instance().MpiSize - 1) {
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

template<typename TModel, typename TMaterial>
void BorderConditions<TModel, CubicGrid, TMaterial>::handleSide() const {
	auto borderIter = mesh->slice(direction, 0);
	if (onTheRight)
		borderIter = mesh->slice(direction, mesh->sizes(direction) - 1);
	while (borderIter != borderIter.end()) {
		for (const auto& condition : conditions) {
			if (condition.area->contains(mesh->coords(borderIter))) {
				handleBorderPoint(borderIter, condition.values);
			}
		}
		++borderIter;
	}
}

template<typename TModel, typename TMaterial>
void BorderConditions<TModel, CubicGrid, TMaterial>::handleBorderPoint
(const Iterator& borderIter, const Map& values) const {
	
	int innerSign = onTheRight ? -1 : 1;
	for (int a = 1; a <= mesh->borderSize; a++) {
		auto realIter = borderIter; realIter(direction) += innerSign * a;
		auto virtIter = borderIter; virtIter(direction) -= innerSign * a;
	
		mesh->_pde(virtIter) = mesh->pde(realIter);
		for (const auto& q : values) {
			const auto& quantity = q.first;
			const auto& timeDependency = q.second;
			real realValue = PdeVariables::QUANTITIES.at(quantity).Get(mesh->pde(realIter));
			real virtValue = - realValue + 2 * timeDependency(Engine::Instance().getCurrentTime());
			PdeVariables::QUANTITIES.at(quantity).Set(virtValue, mesh->_pde(virtIter));
		}
	}
}

template<typename TModel, typename TMaterial>
void BorderConditions<TModel, CubicGrid, TMaterial>::applyBorderAfterStage
(Mesh* mesh_, const real timeStep_, const int stage) {
	// handling inner surfaces
	mesh = mesh_; timeStep = timeStep_; direction = stage;
	for (const auto& innerSurface : innerSurfaces) {
		if (innerSurface.direction == stage) {
			allocateHelpMesh();
			auto sliceIter = mesh->slice(direction, innerSurface.index);
			while (sliceIter != sliceIter.end()) {
				if (innerSurface.area->contains(mesh->coords(sliceIter))) {
					handleInnerSurfacePoint(sliceIter, innerSurface.values, innerSurface.normal);
				}
				++sliceIter;
			}
			delete helpMesh;
		}
	}
}

template<typename TModel, typename TMaterial>
void BorderConditions<TModel, CubicGrid, TMaterial>::allocateHelpMesh() {
	Task::CubicGrid helpCubicGridTask;
	helpCubicGridTask.dimensionality = 1;
	helpCubicGridTask.borderSize = mesh->borderSize;
	helpCubicGridTask.h = mesh->h;
	helpCubicGridTask.sizes = {mesh->borderSize, 1, 1};
	Task helpTask; helpTask.cubicGrid = helpCubicGridTask;
	helpMesh = new Mesh(helpTask);
	helpMesh->allocate();
}

template<typename TModel, typename TMaterial>
void BorderConditions<TModel, CubicGrid, TMaterial>::handleInnerSurfacePoint
(const Iterator& iter, const Map& values, const int surfaceNormal) {
	// copy values to helpMesh
	for (int i = 0; i < 2 * mesh->borderSize; i++) {
		Iterator tmpIter = iter;
		tmpIter(direction) += i * surfaceNormal;
		helpMesh->_pde({i, 0, 0}) = mesh->pde(tmpIter);
		helpMesh->_matrices({i, 0, 0}) = mesh->_matrices(tmpIter);
	}

	// apply border conditions before stage on the helpMesh
	onTheRight = false;
	Mesh* tmpMesh = mesh; int tmpDirection = direction;
	mesh = helpMesh; direction = 0;
	handleBorderPoint({0, 0, 0}, values);
	mesh = tmpMesh; direction = tmpDirection;

	// calculate stage on the helpMesh
	GridCharacteristicMethod<Mesh>::stage(0, timeStep * surfaceNormal, helpMesh);

	// copy calculated values to real mesh
	for (int i = 0; i < mesh->borderSize; i++) {
		Iterator tmpIter = iter;
		tmpIter(direction) += i * surfaceNormal;
		mesh->_pdeNew(tmpIter) = helpMesh->pdeNew({i, 0, 0});
	}
}


template class BorderConditions<Elastic1DModel, CubicGrid, IsotropicMaterial>;
template class BorderConditions<Elastic2DModel, CubicGrid, IsotropicMaterial>;
template class BorderConditions<Elastic3DModel, CubicGrid, IsotropicMaterial>;
template class BorderConditions<SuperDuperModel, CubicGrid, IsotropicMaterial>;
template class BorderConditions<Elastic3DModel, CubicGrid, OrthotropicMaterial>;
template class BorderConditions<SuperDuperModel, CubicGrid, OrthotropicMaterial>;


