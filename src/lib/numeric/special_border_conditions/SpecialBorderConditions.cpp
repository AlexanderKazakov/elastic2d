#include <lib/numeric/special_border_conditions/SpecialBorderConditions.hpp>
#include <lib/rheology/models/Model.hpp>
#include <lib/numeric/gcm/grid_characteristic_methods.hpp>

using namespace gcm;


template<typename TModel, typename TMaterial, int Dimensionality>
void SpecialBorderConditions<TModel, CubicGrid<Dimensionality>, TMaterial>::
beforeStatement(const Statement& statement) {
	for (const auto& bc : statement.cubicGridBorderConditions) {
		conditions.push_back(Condition(bc.area, bc.values));
	}
	for (const auto& fr : statement.fractures) {
		innerSurfaces.push_back(
				InnerSurface(fr.direction, fr.index,    -1, fr.area, fr.values));
		innerSurfaces.push_back(
				InnerSurface(fr.direction, fr.index + 1, 1, fr.area, fr.values));
	}
}


template<typename TModel, typename TMaterial, int Dimensionality>
void SpecialBorderConditions<TModel, CubicGrid<Dimensionality>, TMaterial>::
applyBorderBeforeStage(Mesh* mesh, const real /*timeStep*/, const int stage) const {
	/// handling borders
	if (stage == 0) {
	// special for x-axis (because MPI partitioning along x-axis)
		if (Mpi::ForceSequence() || Mpi::Rank() == 0) {
			handleSide(mesh, stage, false);
		}
		if (Mpi::ForceSequence() || Mpi::Rank() == Mpi::Size() - 1) {
			handleSide(mesh, stage, true);
		}
		return;
	}
	// for other axes
	handleSide(mesh, stage, false);
	handleSide(mesh, stage, true);
}


template<typename TModel, typename TMaterial, int Dimensionality>
void SpecialBorderConditions<TModel, CubicGrid<Dimensionality>, TMaterial>::
handleSide(Mesh* mesh, const int direction, const bool onTheRight) const {
	auto borderIter = onTheRight ? mesh->slice(direction, mesh->sizes(direction) - 1)
	                             : mesh->slice(direction, 0);
	
	while (borderIter != borderIter.end()) {
		for (const auto& condition : conditions) {
			if (condition.area->contains(mesh->coords(borderIter))) {
				handleBorderPoint(mesh, borderIter, condition.values, direction, onTheRight);
			}
		}
		++borderIter;
	}
}


template<typename TModel, typename TMaterial, int Dimensionality>
void SpecialBorderConditions<TModel, CubicGrid<Dimensionality>, TMaterial>::
handleBorderPoint(Mesh* mesh, const Iterator& borderIter, const Map& values,
		const int direction, const bool onTheRight) {

	int innerSign = onTheRight ? -1 : 1;
	for (int a = 1; a <= mesh->borderSize; a++) {
		auto realIter = borderIter; realIter(direction) += innerSign * a;
		auto virtIter = borderIter; virtIter(direction) -= innerSign * a;

		mesh->_pde(virtIter) = mesh->pde(realIter);
		for (const auto& q : values) {
			const auto& quantity = q.first;
			const auto& timeDependency = q.second;
			real realValue =
			        PdeVariables::QUANTITIES.at(quantity).Get(mesh->pde(realIter));
			real virtValue = -realValue + 2 * timeDependency(Clock::Time());
			PdeVariables::QUANTITIES.at(quantity).Set(virtValue, mesh->_pde(virtIter));
		}
	}
}


template<typename TModel, typename TMaterial, int Dimensionality>
void SpecialBorderConditions<TModel, CubicGrid<Dimensionality>, TMaterial>::
applyBorderAfterStage(Mesh* mesh, const real timeStep, const int stage) const {
	/// handling inner surfaces
	for (const auto& innerSurface : innerSurfaces) {
		if (innerSurface.direction == stage) {
			HelpMesh* helpMesh = allocateHelpMesh(mesh);
			auto sliceIter = mesh->slice(stage, innerSurface.index);
			while (sliceIter != sliceIter.end()) {
				if (innerSurface.condition.area->contains(mesh->coords(sliceIter))) {
					handleInnerSurfacePoint(mesh, helpMesh, timeStep, stage,
							sliceIter, innerSurface.condition.values, innerSurface.normal);
				}
				++sliceIter;
			}
			delete helpMesh;
		}
	}
}


template<typename TModel, typename TMaterial, int Dimensionality>
void SpecialBorderConditions<TModel, CubicGrid<Dimensionality>, TMaterial>::
handleInnerSurfacePoint(Mesh* mesh, HelpMesh* helpMesh, const real timeStep,
		const int direction, const Iterator& iter, const Map& values, 
		const int surfaceNormal) {
	// copy values to helpMesh
	for (int i = 0; i < 2 * mesh->borderSize; i++) {
		Iterator tmpIter = iter;
		tmpIter(direction) += i * surfaceNormal;
		helpMesh->_pde({i}) = mesh->pde(tmpIter);
		helpMesh->_matrices({i}) = mesh->_matrices(tmpIter);
	}

	// apply border conditions before stage on the helpMesh
	HelpSpecBorderCond::handleBorderPoint(helpMesh, {0}, values, 0, false);

	// calculate stage on the helpMesh
	GridCharacteristicMethod<TModel, CubicGrid<1>, TMaterial>().stage(
			0, timeStep * surfaceNormal, *helpMesh);

	// copy calculated values to real mesh
	for (int i = 0; i < mesh->borderSize; i++) {
		Iterator tmpIter = iter;
		tmpIter(direction) += i * surfaceNormal;
		mesh->_pdeNew(tmpIter) = helpMesh->pdeNew({i});
	}
}



template class SpecialBorderConditions<Elastic1DModel, CubicGrid<1>, IsotropicMaterial>;
template class SpecialBorderConditions<Elastic2DModel, CubicGrid<2>, IsotropicMaterial>;
template class SpecialBorderConditions<Elastic3DModel, CubicGrid<3>, IsotropicMaterial>;
template class SpecialBorderConditions<SuperDuperModel, CubicGrid<3>, IsotropicMaterial>;

template class SpecialBorderConditions<Elastic2DModel, CubicGrid<2>, OrthotropicMaterial>;
template class SpecialBorderConditions<Elastic3DModel, CubicGrid<3>, OrthotropicMaterial>;
template class SpecialBorderConditions<SuperDuperModel, CubicGrid<3>, OrthotropicMaterial>;



