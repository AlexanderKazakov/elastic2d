#include <lib/numeric/special_border_conditions/SpecialBorderConditions.hpp>
#include <lib/rheology/models/models.hpp>
#include <lib/numeric/gcm/GridCharacteristicMethodCubicGrid.hpp>

using namespace gcm;


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


template class SpecialBorderConditions<ElasticModel<1>, CubicGrid<1>, IsotropicMaterial>;
template class SpecialBorderConditions<ElasticModel<2>, CubicGrid<2>, IsotropicMaterial>;
template class SpecialBorderConditions<ElasticModel<3>, CubicGrid<3>, IsotropicMaterial>;

template class SpecialBorderConditions<ElasticModel<2>, CubicGrid<2>, OrthotropicMaterial>;
template class SpecialBorderConditions<ElasticModel<3>, CubicGrid<3>, OrthotropicMaterial>;

template class SpecialBorderConditions<AcousticModel<1>, CubicGrid<1>, IsotropicMaterial>;
template class SpecialBorderConditions<AcousticModel<2>, CubicGrid<2>, IsotropicMaterial>;
template class SpecialBorderConditions<AcousticModel<3>, CubicGrid<3>, IsotropicMaterial>;

