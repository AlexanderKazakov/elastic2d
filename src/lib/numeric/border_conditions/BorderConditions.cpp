#include <lib/numeric/border_conditions/BorderConditions.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::initialize(const Task& task) {
	for (const auto& bc : task.borderConditions) {
		const auto& values = bc.values;
		for (const auto& q : values) {
			assert_eq(PdeVariables::QUANTITIES.count(q.first), 1);
		}
		conditions.push_back(Condition(bc.area, values));
	}
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::applyBorderConditions
(Mesh* mesh_, const real time_) {
	mesh = mesh_; time = time_;
	// special for x-axis (because MPI partition along x-axis)
	direction = 0;
	if (mesh->getRank() == 0) {
		onTheRight = false;
		handleSide();
	}
	if (mesh->getRank() == mesh->getNumberOfWorkers() - 1) {
		onTheRight = true;
		handleSide();
	}
	// for other axes
	for (int d = 1; d < Mesh::DIMENSIONALITY; d++) {
		direction = d;
		onTheRight = false;
		handleSide();
		onTheRight = true;
		handleSide();
	}
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::handleSide() const {
	auto borderIter = mesh->slice(direction, 0);
	if (onTheRight)
		borderIter = mesh->slice(direction, mesh->getSizes()(direction) - 1);
	while (borderIter != borderIter.end()) {
		for (const auto& condition : conditions) {
			if (condition.area->contains(mesh->coords(borderIter))) {
				handlePoint(borderIter, condition.values);
			}
		}
		++borderIter;
	}
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::handlePoint
(const PartIterator& borderIter, const Map& values) const {
	
	int innerSign = onTheRight ? -1 : 1;
	for (int a = 1; a <= mesh->accuracyOrder; a++) {
		auto realIter = borderIter; realIter(direction) += innerSign * a;
		auto virtIter = borderIter; virtIter(direction) -= innerSign * a;
	
		mesh->_pde(virtIter) = mesh->pde(realIter);
		for (const auto& q : values) {
			const auto& quantity = q.first;
			const auto& timeDependency = q.second;
			real realValue = PdeVariables::QUANTITIES.at(quantity).Get(mesh->pde(realIter));
			real virtValue = - realValue + 2 * timeDependency(time);
			PdeVariables::QUANTITIES.at(quantity).Set(virtValue, mesh->_pde(virtIter));
		}
	}
}


template class BorderConditions<Elastic1DModel, CubicGrid>;
template class BorderConditions<Elastic2DModel, CubicGrid>;
template class BorderConditions<Elastic3DModel, CubicGrid>;
template class BorderConditions<OrthotropicElastic3DModel, CubicGrid>;
template class BorderConditions<ContinualDamageElastic2DModel, CubicGrid>;
template class BorderConditions<IdealPlastic2DModel, CubicGrid>;
template class BorderConditions<SuperDuperModel, CubicGrid>;
