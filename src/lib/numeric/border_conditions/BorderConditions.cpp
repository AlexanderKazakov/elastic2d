#include <lib/numeric/border_conditions/BorderConditions.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::initialize(const Task &task) {
	borderConditions = task.borderConditions;
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::applyBorderConditions(Mesh* mesh) const {
	for (const auto& direction : borderConditions) {
		if ((int)direction.first >= Mesh::DIMENSIONALITY) break;
		switch (direction.second.first) { // left
			case BorderCondition::T::FREE_BORDER:
				handleSide(mesh, direction.first, false); break;
			default: break;
		}
		switch (direction.second.second) { // right
			case BorderCondition::T::FREE_BORDER:
				handleSide(mesh, direction.first, true); break;
			default: break;
		}
	}
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::handleSide(Mesh* mesh,
		DIRECTION direction, const bool onTheRight) const {
	if (direction == DIRECTION::X && !onTheRight 
		&& mesh->getRank() != 0) return;
	if (direction == DIRECTION::X &&  onTheRight 
		&& mesh->getRank() != mesh->getNumberOfWorkers() - 1) return;
	
	for (int a = 1; a <= mesh->getAccuracyOrder(); a++) {
		int realSlice = a, virtSlice = -a;
		if (onTheRight) {
			realSlice = mesh->getSizes()((int)direction) - 1 - a;
			virtSlice = mesh->getSizes()((int)direction) - 1 + a;
		}
		auto realIter = mesh->slice(direction, realSlice);
		auto virtIter = mesh->slice(direction, virtSlice);
		handleSlice(mesh, realIter, virtIter);
	}
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::handleSlice(Mesh* mesh,
		PartIterator realIter, PartIterator virtIter) const {
	while (realIter != realIter.end()) {
		mesh->_pde(virtIter) = mesh->pde(realIter);
		for (int j = 0; j < (Mesh::DIMENSIONALITY * (Mesh::DIMENSIONALITY + 1)) / 2; j++) {
			mesh->_pde(virtIter).S[j] *= -1;
		}
		++realIter; ++virtIter;
	}
}


template class BorderConditions<Elastic1DModel, CubicGrid>;
template class BorderConditions<Elastic2DModel, CubicGrid>;
template class BorderConditions<Elastic3DModel, CubicGrid>;
template class BorderConditions<OrthotropicElastic3DModel, CubicGrid>;
template class BorderConditions<ContinualDamageElastic2DModel, CubicGrid>;
template class BorderConditions<IdealPlastic2DModel, CubicGrid>;
template class BorderConditions<SuperDuperModel, CubicGrid>;
