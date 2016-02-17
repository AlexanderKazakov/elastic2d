#include <lib/numeric/border_conditions/StructuredGridBorderConditions.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/grid/DefaultGrid.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TGrid>
void StructuredGridBorderConditions<TGrid>::initialize(const Task &task) {
	borderConditions = task.borderConditions;
}

template<class TGrid>
void StructuredGridBorderConditions<TGrid>::applyBorderConditions(TGrid* mesh) const {
	SUPPRESS_WUNUSED(mesh);
	// TODO
/*	if (mesh->getRank() == 0 &&
	    borderConditions.at(CUBIC_BORDERS::X_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int y = 0; y < mesh->sizes(1); y++) {
			for (int z = 0; z < mesh->sizes(2); z++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < TGrid::DIMENSIONALITY; j++) {
						mesh->_pde(-i, y, z).V[j] = mesh->pde(i, y, z).V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(-i, y, z).S[j] = -mesh->pde(i, y, z).S[j];
					}
				}
			}
		}
	}
	if (mesh->getRank() == mesh->getNumberOfWorkers() - 1 &&
	    borderConditions.at(CUBIC_BORDERS::X_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int y = 0; y < mesh->sizes(1); y++) {
			for (int z = 0; z < mesh->sizes(2); z++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < TGrid::DIMENSIONALITY; j++) {
						mesh->_pde(mesh->sizes(0) - 1 + i, y, z).V[j] = mesh->pde(mesh->sizes(0) - 1 - i, y, z).V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(mesh->sizes(0) - 1 + i, y, z).S[j] = -mesh->pde(mesh->sizes(0) - 1 - i, y, z).S[j];
					}
				}
			}
		}
	}

	if (mesh->sizes(1) != 1 && borderConditions.at(CUBIC_BORDERS::Y_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < mesh->sizes(0); x++) {
			for (int z = 0; z < mesh->sizes(2); z++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < TGrid::DIMENSIONALITY; j++) {
						mesh->_pde(x, -i, z).V[j] = mesh->pde(x, i, z).V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(x, -i, z).S[j] = -mesh->pde(x, i, z).S[j];
					}
				}
			}
		}
	}
	if (mesh->sizes(1) != 1 && borderConditions.at(CUBIC_BORDERS::Y_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < mesh->sizes(0); x++) {
			for (int z = 0; z < mesh->sizes(2); z++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < TGrid::DIMENSIONALITY; j++) {
						mesh->_pde(x, mesh->sizes(1) - 1 + i, z).V[j] = mesh->pde(x, mesh->sizes(1) - 1 - i, z).V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(x, mesh->sizes(1) - 1 + i, z).S[j] = -mesh->pde(x, mesh->sizes(1) - 1 - i, z).S[j];
					}
				}
			}
		}
	}

	if (mesh->sizes(2) != 1 && borderConditions.at(CUBIC_BORDERS::Z_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < mesh->sizes(0); x++) {
			for (int y = 0; y < mesh->sizes(1); y++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < TGrid::DIMENSIONALITY; j++) {
						mesh->_pde(x, y, -i).V[j] = mesh->pde(x, y, i).V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(x, y, -i).S[j] = -mesh->pde(x, y, i).S[j];
					}
				}
			}
		}
	}
	if (mesh->sizes(2) != 1 && borderConditions.at(CUBIC_BORDERS::Z_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < mesh->sizes(0); x++) {
			for (int y = 0; y < mesh->sizes(1); y++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < TGrid::DIMENSIONALITY; j++) {
						mesh->_pde(x, y, mesh->sizes(2) - 1 + i).V[j] = mesh->pde(x, y, mesh->sizes(2) - 1 - i).V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(x, y, mesh->sizes(2) - 1 + i).S[j] = -mesh->pde(x, y, mesh->sizes(2) - 1 - i).S[j];
					}
				}
			}
		}
	}*/

}

/*template<class TGrid>
void StructuredGridBorderConditions<TGrid>::applyBorderConditions(TGrid* mesh) const {
	applyBorderCondition(0, true, mesh, CUBIC_BORDERS::X_LEFT);
	applyBorderCondition(0, false, mesh, CUBIC_BORDERS::X_RIGHT);
	applyBorderCondition(1, true, mesh, CUBIC_BORDERS::Y_LEFT);
	applyBorderCondition(1, false, mesh, CUBIC_BORDERS::Y_RIGHT);
	applyBorderCondition(2, true, mesh, CUBIC_BORDERS::Z_LEFT);
	applyBorderCondition(2, false, mesh, CUBIC_BORDERS::Z_RIGHT);
}*/
/*template<class TGrid>
void StructuredGridBorderConditions<TGrid>::applyBorderCondition
		(const int direction, const bool leftCorner, TGrid* mesh, CUBIC_BORDERS border) const {
	if (direction == 0) {
		if ( ( leftCorner && mesh->getRank() != 0) ||
		     (!leftCorner && mesh->getRank() != mesh->getNumberOfWorkers() - 1) )
			return;
	}
	if (borderConditions.at(border) != BorderCondition::T::FREE_BORDER) return;
	if (mesh->sizes(direction) <= 1) return;

	auto sizes = mesh->sizes; sizes(direction) = mesh->accuracyOrder;
	typename TGrid::ForwardIterator begin({0,0,0}, sizes);
	typename TGrid::ForwardIterator end({sizes(0),0,0}, sizes);

	for (auto it = begin; it != end; ++it) {
		typename TGrid::ForwardIterator from = it;
		typename TGrid::ForwardIterator to = it;
		if (leftCorner) {
			from(direction) =   from(direction) + 1;
			to(direction)   = - to(direction) + 1;
		} else {
			from(direction) = mesh->sizes(direction) - from(direction);
			from(direction) = mesh->sizes(direction) + to(direction);
		}
		for (int j = 0; j < TGrid::DIMENSIONALITY; j++) {
			mesh->_pde(to).V[j] = mesh->pde(from).V[j];
		}
		for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
			mesh->_pde(to).S[j] = -mesh->pde(from).S[j];
		}
	}
}*/

template class StructuredGridBorderConditions<DefaultGrid<Elastic1DModel, StructuredGrid>>;
template class StructuredGridBorderConditions<DefaultGrid<Elastic2DModel, StructuredGrid>>;
template class StructuredGridBorderConditions<DefaultGrid<Elastic3DModel, StructuredGrid>>;
template class StructuredGridBorderConditions<DefaultGrid<OrthotropicElastic3DModel, StructuredGrid>>;
template class StructuredGridBorderConditions<DefaultGrid<ContinualDamageElastic2DModel, StructuredGrid>>;
template class StructuredGridBorderConditions<DefaultGrid<IdealPlastic2DModel, StructuredGrid>>;
template class StructuredGridBorderConditions<DefaultGrid<SuperDuperModel, StructuredGrid>>;
