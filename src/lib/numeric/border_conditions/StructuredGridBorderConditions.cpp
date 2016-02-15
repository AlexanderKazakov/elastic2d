#include <lib/numeric/border_conditions/StructuredGridBorderConditions.hpp>
#include <lib/grid/StructuredGrid.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<class TGrid>
void StructuredGridBorderConditions<TGrid>::initialize(const Task &task) {
	borderConditions = task.borderConditions;
}

template<class TGrid>
void StructuredGridBorderConditions<TGrid>::applyBorderConditions(TGrid* mesh) {
	if (mesh->getRank() == 0 &&
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
	}

}



template class StructuredGridBorderConditions<StructuredGrid<Elastic1DModel>>;
template class StructuredGridBorderConditions<StructuredGrid<Elastic2DModel>>;
template class StructuredGridBorderConditions<StructuredGrid<Elastic3DModel>>;
template class StructuredGridBorderConditions<StructuredGrid<OrthotropicElastic3DModel>>;
template class StructuredGridBorderConditions<StructuredGrid<ContinualDamageElastic2DModel>>;
template class StructuredGridBorderConditions<StructuredGrid<IdealPlastic2DModel>>;
template class StructuredGridBorderConditions<StructuredGrid<SuperDuperModel>>;
