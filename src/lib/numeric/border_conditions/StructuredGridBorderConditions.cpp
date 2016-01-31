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
						(*mesh)(-i, y, z).u.V[j] = mesh->get(i, y, z).u.V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						(*mesh)(-i, y, z).u.S[j] = -mesh->get(i, y, z).u.S[j];
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
						(*mesh)(mesh->sizes(0) - 1 + i, y, z).u.V[j] = mesh->get(mesh->sizes(0) - 1 - i, y, z).u.V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						(*mesh)(mesh->sizes(0) - 1 + i, y, z).u.S[j] = -mesh->get(mesh->sizes(0) - 1 - i, y, z).u.S[j];
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
						(*mesh)(x, -i, z).u.V[j] = mesh->get(x, i, z).u.V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						(*mesh)(x, -i, z).u.S[j] = -mesh->get(x, i, z).u.S[j];
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
						(*mesh)(x, mesh->sizes(1) - 1 + i, z).u.V[j] = mesh->get(x, mesh->sizes(1) - 1 - i, z).u.V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						(*mesh)(x, mesh->sizes(1) - 1 + i, z).u.S[j] = -mesh->get(x, mesh->sizes(1) - 1 - i, z).u.S[j];
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
						(*mesh)(x, y, -i).u.V[j] = mesh->get(x, y, i).u.V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						(*mesh)(x, y, -i).u.S[j] = -mesh->get(x, y, i).u.S[j];
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
						(*mesh)(x, y, mesh->sizes(2) - 1 + i).u.V[j] = mesh->get(x, y, mesh->sizes(2) - 1 - i).u.V[j];
					}
					for (int j = 0; j < (TGrid::DIMENSIONALITY * (TGrid::DIMENSIONALITY + 1)) / 2; j++) {
						(*mesh)(x, y, mesh->sizes(2) - 1 + i).u.S[j] = -mesh->get(x, y, mesh->sizes(2) - 1 - i).u.S[j];
					}
				}
			}
		}
	}

}



template class StructuredGridBorderConditions<DefaultStructuredGrid<Elastic1DModel>>;
template class StructuredGridBorderConditions<DefaultStructuredGrid<Elastic2DModel>>;
template class StructuredGridBorderConditions<DefaultStructuredGrid<Elastic3DModel>>;
template class StructuredGridBorderConditions<DefaultStructuredGrid<OrthotropicElastic3DModel>>;
