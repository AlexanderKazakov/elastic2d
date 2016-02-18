#include <lib/numeric/border_conditions/BorderConditions.hpp>
#include <lib/rheology/models/Model.hpp>

using namespace gcm;

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::initialize(const Task &task) {
	borderConditions = task.borderConditions;
}

template<typename TModel>
void BorderConditions<TModel, CubicGrid>::applyBorderConditions(Mesh* mesh) const {
	// TODO - make it beautiful with some slice iterator
	if (mesh->getRank() == 0 &&
	    borderConditions.at(CUBIC_BORDERS::X_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int y = 0; y < mesh->sizes(1); y++) {
			for (int z = 0; z < mesh->sizes(2); z++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < Mesh::DIMENSIONALITY; j++) {
						mesh->_pde(linal::VectorInt<3>({-i, y, z})).V[j] =
								mesh->pde(linal::VectorInt<3>({i, y, z})).V[j];
					}
					for (int j = 0; j < (Mesh::DIMENSIONALITY * (Mesh::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(linal::VectorInt<3>({-i, y, z})).S[j] =
								-mesh->pde(linal::VectorInt<3>({i, y, z})).S[j];
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
					for (int j = 0; j < Mesh::DIMENSIONALITY; j++) {
						mesh->_pde(linal::VectorInt<3>({mesh->sizes(0) - 1 + i, y, z})).V[j] =
								mesh->pde(linal::VectorInt<3>({mesh->sizes(0) - 1 - i, y, z})).V[j];
					}
					for (int j = 0; j < (Mesh::DIMENSIONALITY * (Mesh::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(linal::VectorInt<3>({mesh->sizes(0) - 1 + i, y, z})).S[j] =
								-mesh->pde(linal::VectorInt<3>({mesh->sizes(0) - 1 - i, y, z})).S[j];
					}
				}
			}
		}
	}

	if (mesh->sizes(1) != 1 && borderConditions.at(CUBIC_BORDERS::Y_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < mesh->sizes(0); x++) {
			for (int z = 0; z < mesh->sizes(2); z++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < Mesh::DIMENSIONALITY; j++) {
						mesh->_pde(linal::VectorInt<3>({x, -i, z})).V[j] =
								mesh->pde(linal::VectorInt<3>({x, i, z})).V[j];
					}
					for (int j = 0; j < (Mesh::DIMENSIONALITY * (Mesh::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(linal::VectorInt<3>({x, -i, z})).S[j] =
								-mesh->pde(linal::VectorInt<3>({x, i, z})).S[j];
					}
				}
			}
		}
	}
	if (mesh->sizes(1) != 1 && borderConditions.at(CUBIC_BORDERS::Y_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < mesh->sizes(0); x++) {
			for (int z = 0; z < mesh->sizes(2); z++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < Mesh::DIMENSIONALITY; j++) {
						mesh->_pde(linal::VectorInt<3>({x, mesh->sizes(1) - 1 + i, z})).V[j] =
								mesh->pde(linal::VectorInt<3>({x, mesh->sizes(1) - 1 - i, z})).V[j];
					}
					for (int j = 0; j < (Mesh::DIMENSIONALITY * (Mesh::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(linal::VectorInt<3>({x, mesh->sizes(1) - 1 + i, z})).S[j] =
								-mesh->pde(linal::VectorInt<3>({x, mesh->sizes(1) - 1 - i, z})).S[j];
					}
				}
			}
		}
	}

	if (mesh->sizes(2) != 1 && borderConditions.at(CUBIC_BORDERS::Z_LEFT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < mesh->sizes(0); x++) {
			for (int y = 0; y < mesh->sizes(1); y++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < Mesh::DIMENSIONALITY; j++) {
						mesh->_pde(linal::VectorInt<3>({x, y, -i})).V[j] =
								mesh->pde(linal::VectorInt<3>({x, y, i})).V[j];
					}
					for (int j = 0; j < (Mesh::DIMENSIONALITY * (Mesh::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(linal::VectorInt<3>({x, y, -i})).S[j] =
								-mesh->pde(linal::VectorInt<3>({x, y, i})).S[j];
					}
				}
			}
		}
	}
	if (mesh->sizes(2) != 1 && borderConditions.at(CUBIC_BORDERS::Z_RIGHT) == BorderCondition::T::FREE_BORDER) {
		for (int x = 0; x < mesh->sizes(0); x++) {
			for (int y = 0; y < mesh->sizes(1); y++) {
				for (int i = 1; i <= mesh->accuracyOrder; i++) {
					for (int j = 0; j < Mesh::DIMENSIONALITY; j++) {
						mesh->_pde(linal::VectorInt<3>({x, y, mesh->sizes(2) - 1 + i})).V[j] =
								mesh->pde(linal::VectorInt<3>({x, y, mesh->sizes(2) - 1 - i})).V[j];
					}
					for (int j = 0; j < (Mesh::DIMENSIONALITY * (Mesh::DIMENSIONALITY + 1)) / 2; j++) {
						mesh->_pde(linal::VectorInt<3>({x, y, mesh->sizes(2) - 1 + i})).S[j] =
								-mesh->pde(linal::VectorInt<3>({x, y, mesh->sizes(2) - 1 - i})).S[j];
					}
				}
			}
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
