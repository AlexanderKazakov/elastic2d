#include <lib/util/task/InitialCondition.hpp>
#include <lib/rheology/models/Model.hpp>


using namespace gcm;

template<typename TModel, typename TMaterial>
void InitialCondition<TModel, TMaterial>::
initialize(const Statement& statement) {

	for (auto& v : statement.initialCondition.vectors) {
		assert_eq(PdeVector::M, v.list.size());
		pdeConditions.push_back(PdeCondition(v.area, PdeVector(v.list)));
	}

	for (auto& wave : statement.initialCondition.waves) {
		assert_lt(wave.direction, TModel::DIMENSIONALITY);

		auto material = std::dynamic_pointer_cast<TMaterial>
		                        (statement.materialConditions.defaultMaterial);
		auto gcmMatricesPtr = std::make_shared<GCM_MATRICES>();
		TModel::constructGcmMatrices(gcmMatricesPtr, material);

		auto A = gcmMatricesPtr->m[wave.direction];
		int columnNumber = TModel::MATERIALS_WAVES_MAP.at(TMaterial::ID).at(wave.waveType);
		PdeVector tmp = A.U1.getColumn(columnNumber);
		real currentValue = PdeVariables::QUANTITIES.at(wave.quantity).Get(tmp);
		assert_ne(currentValue, 0.0);
		tmp *= wave.quantityValue / currentValue;
		pdeConditions.push_back(PdeCondition(wave.area, tmp));
	}

	for (auto& q : statement.initialCondition.quantities) {
		PdeVector tmp;
		linal::clear(tmp);
		PdeVariables::QUANTITIES.at(q.physicalQuantity).Set(q.value, tmp);
		pdeConditions.push_back(PdeCondition(q.area, tmp));
	}
}


template<typename TModel, typename TMaterial>
void InitialCondition<TModel, TMaterial>::
apply(PdeVector& v, const Real3& coords) const {
	linal::clear(v);
	for (const auto& condition : pdeConditions) {
		if (condition.area->contains(coords)) {
			v += condition.pdeVector;
		}
	}
}


template class InitialCondition<Elastic1DModel, IsotropicMaterial>;
template class InitialCondition<Elastic2DModel, IsotropicMaterial>;
template class InitialCondition<Elastic2DModel, OrthotropicMaterial>;
template class InitialCondition<Elastic3DModel, IsotropicMaterial>;
template class InitialCondition<Elastic3DModel, OrthotropicMaterial>;
template class InitialCondition<SuperDuperModel, OrthotropicMaterial>;
template class InitialCondition<SuperDuperModel, IsotropicMaterial>;


