#include <lib/util/task/InitialCondition.hpp>
#include <lib/rheology/models/Model.hpp>


using namespace gcm;

template<class TModel>
void InitialCondition<TModel>::initialize(const Statement &stmt) {

	for (auto& vectorInitCondition : stmt.initialCondition.vectors) {
		assert_eq(PdeVector::M, vectorInitCondition.list.size());
		pdeConditions.push_back(PdeCondition(vectorInitCondition.area, 
			                 PdeVector(vectorInitCondition.list)));
	}

	for (auto& wave : stmt.initialCondition.waves) {
		assert_lt(wave.direction, TModel::DIMENSIONALITY);
		typename TModel::Material material;
		material.initialize(stmt);
		GCM_MATRICES gcmMatrices(material);
		auto A = gcmMatrices.A(wave.direction);
		int columnNumber = GCM_MATRICES::WAVE_COLUMNS.at(wave.waveType);
		PdeVector tmp;
		tmp = A.U1.getColumn(columnNumber);
		real currentValue = PdeVector::QUANTITIES.at(wave.quantity).Get(tmp);
		assert_ne(currentValue, 0.0);
		tmp *= wave.quantityValue / currentValue;
		pdeConditions.push_back(PdeCondition(wave.area, tmp));
	}

	for (auto& quantityInitCondition : stmt.initialCondition.quantities) {
		PdeVector tmp;
		linal::clear(tmp);
		PdeVector::QUANTITIES.at(quantityInitCondition.physicalQuantity).Set
				(quantityInitCondition.value, tmp);
		pdeConditions.push_back(PdeCondition(quantityInitCondition.area, tmp));
	}
}

template<class TModel>
void InitialCondition<TModel>::apply(PdeVector &v, const linal::Vector3& coords) const {
	linal::clear(v);
	for (auto& condition : pdeConditions) {
		if (condition.area->contains(coords)) {
			v += condition.pdeVector;
		}
	}
}


template class InitialCondition<Elastic1DModel>;
template class InitialCondition<Elastic2DModel>;
template class InitialCondition<Elastic3DModel>;
template class InitialCondition<OrthotropicElastic3DModel>;
template class InitialCondition<ContinualDamageElastic2DModel>;
template class InitialCondition<IdealPlastic2DModel>;
template class InitialCondition<SuperDuperModel>;
