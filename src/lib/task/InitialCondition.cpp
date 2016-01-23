#include <lib/task/InitialCondition.hpp>
#include <lib/model/IdealElastic1DModel.hpp>
#include <lib/model/IdealElastic2DModel.hpp>
#include <lib/model/IdealElastic3DModel.hpp>


using namespace gcm;

template<class TModel>
void InitialCondition<TModel>::initialize(const Task &task, std::shared_ptr<GCM_MATRICES> gcmMatrices) {
	for (auto& vector : task.initialCondition.vectors) {
		assert_eq(Node::M, vector.list.size());
		conditions.push_back(Condition(vector.area, Vector(vector.list)));
	}
	for (auto& wave : task.initialCondition.waves) {
		assert_lt(wave.direction, TModel::DIMENSIONALITY);
		auto& A = gcmMatrices->A(wave.direction);
		int columnNumber = GCM_MATRICES::WAVE_COLUMNS.at(wave.waveType);
		Node node;
		node.matrix = gcmMatrices;
		node.u = A.U1.getColumn(columnNumber);
		real currentValue = TModel::QUANTITIES.at(wave.quantity).Get(node);
		assert_ne(currentValue, 0.0);
		node.u *= wave.quantityValue / currentValue;
		conditions.push_back(Condition(wave.area, node.u));
	}
	for (auto& quantity : task.initialCondition.quantities) {
		Node node;
		linal::clear(node.u);
		node.matrix = gcmMatrices;
		TModel::QUANTITIES.at(quantity.physicalQuantity).Set(quantity.value, node);
		conditions.push_back(Condition(quantity.area, node.u));
	}
}

template<class TModel>
void InitialCondition<TModel>::apply(Node &node, const linal::Vector3 &coords) const {
	linal::clear(node.u);
	for (auto& condition : conditions) {
		if (condition.area->contains(coords)) {
			node.u += condition.vector;
		}
	}
}


template class InitialCondition<IdealElastic1DModel>;
template class InitialCondition<IdealElastic2DModel>;
template class InitialCondition<IdealElastic3DModel>;