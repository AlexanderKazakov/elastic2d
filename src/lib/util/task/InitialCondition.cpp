#include <lib/util/task/InitialCondition.hpp>
#include <lib/grid/nodes/Node.hpp>


using namespace gcm;

template<class TNode>
void InitialCondition<TNode>::initialize(const Task &task) {

	for (auto& vectorInitCondition : task.initialCondition.vectors) {
		assert_eq(TNode::M, vectorInitCondition.list.size());
		conditions.push_back(Condition(vectorInitCondition.area, Vector(vectorInitCondition.list)));
	}

	for (auto& wave : task.initialCondition.waves) {
		assert_lt(wave.direction, TNode::DIMENSIONALITY);
		GcmMatrices gcmMatrices(task.material); 
		auto A = gcmMatrices.A(wave.direction);
		int columnNumber = GcmMatrices::WAVE_COLUMNS.at(wave.waveType);
		TNode tmp;
		tmp.u = A.U1.getColumn(columnNumber);
		real currentValue = Vector::QUANTITIES.at(wave.quantity).Get(tmp.u);
		assert_ne(currentValue, 0.0);
		tmp.u *= wave.quantityValue / currentValue;
		conditions.push_back(Condition(wave.area, tmp.u));
	}

	for (auto& quantityInitCondition : task.initialCondition.quantities) {
		TNode tmp;
		linal::clear(tmp.u);
		Vector::QUANTITIES.at(quantityInitCondition.physicalQuantity).Set(quantityInitCondition.value, tmp.u);
		conditions.push_back(Condition(quantityInitCondition.area, tmp.u));
	}
}

template<class TNode>
void InitialCondition<TNode>::apply(TNode &node, const linal::Vector3 &coords) const {
	linal::clear(node.u);
	for (auto& condition : conditions) {
		if (condition.area->contains(coords)) {
			node.u += condition.vector;
		}
	}
}


template class InitialCondition<IdealElastic1DNode>;
template class InitialCondition<IdealElastic2DNode>;
template class InitialCondition<IdealElastic3DNode>;