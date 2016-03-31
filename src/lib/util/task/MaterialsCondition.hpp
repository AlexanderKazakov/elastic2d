#ifndef LIBGCM_MATERIALSCONDITION_HPP
#define LIBGCM_MATERIALSCONDITION_HPP

#include "Task.hpp"

namespace gcm {

template<typename TModel, typename TMaterial>
class MaterialsCondition {
public:
	typedef typename TModel::GCM_MATRICES GCM_MATRICES;
	typedef typename TModel::PdeVector    PdeVector;

	/** Convert task terms of material conditions to own format */
	MaterialsCondition(const Statement& statement) :
		defaultCondition(Condition
		                         (std::dynamic_pointer_cast<TMaterial>(statement.
		                                                               materialConditions.
		                                                               defaultMaterial))) {
		for (const auto& m : statement.materialConditions.materials) {
			conditions.push_back(Condition
			                             (std::dynamic_pointer_cast<TMaterial>(m.
			                                                                   material),
			                             m.area));
		}
	}

	template<typename TNode>
	void apply(TNode node) const {
		node->_material() = defaultCondition.material;
		node->_matrices() = defaultCondition.matrices;
		for (const auto& condition : conditions) {
			if (condition.area->contains(node->coords())) {
				node->_material() = condition.material;
				node->_matrices() = condition.matrices;
			}
		}
	}

	real getMaximalEigenvalue() const {
		real ans = defaultCondition.matrices->getMaximalEigenvalue();
		for (const auto& c : conditions) {
			ans = fmax(ans, c.matrices->getMaximalEigenvalue());
		}
		return ans;
	}

private:
	struct Condition {
		Condition(std::shared_ptr<TMaterial> material_,
		          std::shared_ptr<Area> area_ = std::shared_ptr<Area>()) :
			area(area_), material(material_) {
			matrices = std::make_shared<GCM_MATRICES>();
			TModel::constructGcmMatrices(matrices, material);
		}

		std::shared_ptr<Area> area;
		std::shared_ptr<TMaterial> material;
		std::shared_ptr<GCM_MATRICES> matrices;
	};
	Condition defaultCondition;
	std::vector<Condition> conditions;

};


}


#endif // LIBGCM_MATERIALSCONDITION_HPP
