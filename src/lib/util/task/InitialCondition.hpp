#ifndef LIBGCM_INITIALCONDITION_HPP
#define LIBGCM_INITIALCONDITION_HPP

#include <lib/util/task/Task.hpp>
#include <lib/util/task/MaterialsCondition.hpp>

namespace gcm {

template<typename TModel, typename TGrid, typename TMaterial,
         template<typename, typename, typename> class TMesh>
class InitialCondition {
public:
	typedef TMesh<TModel, TGrid, TMaterial>  Mesh;
	typedef typename Mesh::PdeVariables      PdeVariables;
	typedef typename Mesh::PdeVector         PdeVector;
	typedef typename Mesh::GCM_MATRICES      GCM_MATRICES;
	typedef typename Mesh::GridId            GridId;
	
	typedef MaterialsCondition<TModel, TGrid, TMaterial, TMesh> MC;
	
	
	/**
	 * Apply initial conditions to mesh according to given task
	 */
	static void apply(const Task& task, Mesh* mesh) {
		
		Conditions conditions = convertToLocalFormat(task, mesh->id);
		
		for (const auto& it : *mesh) {
			linal::clear(mesh->_pde(it));
			for (const auto& condition : conditions) {
				if (condition.area->contains(mesh->coords(it))) {
					mesh->_pde(it) += condition.pdeVector;
				}
			}
		}
		
	}
	
	
private:
	
	struct Condition {
		std::shared_ptr<Area> area;
		PdeVector pdeVector;
	};
	
	typedef std::vector<Condition> Conditions;
	
	
	static Conditions convertToLocalFormat(
			const Task& task, const GridId gridId) {
		
		Conditions conditions;
		
		for (auto& v : task.initialCondition.vectors) {
			assert_eq(PdeVector::M, v.list.size());
			conditions.push_back({v.area, PdeVector(v.list)});
		}
		
		
		typename MC::Conditions mcConditions = 
				MC::convertToLocalFormat(task, gridId);
		auto material = mcConditions.front().material;
		auto gcmMatricesPtr = std::make_shared<GCM_MATRICES>();
		TModel::constructGcmMatrices(gcmMatricesPtr, material);
		
		for (auto& wave : task.initialCondition.waves) {
			assert_lt(wave.direction, TModel::DIMENSIONALITY);
			
			auto A = (*gcmMatricesPtr)(wave.direction);
			int columnNumber = TModel::MATERIALS_WAVES_MAP.at(TMaterial::ID).at(wave.waveType);
			
			PdeVector tmp = A.U1.getColumn(columnNumber);
			real currentValue = PdeVariables::QUANTITIES.at(wave.quantity).Get(tmp);
			assert_ne(currentValue, 0);
			
			tmp *= wave.quantityValue / currentValue;
			conditions.push_back({wave.area, tmp});
		}
		
		
		for (auto& q : task.initialCondition.quantities) {
			PdeVector tmp;
			linal::clear(tmp);
			PdeVariables::QUANTITIES.at(q.physicalQuantity).Set(q.value, tmp);
			conditions.push_back({q.area, tmp});
		}
		
		return conditions;
	}
	
};


}


#endif // LIBGCM_INITIALCONDITION_HPP
