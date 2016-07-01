#ifndef LIBGCM_MATERIALSCONDITION_HPP
#define LIBGCM_MATERIALSCONDITION_HPP

#include <limits>

#include <lib/util/task/Task.hpp>

namespace gcm {

/**
 * Class for applying material condition by areas
 */
template<typename TModel, typename TMaterial>
class MaterialsCondition {
public:
	typedef typename TModel::GCM_MATRICES GCM_MATRICES;

	
	/** Convert task terms of material conditions to own format */
	MaterialsCondition(const Statement& statement) {

		conditions.push_back(Condition(
				statement.materialConditions.defaultMaterial,
				std::make_shared<InfiniteArea>()));
		
		for (const auto& m : statement.materialConditions.materials) {
			conditions.push_back(Condition(m.material, m.area));
		}
		
	}

	
	template<typename TNode>
	void apply(TNode node) const {
		for (const auto& condition : conditions) {
			if (condition.area->contains(node->coords())) {
				node->_material() = condition.material;
				node->_matrices() = condition.matrices;
			}
		}
	}

	
	real getMaximalEigenvalue() const {
		real ans = 0;
		for (const auto& c : conditions) {
			ans = fmax(ans, c.matrices->getMaximalEigenvalue());
		}
		return ans;
	}
	

private:
	
	struct Condition {
		Condition(Statement::MaterialCondition::Material material_,
				std::shared_ptr<Area> area_) {
			area = area_;
			material = std::dynamic_pointer_cast<TMaterial>(material_);
			assert_true(material);
			matrices = std::make_shared<GCM_MATRICES>();
			TModel::constructGcmMatrices(matrices, material);
		}

		std::shared_ptr<Area> area;
		std::shared_ptr<TMaterial> material;
		std::shared_ptr<GCM_MATRICES> matrices;
	};
	
	std::vector<Condition> conditions;

};



/**
 * Class for applying material condition by cells.
 * By now, supported for Cgal3DGrid only
 */
template<template<typename, typename, typename> class TMesh,
         typename TModel,
         typename TGrid,
         typename TMaterial>
class MaterialConditionByCells {
public:
	typedef TMesh<TModel, TGrid, TMaterial> Mesh;
	static void apply(Mesh& /*mesh*/, 
			const Statement::MaterialCondition::MaterialMap& /*materialMap*/) {
		THROW_UNSUPPORTED("Unsupported for this type of grid");
	}
};

class Cgal3DGrid;
template<template<typename, typename, typename> class TMesh,
         typename TModel,
         typename TMaterial>
class MaterialConditionByCells<TMesh, TModel, Cgal3DGrid, TMaterial> {
public:
	typedef TMesh<TModel, Cgal3DGrid, TMaterial>   Mesh;
	
	typedef typename TModel::GCM_MATRICES          GCM_MATRICES;
	typedef std::shared_ptr<TMaterial>             MaterialPtr;
	typedef std::shared_ptr<GCM_MATRICES>          MatricesPtr;
	
	typedef Statement::MaterialCondition::Priority     Priority;
	typedef Statement::MaterialCondition::MaterialFlag MaterialFlag;
	
	struct MaterialProps {
		MaterialPtr material;
		MatricesPtr matrices;
		Priority priority;
	};
	typedef std::map<MaterialFlag, MaterialProps> MaterialMap;
	
	
	static void apply(Mesh& mesh,
			const Statement::MaterialCondition::MaterialMap& task) {
		
		MaterialMap materials = initMaterials(task);
		
		/// set materials to nodes
		for (auto& it : mesh) {
			
			std::map<MaterialFlag, int> occurence;
			for (const auto& m : materials) {
				occurence.insert({m.first, 0});
			}
			
			for (const auto cell : mesh.neighborCells(it)) {
				MaterialFlag flag = mesh.cellInfo(cell);
				occurence.at(flag) += 1;
			}
			
			int maxWeight = -1;
			MaterialFlag chosen;
			for (auto& m : occurence) {
				MaterialFlag flag = m.first;
				Priority priority = materials.at(flag).priority;
				int weight = priority * m.second;
				if (weight > maxWeight) {
					maxWeight = weight;
					chosen = flag;
				}
			}
			assert_gt(maxWeight, 0);
			
			mesh._material(it) = materials.at(chosen).material;
			mesh._matrices(it) = materials.at(chosen).matrices;
		}

		/// set maximal in modulus eigenvalue to mesh
		real maxEval = 0;
		for (const auto& m : materials) {
			maxEval = fmax(maxEval, m.second.matrices->getMaximalEigenvalue());
		}
		mesh.maximalEigenvalue = maxEval;
	}
	
	
private:
	
	static MaterialMap initMaterials(
			const Statement::MaterialCondition::MaterialMap& task) {
	/// convert information from task to local format 
	/// and instantiate necessary structures
		
		assert_false(task.empty());
		
		MaterialMap ans;
		
		for (const auto& m : task) {
			MaterialFlag flag = m.first;
			
			MaterialPtr material = std::dynamic_pointer_cast<TMaterial>(m.second.first);
			assert_true(material);
			
			Priority priority = m.second.second;
			assert_ge(priority, 0);
			
			MatricesPtr matrices = std::make_shared<GCM_MATRICES>();
			TModel::constructGcmMatrices(matrices, material);
			
			ans.insert({
					flag, {material, matrices, priority}
			});
		}
		
		return ans;
	}
	
};


}


#endif // LIBGCM_MATERIALSCONDITION_HPP
