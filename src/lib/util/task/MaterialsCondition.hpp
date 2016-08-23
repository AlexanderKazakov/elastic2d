#ifndef LIBGCM_MATERIALSCONDITION_HPP
#define LIBGCM_MATERIALSCONDITION_HPP

#include <lib/util/task/Task.hpp>

namespace gcm {


/**
 * Class for applying material conditions (setting rheology to mesh nodes)
 */
template<typename TModel, typename TGrid, typename TMaterial,
         template<typename, typename, typename> class TMesh>
class MaterialsCondition {
public:
	typedef TMesh<TModel, TGrid, TMaterial>    Mesh;	
	typedef typename Mesh::GCM_MATRICES        GCM_MATRICES;
	typedef typename Mesh::GridId              GridId;
	
	
	/**
	 * Set materials, gcm matrices and maximal eigenvalue to mesh
	 * according to given task
	 */
	static void apply(const Task& task, Mesh* mesh) {
		
		Conditions conditions = convertToLocalFormat(task, mesh->id);
		
		for (const auto& it : *mesh) {
			for (const auto& condition : conditions) {
				if (condition.area->contains(mesh->coords(it))) {
					mesh->_material(it) = condition.material;
					mesh->_matrices(it) = condition.matrices;
				}
			}
		}
		
		mesh->maximalEigenvalue = getMaximalEigenvalue(conditions);
	}
	
	
	/// Local format of task conditions
	struct Condition {
		Condition(Task::MaterialCondition::Material material_,
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
	
	typedef std::vector<Condition> Conditions;
	
	
	
	static Conditions convertToLocalFormat(
			const Task& task, const GridId id) {
		
		Conditions conditions;
		
		switch (task.materialConditions.type) {
			case Task::MaterialCondition::Type::BY_AREAS:
			{
				conditions.push_back(Condition(
						task.materialConditions.byAreas.defaultMaterial,
						std::make_shared<InfiniteArea>()));
				
				for (const auto& m : task.materialConditions.byAreas.materials) {
					conditions.push_back(Condition(m.material, m.area));
				}
				
				break;
			}
			case Task::MaterialCondition::Type::BY_BODIES:
			{
				const auto material = task.
						materialConditions.byBodies.bodyMaterialMap.at(id);
				
				conditions.push_back(Condition(material,
						std::make_shared<InfiniteArea>()));
				
				break;
			}
			default:
			{
				THROW_UNSUPPORTED("Unknown type of material condition");
			}
		}
		
		return conditions;
	}
	
	
	static real getMaximalEigenvalue(const Conditions& conditions) {
		real ans = 0;
		for (const auto& c : conditions) {
			ans = fmax(ans, c.matrices->getMaximalEigenvalue());
		}
		return ans;
	}
	
	
};


}


#endif // LIBGCM_MATERIALSCONDITION_HPP
