#ifndef LIBGCM_SPECIALBORDERCONDITIONS_HPP
#define LIBGCM_SPECIALBORDERCONDITIONS_HPP

#include <lib/mesh/DefaultMesh.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/mesh/grid/CubicGlobalScene.hpp>


namespace gcm {

/**
 * Border calculations, specific for some grids, models, etc
 */
template<typename TModel, typename TGrid, typename TMaterial> 
struct SpecialBorderConditions {
	typedef DefaultMesh<TModel, TGrid, TMaterial> Mesh;
	
	SpecialBorderConditions(const Task&) { }
	void applyBorderBeforeStage(Mesh*, const real, const int) { }
	void applyBorderAfterStage(Mesh*, const real, const int) { }
};


template<typename TModel, typename TMaterial, int Dimensionality>
struct SpecialBorderConditions<TModel, CubicGrid<Dimensionality>, TMaterial> {
	typedef CubicGrid<Dimensionality>                       Grid;
	typedef DefaultMesh<TModel, Grid, TMaterial>            Mesh;
	typedef typename Mesh::PdeVector                        PdeVector;
	typedef typename Mesh::Iterator                         Iterator;
	typedef typename TModel::PdeVariables                   PdeVariables;
	
	typedef std::function<real(real)>                       TimeDependency;
	typedef std::map<PhysicalQuantities::T, TimeDependency> Map;
	
	struct Condition {
		Condition(const std::shared_ptr<Area> area_, const Map& values_) :
				area(area_), values(values_) {
			for (const auto& q : values) {
				assert_eq(PdeVariables::QUANTITIES.count(q.first), 1);
			}
		}
		std::shared_ptr<Area> area; ///< area of surface
		Map values; ///< fixed values
	};
	
	struct InnerSurface {
		InnerSurface(const int direction_, const int index_, const int normal_,
		             std::shared_ptr<Area> area, const Map& values) :
				direction(direction_), index(index_), normal(normal_),
				condition(area, values) { }
		int direction;       ///< crossing axis
		int index;           ///< index at crossing axis
		int normal;          ///< -1 or 1 only
		Condition condition; ///< border condition od the surface
	};
	
	SpecialBorderConditions(const Task& task) {
		for (const auto& bc : task.cubicGridBorderConditions) {
			conditions.push_back(Condition(bc.area, bc.values));
		}
	}

	void applyBorderBeforeStage(
			Mesh* mesh, const real timeStep, const int stage) const;
	void applyBorderAfterStage(Mesh*, const real, const int) const { }
	
	void handleSide(Mesh* mesh, const int direction, const bool onTheRight) const;
	
	static void handleBorderPoint(Mesh* mesh, const Iterator& borderIter,
			const Map& values, const int direction, const bool onTheRight);
	
	
private:
	/// list of border conditions applied in sequence (overwriting previous)
	std::vector<Condition> conditions;
	
};


}

#endif // LIBGCM_SPECIALBORDERCONDITIONS_HPP
