#ifndef LIBGCM_SPECIALBORDERCONDITIONS_HPP
#define LIBGCM_SPECIALBORDERCONDITIONS_HPP

#include <lib/mesh/DefaultMesh.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>


namespace gcm {
/**
 * Border calculations, specific for some grids, models, etc
 */
template<typename TModel, typename TGrid, typename TMaterial> 
struct SpecialBorderConditions {
	typedef DefaultMesh<TModel, TGrid, TMaterial> Mesh;

	void beforeStatement(const Statement&) { }
	void applyBorderBeforeStage(Mesh*, const real, const int) { }
	void applyBorderAfterStage(Mesh*, const real, const int) { }
};


template<typename TModel, typename TMaterial, int Dimensionality>
struct SpecialBorderConditions<TModel, CubicGrid<Dimensionality>, TMaterial> {
	typedef CubicGrid<Dimensionality>                       Grid;
	typedef DefaultMesh<TModel, Grid, TMaterial>            Mesh;
	typedef CubicGrid<1>                                    HelpGrid;
	typedef DefaultMesh<TModel, HelpGrid, TMaterial>        HelpMesh;
	typedef SpecialBorderConditions<TModel, HelpGrid, TMaterial>
	                                                        HelpSpecBorderCond;

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

	void beforeStatement(const Statement& statement);
	void applyBorderBeforeStage(Mesh* mesh, const real timeStep, const int stage) const;
	void applyBorderAfterStage(Mesh* mesh, const real timeStep, const int stage) const;

	void handleSide(Mesh* mesh, const int direction, const bool onTheRight) const;

	static void handleBorderPoint(Mesh* mesh, const Iterator& borderIter,
			const Map& values, const int direction, const bool onTheRight);

	static void handleInnerSurfacePoint(Mesh* mesh, HelpMesh* helpMesh,
			const real timeStep, const int direction, const Iterator& iter, 
			const Map& values, const int surfaceNormal);
	
	static HelpMesh* allocateHelpMesh(const Mesh* mesh, const int stage) {
		Task helpTask;
		helpTask.cubicGrid.forceSequence = true;
		helpTask.cubicGrid.borderSize = mesh->borderSize;
		helpTask.cubicGrid.h = {mesh->h(stage)};
		helpTask.cubicGrid.sizes = {mesh->borderSize};
		
		auto helpMesh = new HelpMesh(helpTask);
		helpMesh->allocate();
		return helpMesh;
	}
	
	
private:
	/// list of conditions that applied in sequence (overwriting previous)
	std::vector<Condition> conditions;
	/// list of inner surfaces
	std::vector<InnerSurface> innerSurfaces;

};


}

#endif // LIBGCM_SPECIALBORDERCONDITIONS_HPP
