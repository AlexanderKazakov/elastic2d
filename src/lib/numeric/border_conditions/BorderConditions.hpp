#ifndef LIBGCM_BORDERCONDITIONS_HPP
#define LIBGCM_BORDERCONDITIONS_HPP

#include <lib/linal/linal.hpp>
#include <lib/mesh/DefaultMesh.hpp>
#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/mesh/grid/Cgal2DGrid.hpp>


namespace gcm {
/**
 * Work with border conditions on different mesh types
 */
template<typename TModel, typename TGrid, typename TMaterial> 
struct BorderConditions;

template<typename TModel, typename TMaterial>
struct BorderConditions<TModel, CubicGrid, TMaterial> {
	typedef DefaultMesh<TModel, CubicGrid, TMaterial>       Mesh;
	typedef typename Mesh::PdeVector                        PdeVector;
	typedef typename Mesh::Iterator                         Iterator;
	typedef typename TModel::PdeVariables                   PdeVariables;

	typedef std::function<real(real)>                       TimeDependency;
	typedef std::map<PhysicalQuantities::T, TimeDependency> Map;
	struct Condition {
		Condition(const std::shared_ptr<Area> area_, const Map& values_) :
			area(area_), values(values_) { }
		std::shared_ptr<Area> area;
		Map values;
	};
	struct InnerSurface {
		InnerSurface(const int direction_, const int index_, const int normal_,
		             std::shared_ptr<Area> area_, const Map& values_) :
			direction(direction_), index(index_), normal(normal_),
			area(area_), values(values_) { }
		int direction;     ///< crossing axis
		int index;         ///< index at crossing axis
		int normal;        ///< -1 or 1 only
		std::shared_ptr<Area> area; ///< area of surface
		Map values;        ///< fixed values
	};

	BorderConditions(const Task& task);
	void beforeStatement(const Statement& statement);
	void applyBorderBeforeStage(Mesh* mesh_, const real timeStep_, const int stage);
	void applyBorderAfterStage(Mesh* mesh_, const real timeStep_, const int stage);

	private:
		/// list of conditions that applied in sequence (overwriting previous)
		std::vector<Condition> conditions;
		/// list of inner surfaces
		std::vector<InnerSurface> innerSurfaces;
		/// temporary values - just to not send between all the functions
		Mesh* mesh;
		real timeStep;
		int direction;
		bool onTheRight;

		// TODO - wtf? replace
		Int3 sizes;    // mesh sizes
		Real3 startR;  // mesh startR
		Real3 lengths; // mesh lengths

		/// auxiliary mesh for fracture calculation
		Mesh* helpMesh;

		void handleSide() const;
		void handleBorderPoint(const Iterator& borderIter, const Map& values) const;
		void allocateHelpMesh();
		void handleInnerSurfacePoint(const Iterator& borderIter, const Map& values,
		                             const int surfaceNormal);

};


template<typename TModel, typename TMaterial>
struct BorderConditions<TModel, Cgal2DGrid, TMaterial> {
	typedef DefaultMesh<TModel, Cgal2DGrid, TMaterial> Mesh;
	typedef typename Mesh::PdeVector                   PdeVector;
	typedef typename Mesh::Iterator                    Iterator;

	BorderConditions(const Task&) { }
	void beforeStatement(const Statement&) { }
	void applyBorderBeforeStage(Mesh*, const real, const int) { }
	void applyBorderAfterStage(Mesh* mesh, const real timeStep, const int stage);

private:
	USE_AND_INIT_LOGGER("gcm.Cgal2DGridBorderConditions")
};


}

#endif // LIBGCM_BORDERCONDITIONS_HPP
