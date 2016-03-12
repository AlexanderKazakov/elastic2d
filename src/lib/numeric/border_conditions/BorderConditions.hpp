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
	template<typename TModel, typename TGrid>
	struct BorderConditions { };

	template<typename TModel>
	struct BorderConditions<TModel, CubicGrid> {
		typedef DefaultMesh<TModel, CubicGrid>       Mesh;
		typedef typename Mesh::Matrix                Matrix;
		typedef typename Mesh::PdeVector             PdeVector;
		typedef typename Mesh::Iterator              Iterator;
		typedef typename TModel::PdeVariables        PdeVariables;

		typedef std::function<real(real)>                       TimeDependency;
		typedef std::map<PhysicalQuantities::T, TimeDependency> Map;
		struct Condition {
			Condition(const std::shared_ptr<Area> area_, const Map& values_) :
					area(area_), values(values_) { }
			std::shared_ptr<Area> area;
			Map values;
		};
		struct Fracture {
			Fracture(const int direction_, const int index_, const int normal_,
			         std::shared_ptr<Area> area_, const Map& values_) :
					direction(direction_), index(index_), normal(normal_), 
					area(area_), values(values_) { }
			int direction; // crossing axis
			int index; // index at crossing axis
			int normal; // -1 or 1 only
			std::shared_ptr<Area> area; 
			Map values; // fixed values
		};
		
		void initialize(const Task& task);
		void beforeStatement(const Statement &statement);
		void applyBorderBeforeStage(Mesh* mesh_, const real currentTime_, 
		                            const real timeStep_, const int stage);
		void applyBorderAfterStage(Mesh* mesh_, const real currentTime_, 
		                           const real timeStep_, const int stage);
	
	private:
		// list of conditions that applied in sequence (overwriting previous)
		std::vector<Condition> conditions;
		// list of inner fractures
		std::vector<Fracture> fractures;
		// temporary values - just to not send between all the functions
		Mesh* mesh;
		linal::Int3 sizes; // mesh sizes
		linal::Vector<3> startR; // mesh startR
		linal::Vector<3> lengths; // mesh lengths
		real currentTime;
		real timeStep;
		int direction;
		bool onTheRight;
		// auxiliary mesh for fracture calculation
		Mesh* helpMesh;
		
		void handleSide() const;
		void handleBorderPoint(const Iterator& borderIter, const Map& values) const;
		void allocateHelpMesh();
		void handleFracturePoint(const Iterator& borderIter, const Map& values,
		                         const int fracNormal);
	};

	template<typename TModel>
	struct BorderConditions<TModel, Cgal2DGrid> {
		typedef DefaultMesh<TModel, Cgal2DGrid>      Mesh;
		typedef typename Mesh::Matrix                Matrix;
		typedef typename Mesh::PdeVector             PdeVector;
		typedef typename Mesh::Iterator              Iterator;

		void initialize(const Task&) { }
		void beforeStatement(const Statement&) { }
		void applyBorderBeforeStage(Mesh*, const real, const real, const int) { }
		void applyBorderAfterStage(Mesh*, const real, const real, const int) { }
	};
}

#endif // LIBGCM_BORDERCONDITIONS_HPP
