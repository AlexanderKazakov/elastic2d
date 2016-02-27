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
		typedef typename Mesh::PartIterator          PartIterator;
		typedef typename TModel::PdeVariables        PdeVariables;

		typedef std::function<real(real)>                       TimeDependency;
		typedef std::map<PhysicalQuantities::T, TimeDependency> Map;
		struct Condition {
			Condition(const std::shared_ptr<Area> area_, const Map& values_) :
					area(area_), values(values_) { };
			std::shared_ptr<Area> area;
			Map values;
		};
		
		void initialize(const Task& task);
		void applyBorderConditions(Mesh* mesh_, const real time_);
	
	private:
		// list of conditions that applied in sequence (overwriting previous)
		std::vector<Condition> conditions;
		// temporary values - just to not send between all the functions
		Mesh* mesh;
		real time;
		int direction;
		bool onTheRight;
		
		void handleSide() const;
		void handlePoint(const PartIterator& borderIter, const Map& values) const;
	};

	template<typename TModel>
	struct BorderConditions<TModel, Cgal2DGrid> {
		typedef DefaultMesh<TModel, Cgal2DGrid>      Mesh;
		typedef typename Mesh::Matrix                Matrix;
		typedef typename Mesh::PdeVector             PdeVector;
		typedef typename Mesh::Iterator              Iterator;

		void initialize(const Task&) { };
		void applyBorderConditions(Mesh*, const real) const { };
	};
}

#endif // LIBGCM_BORDERCONDITIONS_HPP
