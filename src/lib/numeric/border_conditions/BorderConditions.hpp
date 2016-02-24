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

		void initialize(const Task& task);
		void applyBorderConditions(Mesh* mesh) const;
	private:
		std::map<DIRECTION, 
			std::pair<BorderCondition::T, BorderCondition::T>> borderConditions;

		void handleSide(Mesh* mesh, DIRECTION direction, 
		                const bool onTheRight) const;
		void handleSlice(Mesh* mesh, PartIterator realIter,
		                             PartIterator virtIter) const;
	};

	template<typename TModel>
	struct BorderConditions<TModel, Cgal2DGrid> {
		typedef DefaultMesh<TModel, Cgal2DGrid>      Mesh;
		typedef typename Mesh::Matrix                Matrix;
		typedef typename Mesh::PdeVector             PdeVector;
		typedef typename Mesh::Iterator              Iterator;

		void initialize(const Task& task) {
			SUPPRESS_WUNUSED(task);
		};
		void applyBorderConditions(Mesh* mesh) const {
			SUPPRESS_WUNUSED(mesh);
		};
	};
}

#endif // LIBGCM_BORDERCONDITIONS_HPP
