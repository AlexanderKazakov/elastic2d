#ifndef LIBGCM_BORDERCONDITIONS_HPP
#define LIBGCM_BORDERCONDITIONS_HPP

#include <lib/linal/linal.hpp>
#include <lib/mesh/DefaultMesh.hpp>
#include <lib/mesh/CubicGrid.hpp>
#include <lib/mesh/Cgal2DGrid.hpp>


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

		std::map<CUBIC_BORDERS, BorderCondition::T> borderConditions;

	public:
		void initialize(const Task& task);
		void applyBorderConditions(Mesh* mesh) const;
		
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
