#ifndef LIBGCM_MESHMOVER_HPP
#define LIBGCM_MESHMOVER_HPP

#include <lib/mesh/grid/Cgal2DGrid.hpp>
#include <lib/mesh/DefaultMesh.hpp>


namespace gcm {
	/**
	 * Responsible for mesh motion
	 */
	template<typename TModel, typename TGrid>
	struct MeshMover {
		typedef DefaultMesh<TModel, TGrid>       Mesh;
		
		static void moveMesh(Mesh& mesh, const real timeStep) {
			// by default, mesh is not movable
			SUPPRESS_WUNUSED(mesh);
			SUPPRESS_WUNUSED(timeStep);
		};
		
		USE_AND_INIT_LOGGER("gcm.MeshMover");
	};

	template<typename TModel>
	struct MeshMover<TModel, Cgal2DGrid> {
		typedef DefaultMesh<TModel, Cgal2DGrid>       Mesh;
		typedef typename Mesh::PdeVector              PdeVector;

		static void moveMesh(Mesh& mesh, const real timeStep) {
			LOG_DEBUG("Start mesh motion");
			for (auto& it : mesh) {
				linal::Vector2 dx = {mesh.pde(it).V[0] * timeStep, 
				                     mesh.pde(it).V[1] * timeStep}; // todo Getter?
				mesh.move(it, dx);
			}
		};
		
		USE_AND_INIT_LOGGER("gcm.MeshMover");
	};
}

#endif // LIBGCM_MESHMOVER_HPP
