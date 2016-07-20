#ifndef LIBGCM_MESHMOVER_HPP
#define LIBGCM_MESHMOVER_HPP

#include <lib/mesh/grid/cgal/Cgal2DGrid.hpp>
#include <lib/mesh/DefaultMesh.hpp>


namespace gcm {
/**
 * Responsible for mesh motion
 */
template<typename TModel, typename TGrid, typename TMaterial>
struct MeshMover {
	typedef DefaultMesh<TModel, TGrid, TMaterial> Mesh;

	static void moveMesh(Mesh&, const real) { }

};

template<typename TModel, typename TMaterial>
struct MeshMover<TModel, Cgal2DGrid, TMaterial> {
	typedef DefaultMesh<TModel, Cgal2DGrid, TMaterial> Mesh;
	typedef typename Mesh::PdeVector                   PdeVector;

	static void moveMesh(Mesh& mesh, const real timeStep) {
		if (!mesh.movable) { return; }
		
		USE_AND_INIT_LOGGER("gcm.MeshMover")
		LOG_DEBUG("Start mesh motion");
		for (auto& it : mesh) {
			Real2 dx = { // TODO - use velocity getter indep from dimensions
				mesh.pdeVars(it).velocity(0) * timeStep,
				mesh.pdeVars(it).velocity(1) * timeStep
			};
			mesh.move(it, dx);
		}
	}

};
}

#endif // LIBGCM_MESHMOVER_HPP
