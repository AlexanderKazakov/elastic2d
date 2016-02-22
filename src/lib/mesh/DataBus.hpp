#ifndef LIBGCM_DATABUS_HPP
#define LIBGCM_DATABUS_HPP


#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/mesh/grid/Cgal2DGrid.hpp>
#include <lib/mesh/DefaultMesh.hpp>

namespace gcm {
	/**
	 * Responsible for mpi connection
	 */
	template<typename TModel, typename TGrid>
	struct DataBus {
		USE_AND_INIT_LOGGER("gcm.DataBus");
	};

	template<typename TModel>
	struct DataBus<TModel, CubicGrid> {
		typedef DefaultMesh<TModel, CubicGrid>       Mesh;
		typedef typename Mesh::PdeVector             PdeVector;

		static void exchangeNodesWithNeighbors(Mesh* mesh) {
			LOG_DEBUG("Start data exchange with neighbor cores");
			if (mesh->numberOfWorkers == 1) return;
		
			int bufferSize = mesh->accuracyOrder
			                 * (mesh->getSizes()(1) + 2 * mesh->accuracyOrder)
			                 * (mesh->getSizes()(2) + 2 * mesh->accuracyOrder);
			size_t size = mesh->pdeVectors.size();
		
			if (mesh->rank == 0) {
				MPI_Sendrecv(&(mesh->pdeVectors[size - 2 * bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, mesh->rank + 1, 1,
							 &(mesh->pdeVectors[size - bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, mesh->rank + 1, 1,
							 MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		
			} else if (mesh->rank == mesh->numberOfWorkers - 1) {
				MPI_Sendrecv(&(mesh->pdeVectors[(size_t)bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, mesh->rank - 1, 1,
							 &(mesh->pdeVectors[0]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, mesh->rank - 1, 1,
							 MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		
			} else {
				MPI_Sendrecv(&(mesh->pdeVectors[size - 2 * bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, mesh->rank + 1, 1,
							 &(mesh->pdeVectors[size - bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, mesh->rank + 1, 1,
							 MPI::COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(mesh->pdeVectors[(size_t)bufferSize]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, mesh->rank - 1, 1,
							 &(mesh->pdeVectors[0]), (int) sizeof(PdeVector) * bufferSize, MPI_BYTE, mesh->rank - 1, 1,
							 MPI::COMM_WORLD, MPI_STATUS_IGNORE);
			}
		};
		
		USE_AND_INIT_LOGGER("gcm.DataBus");
	};

	template<typename TModel>
	struct DataBus<TModel, Cgal2DGrid> {
		typedef DefaultMesh<TModel, Cgal2DGrid>      Mesh;
		typedef typename Mesh::PdeVector             PdeVector;

		static void exchangeNodesWithNeighbors(Mesh* mesh) {
			SUPPRESS_WUNUSED(mesh); // TODO
		};

		USE_AND_INIT_LOGGER("gcm.DataBus");
	};

}

#endif // LIBGCM_DATABUS_HPP
