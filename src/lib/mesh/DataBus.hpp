#ifndef LIBGCM_DATABUS_HPP
#define LIBGCM_DATABUS_HPP


#include <lib/mesh/grid/CubicGrid.hpp>
#include <lib/mesh/grid/Cgal2DGrid.hpp>
#include <lib/mesh/DefaultMesh.hpp>

namespace gcm {
	/**
	 * Responsible for mpi connection
	 */
	template<typename TModel, typename TGrid, typename TMaterial>
	struct DataBus;

	template<typename TModel, typename TMaterial>
	struct DataBus<TModel, CubicGrid, TMaterial> {
		typedef DefaultMesh<TModel, CubicGrid, TMaterial>       Mesh;

		static void exchangeNodesWithNeighbors(Mesh* mesh) {
			LOG_DEBUG("Start data exchange with neighbor cores");
			if (mesh->numberOfWorkers == 1) return;

			exchangeSomethingWithNeighbors(mesh, mesh->pdeVectors);
			if (TModel::InternalOde::NonTrivial) {
				exchangeSomethingWithNeighbors(mesh, mesh->odeValues);
			}
		}

		template<typename Smth>
		static void exchangeSomethingWithNeighbors(Mesh* mesh, std::vector<Smth>& vec);
		
		USE_AND_INIT_LOGGER("gcm.DataBus")
	};

	template<typename TModel, typename TMaterial>
	struct DataBus<TModel, Cgal2DGrid, TMaterial> {
		typedef DefaultMesh<TModel, Cgal2DGrid, TMaterial>      Mesh;

		static void exchangeNodesWithNeighbors(Mesh*) { }

		USE_AND_INIT_LOGGER("gcm.DataBus")
	};



	template<typename TModel, typename TMaterial>
	template<typename Smth>
	void DataBus<TModel, CubicGrid, TMaterial>::
	exchangeSomethingWithNeighbors(Mesh* mesh, std::vector<Smth>& vec) {

		int bufferSize = mesh->borderSize * mesh->indexMaker(0);
		size_t size = vec.size();

		if (mesh->rank == 0) {
			MPI_Sendrecv(&(vec[size - 2 * bufferSize]), (int) sizeof(Smth) * bufferSize, MPI_BYTE, mesh->rank + 1, 1,
						 &(vec[size - bufferSize]), (int) sizeof(Smth) * bufferSize, MPI_BYTE, mesh->rank + 1, 1,
						 MPI::COMM_WORLD, MPI_STATUS_IGNORE);

		} else if (mesh->rank == mesh->numberOfWorkers - 1) {
			MPI_Sendrecv(&(vec[(size_t)bufferSize]), (int) sizeof(Smth) * bufferSize, MPI_BYTE, mesh->rank - 1, 1,
						 &(vec[0]), (int) sizeof(Smth) * bufferSize, MPI_BYTE, mesh->rank - 1, 1,
						 MPI::COMM_WORLD, MPI_STATUS_IGNORE);

		} else {
			MPI_Sendrecv(&(vec[size - 2 * bufferSize]), (int) sizeof(Smth) * bufferSize, MPI_BYTE, mesh->rank + 1, 1,
						 &(vec[size - bufferSize]), (int) sizeof(Smth) * bufferSize, MPI_BYTE, mesh->rank + 1, 1,
						 MPI::COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv(&(vec[(size_t)bufferSize]), (int) sizeof(Smth) * bufferSize, MPI_BYTE, mesh->rank - 1, 1,
						 &(vec[0]), (int) sizeof(Smth) * bufferSize, MPI_BYTE, mesh->rank - 1, 1,
						 MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

}

#endif // LIBGCM_DATABUS_HPP
