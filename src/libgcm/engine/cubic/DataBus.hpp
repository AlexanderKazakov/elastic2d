#ifndef LIBGCM_CUBIC_DATABUS_HPP
#define LIBGCM_CUBIC_DATABUS_HPP


namespace gcm {
namespace cubic {

template<typename Mesh>
struct DataBus {
	
	static void exchangeNodesWithNeighbors(Mesh* mesh) {
		if (Mpi::ForceSequence() || Mpi::Size() == 1) { return; }
		exchangeSomethingWithNeighbors(mesh, mesh->pdeVariables);
	}
	
	template<typename Smth>
	static void exchangeSomethingWithNeighbors(Mesh* mesh, std::vector<Smth>& vec) {
		
		size_t bufferSize = (size_t) (mesh->borderSize * mesh->indexMaker(0));
		size_t size = vec.size();
		
		if (Mpi::Rank() == 0) {
			MPI_Sendrecv(&(vec[size - 2 * bufferSize]), (int) (sizeof(Smth) * bufferSize),
						 MPI_BYTE,
						 Mpi::Rank() + 1, 1,
						 &(vec[size - bufferSize]), (int) (sizeof(Smth) * bufferSize), MPI_BYTE,
						 Mpi::Rank() + 1, 1,
						 MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		
		} else if (Mpi::Rank() == Mpi::Size() - 1) {
			MPI_Sendrecv(&(vec[bufferSize]), (int) (sizeof(Smth) * bufferSize), MPI_BYTE,
						 Mpi::Rank() - 1, 1,
						 &(vec[0]), (int) (sizeof(Smth) * bufferSize), MPI_BYTE,
						 Mpi::Rank() - 1, 1,
						 MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		
		} else {
			MPI_Sendrecv(&(vec[size - 2 * bufferSize]), (int) (sizeof(Smth) * bufferSize),
						 MPI_BYTE,
						 Mpi::Rank() + 1, 1,
						 &(vec[size - bufferSize]), (int) (sizeof(Smth) * bufferSize), MPI_BYTE,
						 Mpi::Rank() + 1, 1,
						 MPI::COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv(&(vec[bufferSize]), (int) (sizeof(Smth) * bufferSize), MPI_BYTE,
						 Mpi::Rank() - 1, 1,
						 &(vec[0]), (int) (sizeof(Smth) * bufferSize), MPI_BYTE,
						 Mpi::Rank() - 1, 1,
						 MPI::COMM_WORLD, MPI_STATUS_IGNORE);
		}
		
	}
};


} // namespace cubic
} // namespace gcm


#endif // LIBGCM_CUBIC_DATABUS_HPP
