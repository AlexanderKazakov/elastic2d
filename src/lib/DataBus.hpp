#ifndef LIBGCM_DATABUS_HPP
#define LIBGCM_DATABUS_HPP

#include <mpi.h>

namespace gcm {
	class DataBus {
	public:
		// Special type for node for MPI connection
		static MPI::Datatype MPI_NODE;

		// Creates MPI::Datatype for node - MPI_NODE
		static void createStaticTypes();

	};
}

#endif /* LIBGCM_DATABUS_HPP */
