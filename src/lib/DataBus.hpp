#ifndef ELASTIC2D_DATABUS_HPP
#define ELASTIC2D_DATABUS_HPP

#include <mpi.h>


class DataBus {
public:
	DataBus();

	// Special type for node for MPI connection
	static MPI::Datatype MPI_NODE;
	// Creates MPI::Datatype for node - MPI_NODE
	static void createStaticTypes();

private:

	// Rank of current core
	int rank;
	// Number of cores
	int numberOfWorkers;


};

#endif /* ELASTIC2D_DATABUS_HPP */
