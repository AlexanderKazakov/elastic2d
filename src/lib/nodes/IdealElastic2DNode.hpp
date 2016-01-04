#ifndef LIBGCM_IDEALELASTIC2DNODE_HPP
#define LIBGCM_IDEALELASTIC2DNODE_HPP

#include <mpi.h>
#include "lib/gcm_matrices/GcmMatrices.hpp"

#include "lib/nodes/Node.hpp"


namespace gcm {
	class IdealElastic2DContainer {
	public:
		static const int SIZE = 5;
		union {
			real values[SIZE];
			struct {
				real Vx;
				real Vy;
				real Sxx;
				real Sxy;
				real Syy;
			};
		};
	};

	class IdealElastic2DNode : public Node<5, IdealElastic2DContainer> {
	public:
		static MPI::Datatype MPI_NODE_TYPE; /// Special type for node for MPI connection

		std::shared_ptr<GcmMatrices<5,2>> matrix; /// pointer to GcmMatrices for the node
	};
}


#endif // LIBGCM_IDEALELASTIC2DNODE_HPP
