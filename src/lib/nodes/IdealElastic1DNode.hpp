#ifndef LIBGCM_IDEALELASTIC1DNODE_HPP
#define LIBGCM_IDEALELASTIC1DNODE_HPP

#include <mpi.h>
#include "lib/gcm_matrices/GcmMatrices.hpp"

#include "lib/nodes/Node.hpp"


namespace gcm {
	class IdealElastic1DContainer {
	public:
		static const int SIZE = 2;
		union {
			real values[SIZE];
			struct {
				real Vx;
				real Sxx;
			};
		};
	};

	class IdealElastic1DNode : public Node<2, IdealElastic1DContainer> {
	public:
		static MPI::Datatype MPI_NODE_TYPE; /// Special type for node for MPI connection

		std::shared_ptr<GcmMatrices<M,1>> matrix; /// pointer to GcmMatrices for the node

		// TODO - can it be moved to Node<> or further?
		template<typename Container>
		IdealElastic1DNode& operator=(const gcm::linal::Matrix<M, 1, Container>& vector) {
			*(static_cast<Node<M, IdealElastic1DContainer>*>(this)) = vector;
			return (*this);
		};
		// TODO (end)
	};
}


#endif // LIBGCM_IDEALELASTIC1DNODE_HPP
