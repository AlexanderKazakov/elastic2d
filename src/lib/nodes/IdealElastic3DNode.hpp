#ifndef LIBGCM_IDEALELASTIC3DNODE_HPP
#define LIBGCM_IDEALELASTIC3DNODE_HPP

#include <mpi.h>
#include "lib/gcm_matrices/GcmMatrices.hpp"

#include "lib/nodes/Node.hpp"


namespace gcm {
	class IdealElastic3DContainer {
	public:
		static const int SIZE = 9;
		union {
			real values[SIZE];
			struct {
				real Vx;
				real Vy;
				real Vz;
				real Sxx;
				real Sxy;
				real Sxz;
				real Syy;
				real Syz;
				real Szz;
			};
		};
	};

	class IdealElastic3DNode : public Node<9, IdealElastic3DContainer> {
	public:
		static MPI::Datatype MPI_NODE_TYPE; /// Special type for node for MPI connection

		std::shared_ptr<GcmMatrices<M,3>> matrix; /// pointer to GcmMatrices for the node

		// TODO - can it be moved to Node<> or further?
		template<typename Container>
		IdealElastic3DNode& operator=(const gcm::linal::Matrix<M, 1, Container>& vector) {
			*(static_cast<Node<M, IdealElastic3DContainer>*>(this)) = vector;
			return (*this);
		};
		// TODO (end)
	};
}


#endif // LIBGCM_IDEALELASTIC3DNODE_HPP
