#ifndef LIBGCM_IDEALELASTIC1DNODE_HPP
#define LIBGCM_IDEALELASTIC1DNODE_HPP

#include <mpi.h>

#include "lib/gcm_matrices/GcmMatrices.hpp"
#include "lib/nodes/Node.hpp"


namespace gcm {
	class IdealElastic1DContainer {
	public:
		static const int SIZE = 2;
		static const int V_SIZE = 1;
		static const int S_SIZE = 1;
		union {
			real values[SIZE];
			struct {
				union {
					real V[V_SIZE];
					real Vx;
				};
				union {
					real S[S_SIZE];
					real Sxx;
				};
			};
		};
	};

	class IdealElastic1DNode : public Node<2, IdealElastic1DContainer> {
	public:
		static MPI::Datatype MPI_NODE_TYPE; // Special type for node for MPI connection

		std::shared_ptr<GcmMatrices<M,1>> matrix; // pointer to GcmMatrices for the node

		template<typename Container>
		IdealElastic1DNode& operator=(const linal::Vector<M, Container>& vector) {
			*(static_cast<Node<M, IdealElastic1DContainer>*>(this)) = vector;
			return (*this);
		};

		real getPressure() const {
			return - Sxx;
		};

		void setPressure(const real& pressure) {
			Sxx = - pressure;
		};

	};
}


#endif // LIBGCM_IDEALELASTIC1DNODE_HPP
