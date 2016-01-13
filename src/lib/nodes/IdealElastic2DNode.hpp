#ifndef LIBGCM_IDEALELASTIC2DNODE_HPP
#define LIBGCM_IDEALELASTIC2DNODE_HPP

#include <mpi.h>
#include "lib/gcm_matrices/GcmMatrices.hpp"

#include "lib/nodes/Node.hpp"


namespace gcm {
	class IdealElastic2DContainer {
	public:
		static const int SIZE = 5;
		static const int V_SIZE = 2;
		static const int S_SIZE = 3;
		union {
			real values[SIZE];
			struct {
				union {
					real V[V_SIZE];
					struct {
						real Vx;
						real Vy;
					};
				};
				union {
					real S[S_SIZE];
					struct {
						real Sxx;
						real Sxy;
						real Syy;
					};
				};
			};
		};
	};

	// TODO - why size of the node is 56 bytes?
	class IdealElastic2DNode : public Node<5, IdealElastic2DContainer> {
	public:
		static MPI::Datatype MPI_NODE_TYPE; // Special type for node for MPI connection

		std::shared_ptr<GcmMatrices<M,2>> matrix; // pointer to GcmMatrices for the node

		real getPressure() const {
			return - (u.Sxx + u.Syy) / 2;
		}

		void setPressure(const real& pressure) {
			u.Sxx = - pressure; u.Syy = - pressure;
		}
	};
}


#endif // LIBGCM_IDEALELASTIC2DNODE_HPP
