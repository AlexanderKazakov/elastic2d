#ifndef LIBGCM_IDEALELASTIC3DNODE_HPP
#define LIBGCM_IDEALELASTIC3DNODE_HPP

#include <mpi.h>
#include "lib/gcm_matrices/GcmMatrices.hpp"

#include "lib/nodes/Node.hpp"


namespace gcm {
	class IdealElastic3DContainer {
	public:
		static const int SIZE = 9;
		static const int V_SIZE = 3;
		static const int S_SIZE = 6;
		union {
			real values[SIZE];
			struct {
				union {
					real V[V_SIZE];
					struct {
						real Vx;
						real Vy;
						real Vz;
					};
				};
				union {
					real S[S_SIZE];
					struct {
						real Sxx;
						real Sxy;
						real Sxz;
						real Syy;
						real Syz;
						real Szz;
					};
				};
			};
		};
	};

	class IdealElastic3DNode : public Node<9, IdealElastic3DContainer> {
	public:
		static MPI::Datatype MPI_NODE_TYPE; // Special type for node for MPI connection

		std::shared_ptr<GcmMatrices<M,3>> matrix; // pointer to GcmMatrices for the node

		real getPressure() const {
			return - (u.Sxx + u.Syy + u.Szz) / 3;
		};

		void setPressure(const real& pressure) {
			u.Sxx = - pressure; u.Syy = - pressure; u.Szz = - pressure;
		};
	};
}


#endif // LIBGCM_IDEALELASTIC3DNODE_HPP
