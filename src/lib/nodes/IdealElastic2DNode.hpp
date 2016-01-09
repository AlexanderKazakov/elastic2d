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
			return - (Sxx + Syy) / 2;
		}

		void setPressure(const real& pressure) {
			Sxx = - pressure; Syy = - pressure;
		}

		/*static struct {*/
		static const std::map<const std::string /* name */, const std::pair<const int /* index */, const int /* size */>> VECTORS;
		static const std::map<const std::string /* name */, const int /* index */> SCALARS;
		/*} Map; // for snapshotters*/


		// TODO - can it be moved to Node<> or further?
		template<typename Container>
		IdealElastic2DNode& operator=(const gcm::linal::Matrix<M, 1, Container>& vector) {
			*(static_cast<Node<M, IdealElastic2DContainer>*>(this)) = vector;
			return (*this);
		};
		// TODO (end)
	};
}


#endif // LIBGCM_IDEALELASTIC2DNODE_HPP
