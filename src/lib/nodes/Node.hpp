#ifndef LIBGCM_NODE_HPP
#define LIBGCM_NODE_HPP

#include <memory>

#include "lib/linal/Vector.hpp"
#include "lib/gcm_matrices/GcmMatrices.hpp"

namespace gcm {

	template<int TM, int TN> class DefaultMatrixContainer;

	template<int M, typename Container=DefaultMatrixContainer<M, 1>>
	class Node : public linal::Vector<M, Container> {
	public:
		typedef linal::Vector<M, Container> Vector;
		typedef linal::Matrix<M, M> Matrix;
	};
}

#endif //LIBGCM_NODE_HPP
