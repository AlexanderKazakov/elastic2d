#ifndef LIBGCM_NODE_HPP
#define LIBGCM_NODE_HPP

#include <memory>

#include "lib/linal/Vector.hpp"
#include "lib/gcm_matrices/GcmMatrices.hpp"

namespace gcm {

	template<int TM, int TN> class DefaultMatrixContainer;

	template<int TM, typename Container=DefaultMatrixContainer<TM, 1>>
	class Node : public linal::Vector<TM, Container> {
	public:
		typedef linal::Vector<TM, Container> Vector;
		typedef linal::Matrix<TM, TM> Matrix;

	};
}

#endif //LIBGCM_NODE_HPP
