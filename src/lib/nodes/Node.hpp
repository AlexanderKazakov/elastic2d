#ifndef LIBGCM_NODE_HPP
#define LIBGCM_NODE_HPP

#include <memory>

#include "lib/linal/Vector.hpp"
#include "PDEMatrices2D.hpp"

namespace gcm {
	template<int M, typename Container=DefaultMatrixContainer<M, 1>>
	class Node : public linal::Vector<M, Container> {
	public:
		typedef linal::Vector<M, Container> Vector;
		std::shared_ptr<PDEMatrices2D> matrix;
	};
}

#endif //LIBGCM_NODE_HPP
