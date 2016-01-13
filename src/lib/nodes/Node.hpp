#ifndef LIBGCM_NODE_HPP
#define LIBGCM_NODE_HPP

#include <memory>

#include "lib/linal/Linal.hpp"
#include "lib/gcm_matrices/GcmMatrices.hpp"

namespace gcm {
	template<int TM, int TN> class DefaultMatrixContainer;

	/**
	 * Base class for all nodes for structured grids (it doesn't contain coordinates)
	 * @tparam TM number of PDE variables
	 * @tparam Container container used to store node values
	 */
	template<int TM, typename Container=DefaultMatrixContainer<TM, 1>>
	class Node : public linal::Vector<TM, Container> {
	public:
		typedef linal::Vector<TM, Container> Vector;
		typedef linal::Matrix<TM, TM> Matrix;

		// FIXME - either replace inheritance from Vector to ownership, or keep this shit
		template<typename Container2>
		Node<TM, Container>& operator=(const linal::Vector<TM, Container2>& vector) {
			*(static_cast<linal::Matrix<TM, 1, Container>*>(this)) = vector;
			return (*this);
		};

	};
}

#endif //LIBGCM_NODE_HPP
