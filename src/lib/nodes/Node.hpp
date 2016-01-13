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
	class Node {
	public:
		typedef linal::Vector<TM, Container> Vector;
		typedef linal::Matrix<TM, TM> Matrix;

		static const int M = Vector::M;

		/** The nodal data (PDE variables) */
		Vector u;
	};
}

#endif //LIBGCM_NODE_HPP
