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

		template<typename MatrixContainer>
		Node<TM, Container>& operator=(const gcm::linal::Matrix<TM, 1, MatrixContainer>& vector) {
			*(static_cast<linal::Matrix<TM, 1, Container>*>(this)) = vector;
			return (*this);
		};

	};
}

#endif //LIBGCM_NODE_HPP
