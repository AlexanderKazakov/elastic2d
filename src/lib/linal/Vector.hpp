#ifndef LIBGCM_LINAL_VECTOR_HPP
#define LIBGCM_LINAL_VECTOR_HPP

#include <lib/linal/Matrix.hpp>

namespace gcm {
	namespace linal {
		/**
		 * Generic vector - just a matrix with one column
		 */
		template<int TM, typename Container = DefaultMatrixContainer<TM, 1>>
		using Vector = Matrix<TM, 1, Container>;
	};
};

#endif // LIBGCM_LINAL_VECTOR_HPP
