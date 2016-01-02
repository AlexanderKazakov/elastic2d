#ifndef GCM_VECTOR_HPP
#define GCM_VECTOR_HPP

#include <iostream>

#include "lib/linal/Matrix.hpp"

namespace gcm {
	namespace linal {
		template<int M, typename Container=DefaultMatrixContainer<M, 1>>
		class Vector : public Matrix<M, 1, Container> {
		public:
			/**
			 * Returns vector component.
			 *
			 * @param i index.
			 *
			 * @return Corresponding vector component.
			 */
			real operator()(const int i) const;

			/**
			 * Returns reference to vector component, used to modify vector content.
			 *
			 * @param i index.
			 *
			 * @return Reference to corresponding vector component.
			 */
			real &operator()(const int i);
		};

		template<int M, typename Container>
		inline
		gcm::real Vector<M, Container>::operator()(const int i) const
		{
			return this->values[getIndex(i, 1)];
		};

		template<int M, typename Container>
		inline
		gcm::real& Vector<M, Container>::operator()(const int i)
		{
			return this->values[getIndex(i, 1)];
		};
	}
};

namespace std {
	template<int M>
	inline std::ostream &operator<<(std::ostream &os, const gcm::linal::Vector<M> &vector) {

		os << "Vector:\n";
		for (int i = 0; i < M; i++) {
			os << vector(i) << "\n";
		}

		return os;
	};
};

#endif //GCM_VECTOR_HPP