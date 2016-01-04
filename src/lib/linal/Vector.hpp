#ifndef GCM_VECTOR_HPP
#define GCM_VECTOR_HPP

#include <iostream>

#include "lib/linal/Matrix.hpp"

namespace gcm {
	namespace linal {
		template<int TM, typename Container>
		class Vector : public Matrix<TM, 1, Container> {
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

		template<int TM, typename Container>
		inline
		gcm::real Vector<TM, Container>::operator()(const int i) const
		{
			return this->values[i];
		};

		template<int TM, typename Container>
		inline
		gcm::real& Vector<TM, Container>::operator()(const int i)
		{
			return this->values[i];
		};
	}
};

namespace std {
	template<int TM>
	inline std::ostream &operator<<(std::ostream &os, const gcm::linal::Vector<TM> &vector) {

		os << "Vector:\n";
		for (int i = 0; i < TM; i++) {
			os << vector(i) << "\n";
		}

		return os;
	};
};

#endif //GCM_VECTOR_HPP