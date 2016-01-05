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

			// TODO - why it's not inherited from Matrix?
			Vector(std::initializer_list<real> values);

			Vector() {};

			template<typename Container2>
			Vector(const Vector<TM, Container2>& v2) {
				(*this) = v2;
			};
			template<typename Container2>
			Vector(const Matrix<TM, 1, Container2>& v2) {
				(*this) = v2;
			};

			template<typename Container2>
			Vector<TM, Container>& operator=(const Vector<TM, Container2>& v2)
			{
				for (int i = 0; i < TM; i++)
					(*this)(i) = v2(i);

				return *this;
			};
			template<typename Container2>
			Vector<TM, Container>& operator=(const Matrix<TM, 1, Container2>& v2)
			{
				for (int i = 0; i < TM; i++)
					(*this)(i) = v2(i, 0);

				return *this;
			};
			// TODO (end)
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

		template<int TM, typename Container>
		Vector<TM, Container>::Vector(std::initializer_list<real> values): Vector()
		{
			int i = 0;
			for (auto value: values)
				this->values[i++] = value;
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