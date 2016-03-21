#ifndef LIBGCM_LINAL_VECTORINT_HPP
#define LIBGCM_LINAL_VECTORINT_HPP

#include <lib/linal/LinalRoutines.hpp>

namespace gcm {
	namespace linal {
		template<int M>
		class VectorIntContainer {
		public:
			static const int SIZE = M; // size of storage in units of int
			int values[SIZE];
		};

		template<int TM>
		using VectorInt = Matrix<TM, 1, VectorIntContainer<TM>>;

		typedef linal::VectorInt<1> Int1;
		typedef linal::VectorInt<2> Int2;
		typedef linal::VectorInt<3> Int3;

		/**
		 * Vector of integers
		 * @tparam TM number of components
		 */
		template<int TM>
		class Matrix<TM, 1, VectorIntContainer<TM>> : public VectorIntContainer<TM> {
		public:
			static const int M = TM; // number of strings
			static const int N = 1; // number of columns

			/** Default constructor. */
			Matrix() { }
			/** @param values Values to initialize matrix with, string by string */
			Matrix(std::initializer_list<int> list) {
				this->initialize(list);
			}
			/**
			 * Copy constructor
			 * @param m Matrix to construct from
			 */
			template<typename Container2>
			Matrix(const Matrix<TM, 1, Container2> &m2) {
				(*this) = m2;
			}
			/**
			 * Assignment operator from matrix of equal size and any container
			 * @return reference to modified matrix instance
			 */
			template<typename Container2>
			Matrix<TM, 1, VectorIntContainer<TM>>& operator=(const Matrix<TM, 1, Container2> &m2) {
				static_assert(this->SIZE == m2.SIZE, "Containers must have equal size");
				for (int i = 0; i < M; i++) {
					(*this)(i) = (int) m2(i);
				}
				return *this;
			}
			/** @param values Values to initialize matrix with, string by string */
			void initialize(std::initializer_list<int> list);

			/** Read-only access to matrix component */
			int operator()(const int i, const int j) const {
				return this->values[getIndex(i, j)];
			}
			/** Read/write access to matrix component */
			int &operator()(const int i, const int j) {
				return this->values[getIndex(i, j)];
			}
			/** Read-only access to vector component */
			int operator()(const int i) const {
				return this->values[i];
			}
			/** Read/write access to vector component */
			int &operator()(const int i) {
				return this->values[i];
			}

			bool operator==(const Matrix<TM, 1, VectorIntContainer<TM>>& m2) const {
				for (int i = 0; i < TM; i++)
					if ((*this)(i) != m2(i))
						return false;
				return true;
			}
			bool operator!=(const Matrix<TM, 1, VectorIntContainer<TM>>& m2) const {
				return !((*this) == m2);
			}

		protected:
			/** @return values array index of matrix component */
			int getIndex(const int i, const int j) const {
				return i * N + j;
			}
		};

		template<int TM>
		void Matrix<TM, 1, VectorIntContainer<TM>>::initialize(std::initializer_list<int> list) {
			assert_eq(this->SIZE, list.size());
			int i = 0;
			for (auto value: list)
				this->values[i++] = value;
		}

		/** @return multiplication of all elements */
		template<int TM>
		long int directProduct(const Matrix<TM, 1, VectorIntContainer<TM>>& m) {
			long int ans = 1;
			for (int i = 0; i < TM; i++) {
				ans *= m(i);
			}
			return ans;
		}
	}
}

#endif // LIBGCM_LINAL_VECTORINT_HPP
