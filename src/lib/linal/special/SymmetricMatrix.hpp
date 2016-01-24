#ifndef LIBGCM_LINAL_SYMMETRICMATRIX_HPP
#define LIBGCM_LINAL_SYMMETRICMATRIX_HPP

#include <lib/linal/Matrix.hpp>


namespace gcm {
	namespace linal {
		/**
		 * Special value container for symmetric matrix.
		 * @tparam M strings/columns number in matrix
		 */
		template<int M>
		class SymmetricMatrixContainer {
		public:
			static const int SIZE = M * (M + 1) / 2; // size of storage in units of gcm::real
			real values[SIZE];
		};

		template<int TM>
		using SymmetricMatrix = Matrix<TM, TM, SymmetricMatrixContainer<TM>>;

		/**
		 * Diagonal matrix class (partial specialization)
		 * @tparam TM strings/columns number in matrix
		 */
		template<int TM>
		class Matrix<TM, TM, SymmetricMatrixContainer<TM>> : public SymmetricMatrixContainer<TM> {
		public:
			static const int M = TM; // number of strings
			static const int N = TM; // number of columns

			/** Default constructor. */
			Matrix() { };
			/** @param values Values to initialize matrix with, string by string */
			Matrix(std::initializer_list<real> values) {
				this->initialize(values);
			};
			/** @param values Values to initialize matrix with, string by string */
			void initialize(std::initializer_list<real> values);

			/** Read-only access to matrix component */
			real operator()(const int i, const int j) const {
				return this->values[getIndex(i, j)];
			};
			/** Read/write access to matrix component */
			real &operator()(const int i, const int j) {
				return this->values[getIndex(i, j)];
			};
			/** Read-only access to vector component */
			real operator()(const int i) const {
				return this->values[i];
			};
			/** Read/write access to vector component */
			real &operator()(const int i) {
				return this->values[i];
			};

			/** @return transposed matrix */
			Matrix<TM, TM, SymmetricMatrixContainer<TM>> transpose() const;
			/** Transposes square matrix (modifying matrix contents) */
			void transposeInplace();

			/** @return inverted matrix */
			Matrix<TM, TM, SymmetricMatrixContainer<TM>> invert() const;
			/** Inverses matrix modifying its contents */
			void invertInplace();

			/** @return trace of the matrix */
			real trace() const;
			
			/** @return values array index of matrix component */
			static int getIndex(const int i, const int j) {
				return (i < j) ? i * N - ((i - 1) * i) / 2 + j - i
				               : j * N - ((j - 1) * j) / 2 + i - j;
			};
		};

		template<int TM>
		void Matrix<TM, TM, SymmetricMatrixContainer<TM>>::initialize(std::initializer_list<real> values) {
			assert_eq(this->SIZE, values.size());
			int i = 0;
			for (auto value: values)
				this->values[i++] = value;
		}

		template<int TM>
		Matrix<TM, TM, SymmetricMatrixContainer<TM>> Matrix<TM, TM, SymmetricMatrixContainer<TM>>::transpose() const {
			return (*this);
		}

		template<int TM>
		void Matrix<TM, TM, SymmetricMatrixContainer<TM>>::transposeInplace() { }

		template<int TM>
		Matrix<TM, TM, SymmetricMatrixContainer<TM>> Matrix<TM, TM, SymmetricMatrixContainer<TM>>::invert() const {
			THROW_UNSUPPORTED("TODO: 1) gsl utils to separate file, 2) for symmetric matrices there are more powerful algorithms");
		}

		template<int TM>
		void Matrix<TM, TM, SymmetricMatrixContainer<TM>>::invertInplace() {
			(*this) = invert();
		}

		template<int TM>
		real Matrix<TM, TM, SymmetricMatrixContainer<TM>>::trace() const {
			real trace = 0.0;
			for (int i = 0; i < TM; i++) {
				trace += (*this)(i, i);
			}
			return trace;
		}
	};
};

#endif // LIBGCM_LINAL_SYMMETRICMATRIX_HPP
