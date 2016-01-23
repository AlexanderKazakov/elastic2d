#ifndef LIBGCM_LINAL_MATRIX_HPP
#define LIBGCM_LINAL_MATRIX_HPP

#include <initializer_list>
#include <iostream>
#include <string.h>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


#include <lib/util/Types.hpp>
#include <lib/util/Assertion.hpp>


namespace gcm {
	namespace linal {
		/**
		 * Default implementation for matrix values container.
		 * @tparam TM First matrix dimension (number of strings).
		 * @tparam TN Second matrix dimension (number of columns).
		 */
		template<int TM, int TN>
		class DefaultMatrixContainer {
		public:
			static const int SIZE = TM * TN; // size of storage in units of gcm::real
			real values[SIZE];
		};

		/**
		 * Generic matrix class.
		 *
		 * @tparam TM First matrix dimension (number of strings).
		 * @tparam TN Second matrix dimension (number of columns).
		 * @tparam Container Container class to hold values.
		 */
		template<int TM, int TN, typename Container=DefaultMatrixContainer<TM, TN>>
		class Matrix : public Container {
		public:
			static const int M = TM; // number of strings
			static const int N = TN; // number of columns

			// TODO
			/*Y& operator=(const Y&) = default;	// default copy semantics
			Y(const Y&) = default;*/
			// TODO - rvalues

			/** Default constructor. */
			Matrix() {
				static_assert(this->SIZE >= TM * TN, "Container must have enough memory to store values");
			};
			/** @param values Values to initialize matrix with, string by string */
			Matrix(std::initializer_list<real> values) {
				this->initialize(values);
			};
			/**
			 * Copy constructor
			 * @param m Matrix to construct from
			 */
			template<typename Container2>
			Matrix(const Matrix<TM, TN, Container2> &m2) {
				(*this) = m2;
			};
			/**
			 * Assignment operator from matrix of equal size and any container
			 * @return reference to modified matrix instance
			 */
			template<typename Container2>
			Matrix<TM, TN, Container> &operator=(const Matrix<TM, TN, Container2> &m2) {
				static_assert(this->SIZE == m2.SIZE, "Containers must have equal size");
				memcpy(this->values, m2.values, sizeof(this->values));
				return *this;
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
			Matrix<TN, TM, Container> transpose() const;
			/** Transposes square matrix (modifying matrix contents) */
			void transposeInplace();

			/**
			 * Inverses matrix. Note this method returns Matrix<TN, TM> (not Matrix<TM, TN>) to make compilation fail
			 * if matrix is not square.
			 * @return inverted matrix.
			 */
			Matrix<TN, TM, Container> invert() const;
			/** Inverses matrix modifying its contents */
			void invertInplace();

			/** @return i-th column. */
			template<typename Container2 = DefaultMatrixContainer<TM, 1>>
			Matrix<TM, 1, Container2> getColumn(const int i) const {
				Matrix<TM, 1, Container2> ans;
				for (int j = 0; j < TM; j++) {
					ans(j) = (*this)(j, i);
				}
				return ans;
			};
			// FIXME - when replacing templation error "no known conversion" appears though templated copy constructor is provided
			/** set i-th column */
			template<typename Container2>
			void setColumn(const int i, const Matrix<TM, 1, Container2> &column) {
				for (int j = 0; j < TM; j++) {
					(*this)(j, i) = column(j);
				}
			};

			/** @return in vector diagonal of this matrix multiplied by matrix B */
			Matrix<TM, 1, DefaultMatrixContainer<TM, 1>> diagonalMultiply(const Matrix<TN, TM> &B) const;

			/** @return trace of the matrix */
			real trace() const;

		protected:
			/** @return values array index of matrix component */
			int getIndex(const int i, const int j) const {
				return i * TN + j;
			};
		};

		template<int TM, int TN, typename Container>
		void Matrix<TM, TN, Container>::initialize(std::initializer_list<real> values) {
			assert_eq(this->SIZE, values.size());
			int i = 0;
			for (auto value: values)
				this->values[i++] = value;
		}

		template<int TM, int TN, typename Container>
		Matrix<TN, TM, Container> Matrix<TM, TN, Container>::transpose() const {
			Matrix<TN, TM, Container> result;

			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TN; j++)
					result(j, i) = (*this)(i, j);

			return result;
		}

		template<int TM, int TN, typename Container>
		void Matrix<TM, TN, Container>::transposeInplace() {
			(*this) = transpose();
		}

		template<int TM, int TN, typename Container>
		Matrix<TN, TM, Container> Matrix<TM, TN, Container>::invert() const {
			Matrix<TN, TM, Container> result;

			gsl_matrix *Z1 = gsl_matrix_alloc(TM, TM);
			gsl_matrix *Z = gsl_matrix_alloc(TM, TM);
			gsl_permutation *perm = gsl_permutation_alloc(TM);
			int k;

			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TM; j++)
					gsl_matrix_set(Z1, i, j, (*this)(i, j));

			if (gsl_linalg_LU_decomp(Z1, perm, &k))
				THROW_INVALID_ARG("gsl_linalg_LU_decomp failed");

			if (gsl_linalg_LU_invert(Z1, perm, Z))
				THROW_INVALID_ARG("gsl_linalg_LU_invert failed");

			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TM; j++)
					result(i, j) = gsl_matrix_get(Z, i, j);

			gsl_permutation_free(perm);
			gsl_matrix_free(Z);
			gsl_matrix_free(Z1);

			return result;
		}

		template<int TM, int TN, typename Container>
		void Matrix<TM, TN, Container>::invertInplace() {
			(*this) = invert();
		}

		template<int TM, int TN, typename Container>
		Matrix<TM, 1, DefaultMatrixContainer<TM, 1>> Matrix<TM, TN, Container>::diagonalMultiply(const Matrix<TN, TM> &B) const {
			assert_eq(TM, TN);
			Matrix<TM, 1, DefaultMatrixContainer<TM, 1>> ans;
			for (int i = 0; i < TM; i++) {
				ans(i) = 0;
				for (int j = 0; j < TM; j++) {
					ans(i) += (*this)(i, j) * B(j, i);
				}
			}
			return ans;
		}

		template<int TM, int TN, typename Container>
		real Matrix<TM, TN, Container>::trace() const {
			assert_eq(TM, TN);
			real trace = 0.0;
			for (int i = 0; i < TN; i++) {
				trace += (*this)(i, i);
			}
			return trace;
		}
	};
};

namespace std {
	template<int TM, int TN, typename Container>
	inline std::ostream &operator<<(std::ostream &os, const gcm::linal::Matrix<TM, TN, Container> &matrix) {

		os << "Matrix:\n";
		for (int i = 0; i < TM; i++) {
			for (int j = 0; j < TN; j++) {
				os << matrix(i, j) << "\t";
			}
			os << "\n";
		}

		return os;
	};
};

#endif // LIBGCM_LINAL_MATRIX_HPP
