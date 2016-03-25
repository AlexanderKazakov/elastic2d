#ifndef LIBGCM_LINAL_LINALROUTINES_HPP
#define LIBGCM_LINAL_LINALROUTINES_HPP

#include <limits>

#include <lib/linal/Matrix.hpp>
#include <lib/linal/Vector.hpp>

namespace gcm {
	namespace linal {

		/** @return matrix diagonal into vector */
		template<int TM, typename Container>
		Vector<TM> diag(const Matrix<TM, TM, Container>& m) {
			Vector<TM> ans;
			for (int i = 0; i < TM; i++) {
				ans(i) = m(i, i);
			}
			return ans;
		}

		/** Set all the container's values to zero */
		template<class TContainer>
		TContainer& clear(TContainer &container) {
			memset(container.values, 0, sizeof(TContainer));
			return container;
		}

		/**
		 * Computes negative of matrix B
		 *
		 * @tparam TM First matrix dimension
		 * @tparam TN Second matrix dimension
		 * @tparam Container Matrix container type
		 * @param m Matrix to compute negative for.
		 *
		 * @return Negative of matrix
		 */
		template<int TM, int TN, typename Container>
		Matrix<TM, TN, Container> operator-(const Matrix<TM, TN, Container> &m) {
			Matrix<TM, TN, Container> result;

			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TN; j++)
					result(i, j) = -m(i, j);

			return result;
		}

		/**
		 * Computes summ of two matrices. Generic implementation.
		 *
		 * @tparam TM First matrix dimension.
		 * @tparam TN Second matrix dimension.
		 * @tparam Container1 Container type for first summ item.
		 * @tparam Container2 Container type for second  summ item.
		 * @tparam Container3 Container type for result.
		 * @param m1 First summ item.
		 * @param m2 Second summ item.
		 *
		 * @return Matrix summ (m1+m2).
		 */
		template<int TM, int TN, typename Container1, typename Container2, typename Container3>
		Matrix<TM, TN, Container3> operator+(const Matrix<TM, TN, Container1> &m1,
		                                     const Matrix<TM, TN, Container2> &m2) {
			Matrix<TM, TN, Container3> result;

			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TN; j++)
					result(i, j) = m1(i, j) + m2(i, j);

			return result;
		}

		/**
		 * Computes summ of two matrices. Most useful specialization, made in assumption that result matrix should
		 * use default matrix container.
		 *
		 * @tparam TM First matrix dimension.
		 * @tparam TN Second matrix dimension.
		 * @tparam Container1 Container type for first summ item.
		 * @tparam Container2 Container type for second  summ item.
		 * @param m1 First summ item.
		 * @param m2 Second summ item.
		 *
		 * @return Matrix summ (m1+m2).
		 */
		template<int TM, int TN, typename Container1, typename Container2>
		Matrix<TM, TN, DefaultMatrixContainer<TM, TN>> operator+(const Matrix<TM, TN, Container1> &m1,
		                                                         const Matrix<TM, TN, Container2> &m2) {
			return operator+<TM, TN, Container1, Container2, DefaultMatrixContainer<TM, TN>>(m1, m2);
		}

		/**
		 * Computes summ of two matrices. Specialized implementation, made in assumption that all matrices
		 * use default matrix container.
		 *
		 * @tparam TM First matrix dimension.
		 * @tparam TN Second matrix dimension.
		 * @tparam Container Container type for all matrices.
		 * @param m1 First summ item.
		 * @param m2 Second summ item.
		 *
		 * @return Matrix summ (m1+m2).
		 */
		template<int TM, int TN, typename Container>
		Matrix<TM, TN, Container> operator+(const Matrix<TM, TN, Container> &m1, const Matrix<TM, TN, Container> &m2) {
			return operator+<TM, TN, Container, Container, Container>(m1, m2);
		}

		/**
		 * Computes difference of two matrices. Generic implementation.
		 *
		 * @tparam TM First matrix dimension.
		 * @tparam TN Second matrix dimension.
		 * @tparam Container1 Container type for first difference item.
		 * @tparam Container2 Container type for second  difference item.
		 * @tparam Container3 Container type for result.
		 * @param m1 First difference item.
		 * @param m2 Second difference item.
		 *
		 * @return Matrix difference (m1-m2).
		 */
		template<int TM, int TN, typename Container1, typename Container2, typename Container3>
		Matrix<TM, TN, Container3> operator-(const Matrix<TM, TN, Container1> &m1,
		                                     const Matrix<TM, TN, Container2> &m2) {
			Matrix<TM, TN, Container3> result;

			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TN; j++)
					result(i, j) = m1(i, j) - m2(i, j);

			return result;
		}

		/**
		 * Computes difference of two matrices. Most useful specialization, made in assumption that result matrix should
		 * use default matrix container.
		 *
		 * @tparam TM First matrix dimension.
		 * @tparam TN Second matrix dimension.
		 * @tparam Container1 Container type for first difference item.
		 * @tparam Container2 Container type for second  difference item.
		 * @param m1 First difference item.
		 * @param m2 Second difference item.
		 *
		 * @return Matrix difference (m1-m2).
		 */
		template<int TM, int TN, typename Container1, typename Container2>
		Matrix<TM, TN, DefaultMatrixContainer<TM, TN>> operator-(const Matrix<TM, TN, Container1> &m1,
		                                                         const Matrix<TM, TN, Container2> &m2) {
			return operator-<TM, TN, Container1, Container2, DefaultMatrixContainer<TM, TN>>(m1, m2);
		}

		/**
		 * Computes difference of two matrices. Specialized implementation, made in assumption that all matrices
		 * use default matrix container.
		 *
		 * @tparam TM First matrix dimension.
		 * @tparam TN Second matrix dimension.
		 * @tparam Container Container type for all matrices.
		 * @param m1 First difference item.
		 * @param m2 Second difference item.
		 *
		 * @return Matrix difference (m1-m2).
		 */
		template<int TM, int TN, typename Container>
		Matrix<TM, TN, Container> operator-(const Matrix<TM, TN, Container> &m1, const Matrix<TM, TN, Container> &m2) {
			return operator-<TM, TN, Container, Container, Container>(m1, m2);
		}

		/**
		 * Computes product of two matrices. Generic implementation.
		 *
		 * @tparam TM First dimension of first matrix (TM x TN).
		 * @tparam TN First (second) dimension of second (first) matrix (TM x TN or TN x TK respectively).
		 * @tparam TK Second dimension of second matrix (TN x TK).
		 * @tparam Container1 Container type for first product item.
		 * @tparam Container2 Container type for second  product item.
		 * @tparam Container3 Container type for result.
		 * @param m1 First product item.
		 * @param m2 Second product item.
		 *
		 * @return Matrix product (m1*m2).
		 */
		template<int TM, int TN, int TK, typename Container1, typename Container2, typename Container3>
		Matrix<TM, TK, Container3> operator*(const Matrix<TM, TN, Container1> &m1,
		                                     const Matrix<TN, TK, Container2> &m2) {
			Matrix<TM, TK, Container3> result;

			for (int i = 0; i < TM; i++)

				for (int j = 0; j < TK; j++) {
					result(i, j) = 0;
					for (int n = 0; n < TN; n++)
						result(i, j) += m1(i, n) * m2(n, j);
				}

			return result;
		}

		/**
		 * Computes product of two matrices. Most useful specialization, made in assumption that result matrix should
		 * use default matrix container.
		 *
		 * @tparam TM First dimension of first matrix (TM x TN).
		 * @tparam TN First (second) dimension of second (first) matrix (TM x TN or TN x TK respectively).
		 * @tparam TK Second dimension of second matrix (TN x TK).
		 * @tparam Container1 Container type for first product item.
		 * @tparam Container2 Container type for second  product item.
		 * @param m1 First product item.
		 * @param m2 Second product item.
		 *
		 * @return Matrix product (m1*m2).
		 */
		template<int TM, int TN, int TK, typename Container1, typename Container2>
		Matrix<TM, TK, DefaultMatrixContainer<TM, TK>> operator*(const Matrix<TM, TN, Container1> &m1,
		                                                         const Matrix<TN, TK, Container2> &m2) {
			return operator*<TM, TN, TK, Container1, Container2, DefaultMatrixContainer<TM, TK>>(m1, m2);
		}

		/**
		 * Computes product of two matrices. Most useful specialization, made in assumption that both matrices are
		 * square and have the same container type.
		 *
		 * @tparam TM Matrix size.
		 * @tparam Container Container type for first product item.
		 * @param m1 First product item.
		 * @param m2 Second product item.
		 *
		 * @return Matrix product (m1*m2).
		 */
		template<int TM, typename Container>
		Matrix<TM, TM, Container> operator*(const Matrix<TM, TM, Container> &m1, const Matrix<TM, TM, Container> &m2) {
			return operator*<TM, TM, TM, Container, Container, Container>(m1, m2);
		}

		/**
		 * Performs scalar multiplication.
		 *
		 * @tparam TM First matrix dimesion.
		 * @tparam TN Second matrix dimension.
		 * @tparam Container Container type of matrix.
		 * @param m Matrix to multiply.
		 *
		 * @return Result of scalar multiplication.
		 */
		template<int TM, int TN, typename Container, typename TMultiplier>
		Matrix<TM, TN, Container> operator*(const Matrix<TM, TN, Container> &m, const TMultiplier x) {
			Matrix<TM, TN, Container> result;

			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TN; j++)
					result(i, j) = m(i, j) * x;

			return result;
		}

		/**
		 * Performs scalar multiplication.
		 *
		 * @tparam TM First matrix dimesion.
		 * @tparam TN Second matrix dimension.
		 * @tparam Container Container type of matrix.
		 * @param m Matrix to multiply.
		 *
		 * @return Result of scalar multiplication.
		 */
		template<int TM, int TN, typename Container, typename TMultiplier>
		Matrix<TM, TN, Container> operator*(const TMultiplier x, const Matrix<TM, TN, Container> &m) {
			return m * x;
		}

		/**
		 * Performs scalar division.
		 *
		 * @tparam TM First matrix dimension.
		 * @tparam TN Second matrix dimension.
		 * @tparam Container Container type of matrix.
		 * @param m Matrix to divide.
		 * @param x Scalar to divide by.
		 *
		 * @return Result of scalar division.
		 */
		template<int TM, int TN, typename Container>
		Matrix<TM, TN, Container> operator/(const Matrix<TM, TN, Container> &m, const real x) {
			return m * (1 / x);
		}

		/**
		 * Add m2 to m1. Generic implementation.
		 *
		 * @tparam TM First matrix dimension.
		 * @tparam TN Second matrix dimension.
		 * @tparam Container1 Container type for first summ item.
		 * @tparam Container2 Container type for second  summ item.
		 * @param m1 First summ item.
		 * @param m2 Second summ item.
		 */
		template<int TM, int TN, typename Container1, typename Container2>
		void operator+=(Matrix<TM, TN, Container1> &m1, const Matrix<TM, TN, Container2> &m2) {
			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TN; j++)
					m1(i, j) = m1(i, j) + m2(i, j);
		}

		template<int TM, int TN, typename Container>
		void operator*=(Matrix<TM, TN, Container> &m, const real x) {
			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TN; j++)
					m(i, j) *= x;
		}

		template<int TM, int TN, typename Container>
		void operator/=(Matrix<TM, TN, Container> &m, const real x) {
			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TN; j++)
					m(i, j) /= x;
		}

		template<int TM, int TN, typename Container1, typename Container2>
		bool operator==(const Matrix<TM, TN, Container1> &m1, const Matrix<TM, TN, Container2> &m2) {
			for (int i = 0; i < TM; i++)
				for (int j = 0; j < TN; j++)
					// FIXME Should this constant be replaced by something context-specific?
					if (fabs(m1(i, j) - m2(i, j)) > EQUALITY_TOLERANCE)
						return false;

			return true;
		}

		template<int TM, int TN, typename Container1, typename Container2>
		bool operator!=(const Matrix<TM, TN, Container1> &m1, const Matrix<TM, TN, Container2> &m2) {
			return !(m1 == m2);
		}

		/** @return dot product of specified vectors */
		template<int TM, typename Container1, typename Container2>
		real dotProduct(const Vector<TM, Container1> &v1, const Vector<TM, Container2> &v2) {
			return (v1.transpose() * v2)(0);
		}

		/** @return length of v */
		template<int TM, typename Container>
		real length(const Vector<TM, Container>& v) {
			return sqrt(dotProduct(v, v));
		}

		/** @return co-directional to v vector of length 1.0 */
		template<int TM, typename Container>
		Vector<TM, Container> normalize(const Vector<TM, Container>& v) {
			real l = length(v);
			assert_gt(l, 0.0);
			return v / l;
		}

		/**
		 * Computes component-by-component product of two matrices. Generic implementation.
		 *
		 * @tparam TM, TN matrix dimensions
		 * @tparam Container1 Container type for first product item.
		 * @tparam Container2 Container type for second  product item.
		 * @tparam Container3 Container type for result.
		 * @param m1 First product item.
		 * @param m2 Second product item.
		 *
		 * @return component-by-component product of m1 and m2
		 */
		template<int TM, int TN, typename Container1, typename Container2, typename Container3 = DefaultMatrixContainer<TM, TN>>
		Matrix<TM, TN, Container3> plainMultiply(const Matrix<TM, TN, Container1> &m1, const Matrix<TM, TN, Container2> &m2) {
			Matrix<TM, TN, Container3> result;

			for (int i = 0; i < TM; i++) {
				for (int j = 0; j < TN; j++) {
					result(i, j) = m1(i, j) * m2(i, j);
				}
			}

			return result;
		}

		/**
		 * Computes component-by-component division of two matrices. Generic implementation.
		 * @return component-by-component division of m1 and m2
		 * @warning if some component of m2 == 0 corresponding component of answer will be
		 * std::numeric_limits<real>::max() * sign of m1 component
		 */
		template<int TM, int TN, typename Container1, typename Container2, typename Container3 = DefaultMatrixContainer<TM, TN>>
		Matrix<TM, TN, Container3> plainDivision(const Matrix<TM, TN, Container1> &m1, const Matrix<TM, TN, Container2> &m2) {
			Matrix<TM, TN, Container3> result;

			for (int i = 0; i < TM; i++) {
				for (int j = 0; j < TN; j++) {
					result(i, j) = m1(i, j) / m2(i, j);
					if (m2(i, j) == 0) {
						assert_ne(m1(i, j), 0);
						result(i, j) = (m1(i, j) > 0 ? 1 : -1) * std::numeric_limits<real>::max();
					}
				}
			}

			return result;
		}

		/** @return multiplication of all elements */
		template<int TM, int TN, typename Container>
		real directProduct(const Matrix<TM, TN, Container>& m) {
			real ans = 1;
			for (int i = 0; i < TM; i++) {
				for (int j = 0; j < TN; j++) {
					ans *= m(i, j);
				}
			}
			return ans;
		};

	};
};

#endif // LIBGCM_LINAL_LINALROUTINES_HPP
