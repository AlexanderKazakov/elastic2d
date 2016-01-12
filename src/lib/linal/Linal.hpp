#ifndef LIBGCM_LINAL_HPP
#define LIBGCM_LINAL_HPP

#include "lib/linal/Matrix.hpp"

namespace gcm {
	namespace linal {

		/**
		 * Generic vector - just a matrix with one column
		 */
		template<int TM, typename Container = DefaultMatrixContainer<TM, 1>>
		using Vector = Matrix<TM, 1, Container>;

		/** Set all the container's values to zero */
		template<class TContainer>
		void clear(TContainer &container) {
			memset(container.values, 0, TContainer::SIZE * sizeof(real));
		};

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
		 * @param x Scalar to multiply by.
		 *
		 * @return Result of scalar multiplication.
		 */
		template<int TM, int TN, typename Container>
		Matrix<TM, TN, Container> operator*(const Matrix<TM, TN, Container> &m, const real x) {
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
		 * @param x Scalar to multiply by.
		 * @param m Matrix to multiply.
		 *
		 * @return Result of scalar multiplication.
		 */
		template<int TM, int TN, typename Container>
		Matrix<TM, TN, Container> operator*(const real x, const Matrix<TM, TN, Container> &m) {
			return m * x;
		}

		/**
		 * Performs scalar division.
		 *
		 * @tparam TM First matrix dimesion.
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
	};
};

#endif // LIBGCM_LINAL_HPP
