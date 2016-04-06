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
TContainer& clear(TContainer& container) {
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
Matrix<TM, TN, Container> operator-(const Matrix<TM, TN, Container>& m) {
	Matrix<TM, TN, Container> result;

	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(i, j) = -m(i, j);
		}
	}

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
Matrix<TM, TN, Container3> operator+(const Matrix<TM, TN, Container1>& m1,
                                     const Matrix<TM, TN, Container2>& m2) {
	Matrix<TM, TN, Container3> result;

	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(i, j) = m1(i, j) + m2(i, j);
		}
	}

	return result;
}


/**
 * Computes summ of two matrices. Most useful specialization, made in assumption that result matrix
 *should
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
Matrix<TM, TN, DefaultMatrixContainer<TM, TN> > operator+(const Matrix<TM, TN, Container1>& m1,
                                                          const Matrix<TM, TN, Container2>& m2) {
	return operator+<TM, TN, Container1, Container2, DefaultMatrixContainer<TM, TN> >(m1, m2);
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
Matrix<TM, TN, Container> operator+(const Matrix<TM, TN, Container>& m1, const Matrix<TM, TN,
                                                                                      Container>&
                                    m2) {
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
Matrix<TM, TN, Container3> operator-(const Matrix<TM, TN, Container1>& m1,
                                     const Matrix<TM, TN, Container2>& m2) {
	Matrix<TM, TN, Container3> result;

	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(i, j) = m1(i, j) - m2(i, j);
		}
	}

	return result;
}


/**
 * Computes difference of two matrices. Most useful specialization, made in assumption that result
 *matrix should
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
Matrix<TM, TN, DefaultMatrixContainer<TM, TN> > operator-(const Matrix<TM, TN, Container1>& m1,
                                                          const Matrix<TM, TN, Container2>& m2) {
	return operator-<TM, TN, Container1, Container2, DefaultMatrixContainer<TM, TN> >(m1, m2);
}


/**
 * Computes difference of two matrices. Specialized implementation, made in assumption that all
 *matrices
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
Matrix<TM, TN, Container> operator-(const Matrix<TM, TN, Container>& m1, const Matrix<TM, TN,
                                                                                      Container>&
                                    m2) {
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
Matrix<TM, TK, Container3> operator*(const Matrix<TM, TN, Container1>& m1,
                                     const Matrix<TN, TK, Container2>& m2) {
	Matrix<TM, TK, Container3> result;

	for (int i = 0; i < TM; i++) {

		for (int j = 0; j < TK; j++) {
			result(i, j) = 0;
			for (int n = 0; n < TN; n++) {
				result(i, j) += m1(i, n) * m2(n, j);
			}
		}
	}

	return result;
}


/**
 * Computes product of two matrices. Most useful specialization, made in assumption that result
 *matrix should
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
Matrix<TM, TK, DefaultMatrixContainer<TM, TK> > operator*(const Matrix<TM, TN, Container1>& m1,
                                                          const Matrix<TN, TK, Container2>& m2) {
	return operator*<TM, TN, TK, Container1, Container2,
	                 DefaultMatrixContainer<TM, TK> >(m1, m2);
}


/**
 * Computes product of two matrices. Most useful specialization, made in assumption that both
 *matrices are
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
Matrix<TM, TM, Container> operator*(const Matrix<TM, TM, Container>& m1, const Matrix<TM, TM,
                                                                                      Container>&
                                    m2) {
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
Matrix<TM, TN, Container> operator*(const Matrix<TM, TN, Container>& m, const TMultiplier x) {
	Matrix<TM, TN, Container> result;

	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(i, j) = m(i, j) * x;
		}
	}

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
Matrix<TM, TN, Container> operator*(const TMultiplier x, const Matrix<TM, TN, Container>& m) {
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
Matrix<TM, TN, Container> operator/(const Matrix<TM, TN, Container>& m, const real x) {
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
void operator+=(Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN, Container2>& m2) {
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			m1(i, j) = m1(i, j) + m2(i, j);
		}
	}
}


template<int TM, int TN, typename Container>
void operator*=(Matrix<TM, TN, Container>& m, const real x) {
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			m(i, j) *= x;
		}
	}
}


template<int TM, int TN, typename Container>
void operator/=(Matrix<TM, TN, Container>& m, const real x) {
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			m(i, j) /= x;
		}
	}
}


template<int TM, int TN, typename Container1, typename Container2>
bool operator==(const Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN, Container2>& m2) {
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			// FIXME Should this constant be replaced by something context-specific?
			if (fabs(m1(i, j) - m2(i, j)) > EQUALITY_TOLERANCE) {
				return false;
			}
		}
	}

	return true;
}


template<int TM, int TN, typename Container1, typename Container2>
bool operator!=(const Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN, Container2>& m2) {
	return !(m1 == m2);
}


/** @return dot product of specified vectors */
template<int TM, typename Container1, typename Container2>
real dotProduct(const Vector<TM, Container1>& v1, const Vector<TM, Container2>& v2) {
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
template<int TM, int TN, typename Container1, typename Container2, typename Container3 =
                 DefaultMatrixContainer<TM, TN> >
Matrix<TM, TN, Container3> plainMultiply(const Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN,
                                                                                            Container2>
                                         & m2) {
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
template<int TM, int TN, typename Container1, typename Container2, typename Container3 =
                 DefaultMatrixContainer<TM, TN> >
Matrix<TM, TN, Container3> plainDivision(const Matrix<TM, TN, Container1>& m1, const Matrix<TM, TN,
                                                                                            Container2>
                                         & m2) {
	Matrix<TM, TN, Container3> result;

	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(i, j) = m1(i, j) / m2(i, j);
			if (m2(i, j) == 0) {
				assert_ne(m1(i, j), 0);
				result(i,
				       j) =
				        (m1(i, j) > 0 ? 1 : -1) * std::numeric_limits<real>::max();
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
}


/** Solve SLE \f$ A \vec{x} = \vec{b} \f$ with NxN matrix, N > 3 */
template<int N, typename MatrixContainer, typename VectorContainer>
Vector<N> solveLinearSystem(const Matrix<N, N, MatrixContainer>& A,
                            const Vector<N, VectorContainer>& b);

/** Solve SLE \f$ A \vec{x} = \vec{b} \f$ with 1x1 matrix */
template<typename MatrixContainer, typename VectorContainer>
Vector<1> solveLinearSystem(const Matrix<1, 1, MatrixContainer>& A,
                            const Vector<1, VectorContainer>& b) {
	return {b(0) / A(0, 0)};
}

/** Solve SLE \f$ A \vec{x} = \vec{b} \f$ with 2x2 matrix */
template<typename MatrixContainer, typename VectorContainer>
Vector<2> solveLinearSystem(const Matrix<2, 2, MatrixContainer>& A,
                            const Vector<2, VectorContainer>& b) {
	real det = determinant(A);
	assert_gt(fabs(det), 0);
	
	real det1 = determinant(b(0), A(0, 1),
	                        b(1), A(1, 1));
	                        
	real det2 = determinant(A(0, 0), b(0),
	                        A(1, 0), b(1));

	return {det1 / det, det2 / det};
}

/** Solve SLE \f$ A \vec{x} = \vec{b} \f$ with 3x3 matrix */
template<typename MatrixContainer, typename VectorContainer>
Vector<3> solveLinearSystem(const Matrix<3, 3, MatrixContainer>& A,
                            const Vector<3, VectorContainer>& b) {
	real det = determinant(A);
	assert_gt(fabs(det), 0);
	
	real det1 = determinant(b(0), A(0, 1), A(0, 2),
	                        b(1), A(1, 1), A(1, 2),
	                        b(2), A(2, 1), A(2, 2));
	                        
	real det2 = determinant(A(0, 0), b(0), A(0, 2),
	                        A(1, 0), b(1), A(1, 2),
	                        A(2, 0), b(2), A(2, 2));

	real det3 = determinant(A(0, 0), A(0, 1), b(0),
	                        A(1, 0), A(1, 1), b(1),
	                        A(2, 0), A(2, 1), b(2));

	return {det1 / det, det2 / det, det3 / det};
}

/** 
 * @name Local basis creation
 * Create local orthogonal basis for given vector n of unit length.
 * @warning vector n MUST BE of unit length
 * 
 * In created basis, given vector n is always at last position
 * and vector at the first position tau_1 is always in XY-plane.
 * {tau_1, tau_2, n} is right-hand orthogonal triple of vectors.
 * 
 * In 3D, it's impossible to go around the whole sphere with local basis,
 * continuos at every point. In current realization, the local basis is 
 * continuos while going around the sphere slice by XY-plane, but it has 
 * discontinuities while going through the points with normal(0, 0, n).
 * 
 *  
 *          Z |      tau_2
 *            |    / 
 *  n         |   / 
 *    \__     |  /
 *       \__  | /
 *          \ |/            Y
 *            .-\-----------
 *           / \/
 *          /   \
 *         /     \
 *        /       \ 
 *     X /          tau_1
 * 
 * 
 * @return orthogonal transfer matrix from created basis to global {X, Y, Z} basis,
 * i.e matrix with vectors of created basis in columns.
 */
 ///@{
inline Matrix11 createLocalBasis(const Real1& n) {
	return Matrix11({n(0)});
}

inline Matrix22 createLocalBasis(const Real2& n) {
	const Real2 tau = perpendicularClockwise(n);
	return Matrix22({tau(0), n(0),
	                 tau(1), n(1)});
}

inline Matrix33 createLocalBasis(const Real3& n) {
	const Real3 tau_1 = normalize(perpendicularClockwise(n));
	const Real3 tau_2 = linal::crossProduct(n, tau_1);
	return Matrix33({tau_1(0), tau_2(0), n(0),
	                 tau_1(1), tau_2(1), n(1),
	                 tau_1(2), tau_2(2), n(2)});
}
 ///@}

}
}

#endif // LIBGCM_LINAL_LINALROUTINES_HPP
