#ifndef LIBGCM_LINAL_FUNCTIONS_HPP
#define LIBGCM_LINAL_FUNCTIONS_HPP

#include <limits>

#include <lib/util/Utils.hpp>
#include <lib/linal/Matrix.hpp>
#include <lib/linal/determinants.hpp>

namespace gcm {
namespace linal {


/**
 * Transpose usual matrix
 */
template<int TM, int TN, 
         typename TElement,
         template<int, typename> class TContainer>
MatrixBase<TN, TM, TElement, NonSymmetric, TContainer>
transpose(const MatrixBase<TM, TN, TElement, NonSymmetric, TContainer>& m) {
	MatrixBase<TN, TM, TElement, NonSymmetric, TContainer> result;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(j, i) = m(i, j);
		}
	}
	return result;
}


/**
 * Transpose symmetric matrix
 */
template<int TM, int TN, 
         typename TElement,
         template<int, typename> class TContainer>
MatrixBase<TN, TM, TElement, Symmetric, TContainer>
transpose(const MatrixBase<TM, TN, TElement, Symmetric, TContainer>& m) {
	return m;
}


/**
 * Transpose diagonal matrix
 */
template<int TM, int TN, 
         typename TElement,
         template<int, typename> class TContainer>
MatrixBase<TN, TM, TElement, Diagonal, TContainer>
transpose(const MatrixBase<TM, TN, TElement, Diagonal, TContainer>& m) {
	return m;
}


/**
 * Invert 1x1 matrix
 */
template<typename TElement,
         template<int, typename> class TContainer>
MatrixBase<1, 1, TElement, NonSymmetric, TContainer>
invert(const MatrixBase<1, 1, TElement, NonSymmetric, TContainer>& m) {
	return MatrixBase<1, 1, TElement, NonSymmetric, TContainer>({1.0 / m(0, 0)});
}


/**
 * Invert 2x2 matrix
 */
template<typename TElement,
         template<int, typename> class TContainer>
MatrixBase<2, 2, TElement, NonSymmetric, TContainer>
invert(const MatrixBase<2, 2, TElement, NonSymmetric, TContainer>& m) {
	return MatrixBase<2, 2, TElement, NonSymmetric, TContainer>(
			{ m(1, 1), -m(0, 1),
			 -m(1, 0),  m(0, 0) }) / determinant(m);
}


/**
 * Invert 3x3 matrix
 */
template<typename TElement,
         template<int, typename> class TContainer>
MatrixBase<3, 3, TElement, NonSymmetric, TContainer>
invert(const MatrixBase<3, 3, TElement, NonSymmetric, TContainer>& m) {
	return MatrixBase<3, 3, TElement, NonSymmetric, TContainer>(
{ m(1,1)*m(2,2)-m(1,2)*m(2,1), m(0,2)*m(2,1)-m(0,1)*m(2,2), m(0,1)*m(1,2)-m(1,1)*m(0,2),
  m(1,2)*m(2,0)-m(1,0)*m(2,2), m(0,0)*m(2,2)-m(0,2)*m(2,0), m(0,2)*m(1,0)-m(0,0)*m(1,2),
  m(1,0)*m(2,1)-m(1,1)*m(2,0), m(0,1)*m(2,0)-m(0,0)*m(2,1), m(0,0)*m(1,1)-m(0,1)*m(1,0) }
			) / determinant(m);
}


/**
 * Invert diagonal matrix
 */
template<int TM, int TN,
         typename TElement,
         template<int, typename> class TContainer>
MatrixBase<TN, TM, TElement, Diagonal, TContainer>
invert(const MatrixBase<TM, TN, TElement, Diagonal, TContainer>& m) {
	MatrixBase<TN, TM, TElement, Diagonal, TContainer> result;
	for (int i = 0; i < TM; i++) {
		assert_ne(m(i, i), 0);
		result(i, i) = 1.0 / m(i, i);
	}
	return result;
}


/**
 * Norm of the matrix consistent with maximum (aka infinity) norm of vectors
 * @see Petrov, Lobanov "Lections on numerical mathematics" page 34
 */
template<int TM,
         typename TElement,
         typename TSymmetry,
         template<int, typename> class TContainer>
TElement normMax(const MatrixBase<TM, TM, TElement, TSymmetry, TContainer>& m) {
	TElement ans = zeros(TElement());
	for (int i = 0; i < TM; i++) {
		TElement ith = zeros(TElement());
		for (int j = 0; j < TM; j++) {
			ith += fabs(m(i, j));
		}
		if (ans < ith) { ans = ith; }
	}
	return ans;
}


/**
 * Condition number of matrix A in maximum (aka infinity) norm
 * \f$  mu = normMax(A) * normMax(invert(A))  \f$
 */
template<int TM,
         typename TElement,
         typename TSymmetry,
         template<int, typename> class TContainer>
TElement conditionNumber(const MatrixBase<TM, TM, TElement, TSymmetry, TContainer>& A) {
	return normMax(A) * normMax(invert(A));
}


/**
 * Computes product of two matrices, but the first one is transposed:
 * C = transpose(m1) * m2.
 * m1: TNxTM.
 * m2: TNxTK.
 * C:  TMxTK.
 * @return equal to ( transpose(m1) * m2 )
 */
template<int TM, int TN, int TK,
         typename TElement1,
         typename TSymmetry1,
         template<int, typename> class TContainer1,
         typename TElement2,
         typename TSymmetry2,
         template<int, typename> class TContainer2,
         template<int, typename> class TContainer3 = DefaultContainer>
MatrixBase<TM, TK,
           typename std::remove_cv<decltype(TElement1() * TElement2())>::type,
           NonSymmetric, TContainer3>
transposeMultiply(const MatrixBase<TN, TM, TElement1, TSymmetry1, TContainer1>& m1,
                  const MatrixBase<TN, TK, TElement2, TSymmetry2, TContainer2>& m2) {
	MatrixBase<TM, TK,
	           typename std::remove_cv<decltype(TElement1() * TElement2())>::type,
	           NonSymmetric, TContainer3> result;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TK; j++) {
			result(i, j) = m1(0, i) * m2(0, j);
			for (int n = 1; n < TN; n++) {
				result(i, j) += m1(n, i) * m2(n, j);
			}
		}
	}
	return result;
}


/**
 * Let matrix \f$ \matrix{C} = \matrix{A} * \matrix{B} \f$.
 * Then diagonal of C is returned. 
 * Calculation of non-diagonal elements of C is not performed.
 * @return diagonal of A * B
 */
template<int TM,
         typename TElement1,
         typename TSymmetry1,
         template<int, typename> class TContainer1,
         typename TElement2,
         typename TSymmetry2,
         template<int, typename> class TContainer2,
         template<int, typename> class TContainer3 = DefaultContainer>
MatrixBase<TM, 1,
           typename std::remove_cv<decltype(TElement1() * TElement2())>::type,
           NonSymmetric, TContainer3>
diagonalMultiply(const MatrixBase<TM, TM, TElement1, TSymmetry1, TContainer1>& A,
                 const MatrixBase<TM, TM, TElement2, TSymmetry2, TContainer2>& B) {

	MatrixBase<TM, 1,
               typename std::remove_cv<decltype(TElement1() * TElement2())>::type,
               NonSymmetric, TContainer3> result;
	for (int i = 0; i < TM; i++) {
		result(i) = A(i, 0) * B(0, i);
		for (int j = 1; j < TM; j++) {
			result(i) += A(i, j) * B(j, i);
		}
	}
	return result;
}


/** 
 * @return summ of diagonal elements of the matrix A
 */
template<typename TMatrix>
typename TMatrix::ElementType
trace(const TMatrix& A) {
	static_assert(TMatrix::M == TMatrix::N, "Square matrices only have trace");
	auto result = A(0, 0);
	for (int i = 1; i < TMatrix::M; i++) {
		result += A(i, i);
	}
	return result;
}


/** 
 * @return matrix A diagonal into vector 
 */
template<typename TMatrix>
typename TMatrix::ColumnType
diag(const TMatrix& m) {
	typename TMatrix::ColumnType result;
	for (int i = 0; i < TMatrix::M; i++) {
		result(i) = m(i, i);
	}
	return result;
}


/** 
 * @return diagonal of the matrix into diagonal matrix
 */
template<int TM,
         typename TElement,
         typename TSymmetry1,
         template<int, typename> class TContainer1,
         template<int, typename> class TContainer2 = DefaultContainer>
MatrixBase<TM, TM, TElement, Diagonal, TContainer2>
Diag(const MatrixBase<TM, TM, TElement, TSymmetry1, TContainer1>& m) {
	MatrixBase<TM, TM, TElement, Diagonal, TContainer2> ans;
	for (int i = 0; i < TM; i++) {
		ans(i) = m(i, i);
	}
	return ans;
}


/** 
 * @return dot (aka scalar) product of specified vectors v1, v2
 * @note equal to transpose(v1) * v2
 */
template<int TM,
         typename TElement1,
         template<int, typename> class TContainer1,
         typename TElement2,
         template<int, typename> class TContainer2>
typename std::remove_cv<decltype(TElement1() * TElement2())>::type
dotProduct(const MatrixBase<TM, 1, TElement1, NonSymmetric, TContainer1>& v1,
           const MatrixBase<TM, 1, TElement2, NonSymmetric, TContainer2>& v2) {
	auto result = v1(0) * v2(0);
	for (int i = 1; i < TM; i++) {
		result += v1(i) * v2(i);
	}
	return result;
}


/** 
 * @return dot (aka scalar) product of specified vectors v1, v2
 * with specified symmetric Gramian matrix H
 * @note equal to transpose(v1) * H * v2
 */
template<int TM,
         typename TElement1,
         template<int, typename> class TContainer1,
         typename TElementH,
         template<int, typename> class TContainerH,
         typename TElement2,
         template<int, typename> class TContainer2>
typename std::remove_cv<decltype(TElement1() * TElementH() * TElement2())>::type
dotProduct(const MatrixBase<TM,  1, TElement1, NonSymmetric, TContainer1>& v1,
           const MatrixBase<TM, TM, TElementH,    Symmetric, TContainerH>& H,
           const MatrixBase<TM,  1, TElement2, NonSymmetric, TContainer2>& v2) {
	auto v1H = transposeMultiply(v1, H);
	auto result = v1H(0) * v2(0);
	for (int i = 1; i < TM; i++) {
		result += v1H(i) * v2(i);
	}
	return result;
}


/** 
 * @return length of vector v in Euclidian space 
 */
template<typename TVector>
real length(const TVector& v) {
	return sqrt(dotProduct(v, v));
}


/** 
 * @return co-directional to v vector of unit length
 */
template<typename TVector>
TVector normalize(const TVector& v) {
	real l = length(v);
	assert_gt(l, 0.0);
	return v / l;
}


/**
 * Element-by-element multiplication
 */
template<int TM, int TN,
         typename TElement1,
         typename TSymmetry1,
         template<int, typename> class TContainer1,
         typename TElement2,
         typename TSymmetry2,
         template<int, typename> class TContainer2,
         template<int, typename> class TContainer3 = DefaultContainer>
MatrixBase<TM, TN,
           typename std::remove_cv<decltype(TElement1() * TElement2())>::type,
           NonSymmetric, TContainer3>
plainMultiply(const MatrixBase<TM, TN, TElement1, TSymmetry1, TContainer1>& m1,
              const MatrixBase<TM, TN, TElement2, TSymmetry2, TContainer2>& m2) {
	
	MatrixBase<TM, TN,
	           typename std::remove_cv<decltype(TElement1() * TElement2())>::type,
	           NonSymmetric, TContainer3> result;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(i, j) = m1(i, j) * m2(i, j);
		}
	}
	return result;
}


/**
 * Element-by-element division m1 by m2.
 * @note if some component of m2 == 0 corresponding component of answer will be
 * Utils::sign(m1(i, j)) * std::numeric_limits<real>::max()
 */
template<int TM, int TN,
         typename TElement1,
         typename TSymmetry1,
         template<int, typename> class TContainer1,
         typename TElement2,
         typename TSymmetry2,
         template<int, typename> class TContainer2,
         template<int, typename> class TContainer3 = DefaultContainer>
MatrixBase<TM, TN,
           typename std::remove_cv<decltype(TElement1() / TElement2())>::type,
           NonSymmetric, TContainer3>
plainDivision(const MatrixBase<TM, TN, TElement1, TSymmetry1, TContainer1>& m1,
              const MatrixBase<TM, TN, TElement2, TSymmetry2, TContainer2>& m2) {
	
	MatrixBase<TM, TN,
	           typename std::remove_cv<decltype(TElement1() / TElement2())>::type,
	           NonSymmetric, TContainer3> result;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(i, j) = m1(i, j) / m2(i, j);
			if (m2(i, j) == 0) {
				result(i, j) = Utils::sign(m1(i, j)) * std::numeric_limits<real>::max();
			}
		}
	}
	return result;
}


/** 
 * @return multiplication of all matrix elements 
 */
template<typename TMatrix>
typename TMatrix::ElementType
directProduct(const TMatrix& m) {
	typename TMatrix::ElementType result = identity(m(0, 0));
	for (int i = 0; i < TMatrix::M; i++) {
		for (int j = 0; j < TMatrix::N; j++) {
			result *= m(i, j);
		}
	}
	return result;
}


/**
 * Test on APPROXIMATE equality.
 */
template<int TM, int TN,
         typename TElement1,
         typename TSymmetry1,
         template<int, typename> class TContainer1,
         typename TElement2,
         typename TSymmetry2,
         template<int, typename> class TContainer2>
bool approximatelyEqual(
		const MatrixBase<TM, TN, TElement1, TSymmetry1, TContainer1>& m1,
        const MatrixBase<TM, TN, TElement2, TSymmetry2, TContainer2>& m2,
        const real tolerance = EQUALITY_TOLERANCE) {
    // TODO - or check (m1 - m2) with some norm?
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			if ( !Utils::approximatelyEqual(m1(i, j), m2(i, j), tolerance) ) {
				return false;
			}
		}
	}
	return true;
}


/** 
 * @return random matrix (all values are uniformly distributed from min to max)
 */
template<typename TMatrix>
TMatrix
random(const real min = 0, const real max = 1) {
	TMatrix result;
	for (int i = 0; i < TMatrix::M; i++) {
		for (int j = 0; j < TMatrix::N; j++) {
			result(i, j) = Utils::randomReal(min, max);
		}
	}
	return result;
}


/**
 * Direct product of two vectors:
 * p(i, j) = v1(i) * v2(j)
 */
template<int TM, int TN,
         typename TElement1,
         template<int, typename> class TContainer1,
         typename TElement2,
         template<int, typename> class TContainer2,
         template<int, typename> class TContainer3 = DefaultContainer>
MatrixBase<TM, TN,
           typename std::remove_cv<decltype(TElement1() * TElement2())>::type,
           NonSymmetric, TContainer3>
directProduct(const MatrixBase<TM, 1, TElement1, NonSymmetric, TContainer1>& v1,
              const MatrixBase<TN, 1, TElement2, NonSymmetric, TContainer2>& v2) {
	
	MatrixBase<TM, TN,
	           typename std::remove_cv<decltype(TElement1() * TElement2())>::type,
	           NonSymmetric, TContainer3> result;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			result(i, j) = v1(i) * v2(j);
		}
	}
	return result;
}


/**
 * Symmetrized direct product of two vectors:
 * p(i, j) = ( v1(i) * v2(j) + v2(i) * v1(j) ) / 2
 */
template<int TM,
         typename TElement1,
         template<int, typename> class TContainer1,
         typename TElement2,
         template<int, typename> class TContainer2,
         template<int, typename> class TContainer3 = DefaultContainer>
MatrixBase<TM, TM,
           typename std::remove_cv<decltype(TElement1() * TElement2())>::type,
           Symmetric, TContainer3>
symmDirectProduct(const MatrixBase<TM, 1, TElement1, NonSymmetric, TContainer1>& v1,
                  const MatrixBase<TM, 1, TElement2, NonSymmetric, TContainer2>& v2) {
	
	MatrixBase<TM, TM,
	           typename std::remove_cv<decltype(TElement1() * TElement2())>::type,
	           Symmetric, TContainer3> result;
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j <= i; j++) {
			result(i, j) = ( v1(i) * v2(j) + v2(i) * v1(j) ) / 2;
		}
	}
	return result;
}


}
}

#endif // LIBGCM_LINAL_FUNCTIONS_HPP
