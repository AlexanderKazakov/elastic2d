#ifndef LIBGCM_LINAL_LINEARSYSTEMS_HPP
#define LIBGCM_LINAL_LINEARSYSTEMS_HPP

#include <lib/linal/Matrix.hpp>
#include <lib/linal/determinants.hpp>
#include <lib/linal/functions.hpp>

namespace gcm {
namespace linal {

/** 
 * Solve SLE \f$ \matrix{A} * \vec{x} = \vec{b} \f$ with TNxTN 
 * non-symmetric matrix, TN > 3.
 */
template<int TN,
         typename TVectorElement,
         template<int, typename> class TVectorContainer,
         typename TMatrixElement,
         template<int, typename> class TMatrixContainer>
MatrixBase<TN, 1,
           decltype(TVectorElement() / TMatrixElement()),
           NonSymmetric, TVectorContainer>
solveLinearSystem(const MatrixBase<TN, TN,
                                   TMatrixElement,
                                   NonSymmetric, TMatrixContainer>& A,
                  const MatrixBase<TN, 1,
                                   TVectorElement,
                                   NonSymmetric, TVectorContainer>& b);


/** 
 * Solve SLE \f$ \matrix{A} * \vec{x} = \vec{b} \f$ with 1x1 
 * non-symmetric matrix.
 */
template<typename TVectorElement,
         template<int, typename> class TVectorContainer,
         typename TMatrixElement,
         template<int, typename> class TMatrixContainer>
MatrixBase<1, 1,
           decltype(TVectorElement() / TMatrixElement()),
           NonSymmetric, TVectorContainer>
solveLinearSystem(const MatrixBase<1, 1,
                                   TMatrixElement,
                                   NonSymmetric, TMatrixContainer>& A,
                  const MatrixBase<1, 1,
                                   TVectorElement,
                                   NonSymmetric, TVectorContainer>& b) {
	
	assert_ne(A(0, 0), zeros(A(0, 0)));
	return {b(0) / A(0, 0)};
}


/** 
 * Solve SLE \f$ \matrix{A} * \vec{x} = \vec{b} \f$ with 2x2 
 * non-symmetric matrix.
 */
template<typename TVectorElement,
         template<int, typename> class TVectorContainer,
         typename TMatrixElement,
         template<int, typename> class TMatrixContainer>
MatrixBase<2, 1,
           decltype(TVectorElement() / TMatrixElement()),
           NonSymmetric, TVectorContainer>
solveLinearSystem(const MatrixBase<2, 2,
                                   TMatrixElement,
                                   NonSymmetric, TMatrixContainer>& A,
                  const MatrixBase<2, 1,
                                   TVectorElement,
                                   NonSymmetric, TVectorContainer>& b) {

	auto det = determinant(A);
	assert_ne(det, zeros(det));
	
	auto det1 = determinant(b(0), A(0, 1),
	                        b(1), A(1, 1));
	                        
	auto det2 = determinant(A(0, 0), b(0),
	                        A(1, 0), b(1));

	return {det1 / det, det2 / det};
}


/** 
 * Solve SLE \f$ \matrix{A} * \vec{x} = \vec{b} \f$ with 3x3 
 * non-symmetric matrix.
 */
template<typename TVectorElement,
         template<int, typename> class TVectorContainer,
         typename TMatrixElement,
         template<int, typename> class TMatrixContainer>
MatrixBase<3, 1,
           decltype(TVectorElement() / TMatrixElement()),
           NonSymmetric, TVectorContainer>
solveLinearSystem(const MatrixBase<3, 3,
                                   TMatrixElement,
                                   NonSymmetric, TMatrixContainer>& A,
                  const MatrixBase<3, 1,
                                   TVectorElement,
                                   NonSymmetric, TVectorContainer>& b) {

	auto det = determinant(A);
	assert_ne(det, zeros(det));
	
	auto det1 = determinant(b(0), A(0, 1), A(0, 2),
	                        b(1), A(1, 1), A(1, 2),
	                        b(2), A(2, 1), A(2, 2));
	                        
	auto det2 = determinant(A(0, 0), b(0), A(0, 2),
	                        A(1, 0), b(1), A(1, 2),
	                        A(2, 0), b(2), A(2, 2));

	auto det3 = determinant(A(0, 0), A(0, 1), b(0),
	                        A(1, 0), A(1, 1), b(1),
	                        A(2, 0), A(2, 1), b(2));

	return {det1 / det, det2 / det, det3 / det};
}


/** 
 * Solve overdetermined SLE \f$ \matrix{A} * \vec{x} = \vec{b} \f$
 * with TMxTN-matrix A, TM-vector b and TN-vector x, where TM >= TN,
 * by the Linear Least Squares Method with weight matrix W.
 * @param A matrix of SLE
 * @param b right part of SLE
 * @param W diagonal matrix of weights in linear least squares method,
 * W(i) determines the importance of the i-th row of the SLE
 * @return solution of the SLE in terms of least square deviation
 */
template<int TM, int TN,
         typename TVectorElement,
         template<int, typename> class TVectorContainer,
         typename TMatrixElement,
         template<int, typename> class TMatrixContainer>
MatrixBase<TN, 1,
           decltype(TVectorElement() / TMatrixElement()),
           NonSymmetric, TVectorContainer>
linearLeastSquares(
		const MatrixBase<TM, TN, TMatrixElement, NonSymmetric, TMatrixContainer>& A,
		const MatrixBase<TM, 1, TVectorElement, NonSymmetric, TVectorContainer>& b,
		const MatrixBase<TM, TM, TMatrixElement, Diagonal, TMatrixContainer>& W = 
				MatrixBase<TM, TM, TMatrixElement, Diagonal, TMatrixContainer>::Identity()) {
	
	static_assert(TM >= TN, "The system has to be at least determined");
	return solveLinearSystem(transposeMultiply(A, W * A), transposeMultiply(A, W * b));
}


}
}


#endif // LIBGCM_LINAL_LINEARSYSTEMS_HPP
