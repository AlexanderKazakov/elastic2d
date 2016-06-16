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
           typename std::remove_cv<decltype(TVectorElement() / TMatrixElement())>::type,
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
           typename std::remove_cv<decltype(TVectorElement() / TMatrixElement())>::type,
           NonSymmetric, TVectorContainer>
solveLinearSystem(const MatrixBase<1, 1,
                                   TMatrixElement,
                                   NonSymmetric, TMatrixContainer>& A,
                  const MatrixBase<1, 1,
                                   TVectorElement,
                                   NonSymmetric, TVectorContainer>& b) {
	
	if (A(0, 0) == zeros(A(0, 0))) {
		THROW_INVALID_ARG("SLE determinant is zero");
	}
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
           typename std::remove_cv<decltype(TVectorElement() / TMatrixElement())>::type,
           NonSymmetric, TVectorContainer>
solveLinearSystem(const MatrixBase<2, 2,
                                   TMatrixElement,
                                   NonSymmetric, TMatrixContainer>& A,
                  const MatrixBase<2, 1,
                                   TVectorElement,
                                   NonSymmetric, TVectorContainer>& b) {

	auto det = determinant(A);
	if (det == zeros(det)) {
		THROW_INVALID_ARG("SLE determinant is zero");
	}
	
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
           typename std::remove_cv<decltype(TVectorElement() / TMatrixElement())>::type,
           NonSymmetric, TVectorContainer>
solveLinearSystem(const MatrixBase<3, 3,
                                   TMatrixElement,
                                   NonSymmetric, TMatrixContainer>& A,
                  const MatrixBase<3, 1,
                                   TVectorElement,
                                   NonSymmetric, TVectorContainer>& b) {

	auto det = determinant(A);
	if (det == zeros(det)) {
		THROW_INVALID_ARG("SLE determinant is zero");
	}
	
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
           typename std::remove_cv<decltype(TVectorElement() / TMatrixElement())>::type,
           NonSymmetric, TVectorContainer>
linearLeastSquares(
		const MatrixBase<TM, TN, TMatrixElement, NonSymmetric, TMatrixContainer>& A,
		const MatrixBase<TM, 1, TVectorElement, NonSymmetric, TVectorContainer>& b,
		const MatrixBase<TM, TM, TMatrixElement, Diagonal, TMatrixContainer>& W = 
				MatrixBase<TM, TM, TMatrixElement, Diagonal, TMatrixContainer>::Identity()) {
	
	static_assert(TM >= TN, "The system has to be at least determined");
	return solveLinearSystem(transposeMultiply(A, W * A), transposeMultiply(A, W * b));
}


/** 
 * Solve SLE \f$   \matrix{A} * \vec{x} = \vec{0}   \f$ 
 * with 3x3 non-symmetric degenerate matrix with rank 1 or 2.
 * I.e, if rank == 2, numberOfSolutions == 1,
 *      if rank == 1, numberOfSolutions == 2.
 * Size of returned std::vector is equal to numberOfSolutions.
 * The method is stable to numerical inexactness.
 */
template<typename TMatrixElement,
         template<int, typename> class TMatrixContainer>
std::vector<MatrixBase<3, 1,
            TMatrixElement,
            NonSymmetric, TMatrixContainer>>
solveDegenerateLinearSystem(const MatrixBase<3, 3,
                            TMatrixElement,
                            NonSymmetric, TMatrixContainer>& A,
                            const int numberOfSolutions) {
	
	typedef MatrixBase<3, 1,
            TMatrixElement,
            NonSymmetric, TMatrixContainer> Answer;

	switch (numberOfSolutions) {
	case 1: // rank(A) == 2
	{
		int I = 0, J = 1, P = 0, Q = 1;
		// find the most determined matrix minor
		TMatrixElement det = 0;
		for(int i = 0; i < 2; i++)
		for(int j = i+1; j < 3; j++) {
			for(int p = 0; p < 2; p++)
			for(int q = p+1; q < 3; q++) {
			
				if (std::fabs(A(p, i)*A(q, j) - A(q, i)*A(p, j)) > std::fabs(det)) {
					det = A(p, i)*A(q, j) - A(q, i)*A(p, j);
					I = i; J = j; P = p; Q = q;
				}
				
			}
		}
		
		int U = Utils::other012(I, J); //< no I and no J
		
		Answer x;
		x(U) = 1;
		x(I) = (-A(P, U)*A(Q, J) + A(Q, U)*A(P, J)) / det;
		x(J) = (-A(P, I)*A(Q, U) + A(Q, I)*A(P, U)) / det;
		
		return {x};
	}
	
	case 2: // rank(A) == 1
	{
		int I = 0, J = 0;
		// find the most determined matrix minor
		TMatrixElement det = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (std::fabs(det) < std::fabs(A(i, j))) {
					det = A(i, j);
					I = i; J = j;
				}
			}
		}
		
		int p, q;
		// p and q = no J
		if      (J == 0) { p = 1; q = 2; }
		else if (J == 1) { p = 0; q = 2; }
		else             { p = 0; q = 1; }
		
		Answer x, y;
		x(p) = y(q) = 1;
		x(q) = y(p) = 0;
		x(J) = -A(I, p) / det;
		y(J) = -A(I, q) / det;
		
		return {x, y};
	}
	
	default:
		THROW_INVALID_ARG("Invalid number of linearly independent solutions required");
	}
}


}
}


#endif // LIBGCM_LINAL_LINEARSYSTEMS_HPP
