#ifndef LIBGCM_GSLUTILS_HPP
#define LIBGCM_GSLUTILS_HPP

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_poly.h>

#include <libgcm/util/infrastructure/infrastructure.hpp>
#include <libgcm/linal/linal.hpp>

namespace gcm {
namespace gsl_utils {

static constexpr real eps = 1e-2;


/// Convertion between linal and gsl formats
template<int TM, int TN, 
         typename TSymmetry,
         template<int, typename> class TContainer>
void setGslMatrixFromLinalMatrix(
		gsl_matrix* gslMatrix,
		const linal::MatrixBase<TM, TN, real, TSymmetry, TContainer>& linalMatrix) {
	assert_eq(TM, gslMatrix->size1);
	assert_eq(TN, gslMatrix->size2);
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			gsl_matrix_set(gslMatrix, (size_t)i, (size_t)j, linalMatrix(i, j));
		}
	}
}

template<int TM, template<int, typename> class TContainer>
void setGslVectorFromLinalVector(
		gsl_vector* gslVector,
		const linal::MatrixBase<TM, 1, real, linal::NonSymmetric, TContainer>& linalVector) {
	assert_eq(TM, gslVector->size);
	for (int i = 0; i < TM; i++) {
		gsl_vector_set(gslVector, (size_t)i, linalVector(i));
	}
}

template<int TM, int TN, 
         template<int, typename> class TContainer>
void setLinalMatrixFromGslMatrix(
		linal::MatrixBase<TM, TN, real, linal::NonSymmetric, TContainer>& linalMatrix,
		const gsl_matrix* gslMatrix) {
	assert_eq(TM, gslMatrix->size1);
	assert_eq(TN, gslMatrix->size2);
	for (int i = 0; i < TM; i++) {
		for (int j = 0; j < TN; j++) {
			linalMatrix(i, j) = gsl_matrix_get(gslMatrix, (size_t)i, (size_t)j);
		}
	}
}

template<int TM, template<int, typename> class TContainer>
void setLinalVectorFromGslVector(
		linal::MatrixBase<TM, 1, real, linal::NonSymmetric, TContainer>& linalVector,
		const gsl_vector* gslVector) {
	assert_eq(TM, gslVector->size);
	for (int i = 0; i < TM; i++) {
		linalVector(i) = gsl_vector_get(gslVector, (size_t)i);
	}
}
/// @}


/** Result of LU-decomposition in GSL format */
struct LuDecomposition {
	template<int TM, template<int, typename> class TContainer>
	LuDecomposition(
			const linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer>& m) {
		LU = gsl_matrix_alloc(TM, TM);
		setGslMatrixFromLinalMatrix(LU, m);
		permut = gsl_permutation_alloc(TM);
		gsl_set_error_handler_off();
		int gslErrorStatus = gsl_linalg_LU_decomp(LU, permut, &signum);
		if (gslErrorStatus != GSL_SUCCESS) {
			gsl_matrix_free(LU);
			gsl_permutation_free(permut);
			std::string error_message = "gsl_linalg_LU_decomp failed with " +
					std::to_string(gslErrorStatus) + "error status";
			THROW_INVALID_ARG(error_message);
		}
	}
	
	LuDecomposition() = delete;
	LuDecomposition(const LuDecomposition& other) = delete;
	~LuDecomposition() {
		gsl_matrix_free(LU);
		gsl_permutation_free(permut);
	}
	
	gsl_matrix* LU;
	gsl_permutation* permut;
	int signum;
};


template<int TM, template<int, typename> class TContainer>
linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer>
invert(
		const linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer>& m) {
	LuDecomposition lu(m);
	gsl_matrix* inv = gsl_matrix_alloc(TM, TM);
	gsl_set_error_handler_off();
	int gslErrorStatus = gsl_linalg_LU_invert(lu.LU, lu.permut, inv);
	if (gslErrorStatus != GSL_SUCCESS) {
		gsl_matrix_free(inv);
		std::string error_message = "gsl_linalg_LU_invert failed with \"" +
				std::to_string(gslErrorStatus) + "\" error status";
		THROW_INVALID_ARG(error_message);
	}
	linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer> ans;
	setLinalMatrixFromGslMatrix(ans, inv);
	gsl_matrix_free(inv);
	return ans;
}


template<int TM, template<int, typename> class TContainer>
real
determinant(
		const linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer>& m) {
	LuDecomposition lu(m);
	return gsl_linalg_LU_det(lu.LU, lu.signum);
}


template<int TM, template<int, typename> class TContainer>
linal::MatrixBase<TM, 1, real, linal::NonSymmetric, TContainer>
solveLinearSystem(
		const linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer>& A,
		const linal::MatrixBase<TM,  1, real, linal::NonSymmetric, TContainer>& b) {
	LuDecomposition lu(A);
	gsl_vector* gslB = gsl_vector_alloc(TM);
	setGslVectorFromLinalVector(gslB, b);
	gsl_vector* gslX = gsl_vector_alloc(TM);
	gsl_set_error_handler_off();
	int gslErrorStatus = gsl_linalg_LU_solve(lu.LU, lu.permut, gslB, gslX);
	if (gslErrorStatus != GSL_SUCCESS) {
		gsl_vector_free(gslX);
		gsl_vector_free(gslB);
		std::string error_message = "gsl_linalg_LU_solve failed with \"" +
				std::to_string(gslErrorStatus) + "\" error status";
		THROW_INVALID_ARG(error_message);
	}
	linal::MatrixBase<TM, 1, real, linal::NonSymmetric, TContainer> x;
	setLinalVectorFromGslVector(x, gslX);
	gsl_vector_free(gslX);
	gsl_vector_free(gslB);
	return x;
}


/**
 * Solve third order polynomial with *all real* roots
 * \f$   x^3 + p(0) * x^2 + p(1) * x + p(2)   $\f.
 * The first root is unique, the next two can be equal, or all three are equal.
 */
static Real3 solveThirdOrderPolynomial(const Real3 p) {
	double x1 = 0, x2 = 0, x3 = 0;
	int numberOfRoots = gsl_poly_solve_cubic(p(0), p(1), p(2), &x1, &x2, &x3);
	
	if (numberOfRoots != 3) {
		gsl_complex z1, z2, z3;
		gsl_poly_complex_solve_cubic(p(0), p(1), p(2), &z1, &z2, &z3);
		
		assert_lt(std::fabs(z1.dat[1]), std::fabs(z1.dat[0]) * eps);
		assert_lt(std::fabs(z2.dat[1]), std::fabs(z2.dat[0]) * eps);
		assert_lt(std::fabs(z3.dat[1]), std::fabs(z3.dat[0]) * eps);
		
		x1 = z1.dat[0]; x2 = z2.dat[0]; x3 = z3.dat[0];
	}
	
	// now, sort such that two equal roots are always at the end
	if (std::fabs(x1 - x2) < std::fmax(std::fabs(x1), std::fabs(x2)) * eps) {
	
		if (std::fabs(x3 - x2) < std::fmax(std::fabs(x3), std::fabs(x2)) * eps) {
		// all equal
			x1 = x2 = x3 = (x1 + x2 + x3) / 3;
		
		} else {
		// x1 == x2 != x3
			x2 = (x1 + x2) / 2;
			x1 = x3;
			x3 = x2;
		}
		
	} else if (std::fabs(x1 - x3) < std::fmax(std::fabs(x1), std::fabs(x3)) * eps) {
	// x1 == x3 != x2
		x3 = (x1 + x3) / 2;
		x1 = x2;
		x2 = x3;
		
	} else if (std::fabs(x2 - x3) < std::fmax(std::fabs(x2), std::fabs(x3)) * eps) {
	// x1 != x2 == x3
		x2 = x3 = (x2 + x3) / 2;
	}
	
	return {x1, x2, x3};
}


} // namespace gsl_utils
} // namespace gcm

#endif /* LIBGCM_GSLUTILS_HPP */
