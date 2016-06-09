#ifndef LIBGCM_GSLUTILS_HPP
#define LIBGCM_GSLUTILS_HPP

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_poly.h>

#include <lib/util/Logging.hpp>
#include <lib/linal/linal.hpp>

namespace gcm {

class GslUtils {
public:
	static constexpr real eps = 1e-2;

	template<int TM, 
	         template<int, typename> class TContainer>
	static
	linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer>
	invert(const linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer>& m) {
	
		linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer> result;
		
		gsl_set_error_handler_off();
		
		gsl_matrix* Z1 = gsl_matrix_alloc(TM, TM);
		gsl_matrix* Z = gsl_matrix_alloc(TM, TM);
		gsl_permutation* perm = gsl_permutation_alloc(TM);
		int k;

		for (int i = 0; i < TM; i++) {
			for (int j = 0; j < TM; j++) {
				gsl_matrix_set(Z1, (size_t)i, (size_t)j, m(i, j));
			}
		}

		int gsl_err_status = gsl_linalg_LU_decomp(Z1, perm, &k);
		if (gsl_err_status != GSL_SUCCESS) {
			std::string error_message = "gsl_linalg_LU_decomp failed with " +
					std::to_string(gsl_err_status) + "error status";
			THROW_INVALID_ARG(error_message);
		}
		
		gsl_err_status = gsl_linalg_LU_invert(Z1, perm, Z);
		if (gsl_err_status != GSL_SUCCESS) {
			std::string error_message = "gsl_linalg_LU_invert failed with \"" +
					std::to_string(gsl_err_status) + "\" error status";
			USE_AND_INIT_LOGGER("GslUtils")
			LOG_DEBUG(error_message << "\ngiven matrix:" << m);
			THROW_INVALID_ARG(error_message);
		}

		for (int i = 0; i < TM; i++) {
			for (int j = 0; j < TM; j++) {
				result(i, j) = gsl_matrix_get(Z, (size_t)i, (size_t)j);
			}
		}

		gsl_permutation_free(perm);
		gsl_matrix_free(Z);
		gsl_matrix_free(Z1);

		return result;
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
	

	// TODO - special gsl function
	template<int TM, 
	         template<int, typename> class TContainer>
	static
	linal::MatrixBase<TM, TM, real, linal::Symmetric, TContainer>
	invertSymmetric(
		const linal::MatrixBase<TM, TM, real, linal::Symmetric, TContainer>& m);

};


}

#endif /* LIBGCM_GSLUTILS_HPP */
