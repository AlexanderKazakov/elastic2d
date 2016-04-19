#ifndef LIBGCM_GSLUTILS_HPP
#define LIBGCM_GSLUTILS_HPP

#include <gsl/gsl_linalg.h>

#include <lib/linal/Matrix.hpp>

namespace gcm {

class GslUtils {
public:
	template<int TM, 
	         template<int, typename> class TContainer>
	static
	linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer>
	invert(const linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer>& m) {
	
		linal::MatrixBase<TM, TM, real, linal::NonSymmetric, TContainer> result;

		gsl_matrix* Z1 = gsl_matrix_alloc(TM, TM);
		gsl_matrix* Z = gsl_matrix_alloc(TM, TM);
		gsl_permutation* perm = gsl_permutation_alloc(TM);
		int k;

		for (int i = 0; i < TM; i++) {
			for (int j = 0; j < TM; j++) {
				gsl_matrix_set(Z1, (size_t)i, (size_t)j, m(i, j));
			}
		}

		if (gsl_linalg_LU_decomp(Z1, perm, &k)) {
			THROW_INVALID_ARG("gsl_linalg_LU_decomp failed");
		}

		if (gsl_linalg_LU_invert(Z1, perm, Z)) {
			THROW_INVALID_ARG("gsl_linalg_LU_invert failed");
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
