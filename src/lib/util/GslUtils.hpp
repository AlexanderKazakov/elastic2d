#ifndef LIBGCM_GSLUTILS_HPP
#define LIBGCM_GSLUTILS_HPP

#include <gsl/gsl_linalg.h>

#include <lib/util/Exception.hpp>

namespace gcm {
template<int M> class SymmetricMatrix;

class GslUtils {
public:
	template<typename TMatrix>
	static TMatrix invert(const TMatrix& matrix) {
		static_assert(TMatrix::M == TMatrix::N, "Only quadratic matrices can be inverted");
		TMatrix result;

		gsl_matrix* Z1 = gsl_matrix_alloc(TMatrix::M, TMatrix::M);
		gsl_matrix* Z = gsl_matrix_alloc(TMatrix::M, TMatrix::M);
		gsl_permutation* perm = gsl_permutation_alloc(TMatrix::M);
		int k;

		for (int i = 0; i < TMatrix::M; i++) {
			for (int j = 0; j < TMatrix::M; j++) {
				gsl_matrix_set(Z1, (size_t)i, (size_t)j, matrix(i, j));
			}
		}

		if (gsl_linalg_LU_decomp(Z1, perm, &k)) {
			THROW_INVALID_ARG("gsl_linalg_LU_decomp failed");
		}

		if (gsl_linalg_LU_invert(Z1, perm, Z)) {
			THROW_INVALID_ARG("gsl_linalg_LU_invert failed");
		}

		for (int i = 0; i < TMatrix::M; i++) {
			for (int j = 0; j < TMatrix::M; j++) {
				result(i, j) = gsl_matrix_get(Z, (size_t)i, (size_t)j);
			}
		}

		gsl_permutation_free(perm);
		gsl_matrix_free(Z);
		gsl_matrix_free(Z1);

		return result;
	}

	// TODO - special gsl function
	template<int M>
	static SymmetricMatrix<M> invertSymmetric(const SymmetricMatrix<M>& matrix);

};


}

#endif /* LIBGCM_GSLUTILS_HPP */
