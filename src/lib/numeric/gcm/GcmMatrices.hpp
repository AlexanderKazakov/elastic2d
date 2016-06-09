#ifndef LIBGCM_GCMMATRICES_HPP
#define LIBGCM_GCMMATRICES_HPP

#include <array>

#include <lib/linal/linal.hpp>
#include <lib/util/Types.hpp>
#include <lib/util/Enum.hpp>
#include <lib/config.hpp>

namespace gcm {

/**
 * The PDE is like
 * \f[
 * d\vec{u}/dt + A_0 \partial{\vec{u}} / \partial{x_0}
 *             + ...
 *             + A_{Dim} \partial{\vec{u}} / \partial{x_{Dim}} = \vec{f}.
 * \f]
 * The class represents \f$ A_i \f$ with their eigensystems.
 * @tparam M size of PDE
 * @tparam Dim number of GcmMatrixes
 */
template<int TM, int Dim>
struct GcmMatrices {
	static const int M = TM;
	static const int D = Dim;
	typedef linal::Matrix<M, M> Matrix;

	/** Matrix in PDE along some direction with its eigensystem */
	struct GcmMatrix {
		/// matrix for some axis in PDE // TODO - skip this in Release mode?
		/// A = U1 * L * U  and so right eigenvectors are columns of the U1.
		Matrix A;
		/// matrix of left eigenvectors (aka eigenstrings)
		Matrix U;
		/// matrix of right eigenvectors
		Matrix U1;
		///< diagonal eigenvalue matrix
		linal::DiagonalMatrix<M> L;

		real getMaximalEigenvalue() const {
			real ans = 0;
			for (int i = 0; i < M; i++) {
				ans = fmax(ans, fabs(L(i, i)));
			}
			return ans;
		}
		
		void clear() {
			linal::clear(A);
			linal::clear(L);
			linal::clear(U1);
			linal::clear(U);
		}

	};

	GcmMatrix m[D];
	const GcmMatrix& operator()(const int s) const { return m[s]; }
	      GcmMatrix& operator()(const int s)       { return m[s]; }

	/** @return maximal in modulus eigenvalue of matrices from all directions */
	real getMaximalEigenvalue() const {
		real ans = 0;
		for (int i = 0; i < D; i++) {
			ans = fmax(ans, m[i].getMaximalEigenvalue());
		}
		return ans;
	}

	/** @throw gcm::Exception */
	void checkDecomposition(const real eps = EQUALITY_TOLERANCE) const {
		for (int s = 0; s < D; s++) {
			// traces
			assert_near(linal::trace(m[s].A), linal::trace(m[s].L), eps);
			// eigenvectors
			assert_true(linal::approximatelyEqual(m[s].A * m[s].U1, m[s].U1 * m[s].L, eps*1000));
			// eigenraws
			assert_true(linal::approximatelyEqual(m[s].U * m[s].A, m[s].L * m[s].U, eps*1000));
			// inverse matrices
			assert_true(linal::approximatelyEqual(m[s].U * m[s].U1, Matrix::Identity(), eps*100));
		}
	}
	
	void clear() {
		for (int s = 0; s < D; s++) {
			m[s].clear();
		}
	}

};


}

#endif // LIBGCM_GCMMATRICES_HPP
