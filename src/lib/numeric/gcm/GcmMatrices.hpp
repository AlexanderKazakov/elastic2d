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
 * d\vec{u}/dt + \mat{A_0}   \partial{\vec{u}}/\partial{x_0}
 *             + \mat{A_1}   \partial{\vec{u}}/\partial{x_1}
 *             + \mat{A_Dimensionality} \partial{\vec{u}}/\partial{x_Dimensionality} = \vec{f}.
 * The class represents \mat{A_i} with their eigensystems
 * @tparam M size of PDE
 * @tparam Dimensionality number of GcmMatrixes
 */
template<int TM, int Dimensionality>
struct GcmMatrices {
	static const int M = TM;
	static const int DIMENSIONALITY = Dimensionality;
	typedef linal::Matrix<M, M> Matrix;

	/** Matrix in PDE along some direction with its eigensystem */
	struct GcmMatrix {
		Matrix A;                   // matrix for some axis in PDE // TODO - skip this in
			                    // Release mode?
		/* A = U1 * L * U  and so right eigenvectors are columns of the U1 */
		Matrix U;                   // matrix of left eigenvectors (aka eigenstrings)
		Matrix U1;                  // matrix of right eigenvectors
		linal::DiagonalMatrix<M> L; // diagonal eigenvalue matrix

		real getMaximalEigenvalue() const {
			real ans = 0;
			for (int i = 0; i < M; i++) {
				ans = fmax(ans, fabs(L(i, i)));
			}
			return ans;
		}

	};

	GcmMatrix m[Dimensionality];

	/** @return maximal in modulus eigenvalue of matrices from all directions */
	real getMaximalEigenvalue() const {
		real ans = 0;
		for (int i = 0; i < Dimensionality; i++) {
			ans = fmax(ans, m[i].getMaximalEigenvalue());
		}
		return ans;
	}

	void checkDecomposition() const {
		checkTraces();
		checkLeftEigenvectors();
		checkRightEigenvectors();
		checkInverseMatrices();
	}

	void checkTraces() const {
		for (int s = 0; s < Dimensionality; s++) {
			assert_near(m[s].A.trace(), m[s].L.trace(), EQUALITY_TOLERANCE);
		}
	}

	void checkLeftEigenvectors() const {
		for (int s = 0; s < Dimensionality; s++) {
			Matrix AU1 = m[s].A * m[s].U1;
			Matrix U1L = m[s].U1 * m[s].L;
			for (int i = 0; i < M; i++) {
				for (int j = 0; j < M; j++) {
					assert_near(AU1(i, j), U1L(i, j), EQUALITY_TOLERANCE);
				}
			}
		}
	}

	void checkRightEigenvectors() const {
		for (int s = 0; s < Dimensionality; s++) {
			Matrix UA = m[s].U * m[s].A;
			Matrix LU = m[s].L * m[s].U;
			for (int i = 0; i < M; i++) {
				for (int j = 0; j < M; j++) {
					assert_near(UA(i, j), LU(i, j), EQUALITY_TOLERANCE);
				}
			}
		}
	}

	void checkInverseMatrices() const {
		for (int s = 0; s < Dimensionality; s++) {
			Matrix UU1 = m[s].U * m[s].U1;
			for (int i = 0; i < M; i++) {
				for (int j = 0; j < M; j++) {
					assert_near(UU1(i, j), (i == j), EQUALITY_TOLERANCE);
				}
			}
		}
	}

};

}

#endif // LIBGCM_GCMMATRICES_HPP
