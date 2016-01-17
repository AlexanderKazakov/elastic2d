#ifndef LIBGCM_MATRIX_HPP
#define LIBGCM_MATRIX_HPP

#include <map>

#include "lib/linal/Linal.hpp"
#include "lib/util/Types.hpp"
#include "lib/util/Concepts.hpp"
#include "lib/config.hpp"

namespace gcm {
	/**
	 * Represents a matrix in PDE along some direction with its eigensystem
	 * @tparam M size of vector in PDE
	 */
	template<int M>
	class GcmMatrix {
	public:
		linal::Matrix<M,M> A; // matrix for some axis in PDE
		/* A = U1 * L * U  and so right eigenvectors are columns of the U1 */
		linal::Matrix<M,M> U; // matrix of left eigenvectors (aka eigenstrings)
		linal::Matrix<M,M> U1; // matrix of right eigenvectors
		linal::Matrix<M,M> L; // diagonal eigenvalue matrix

		real getMaximalEigenvalue() const {
			real ans = 0.0;
			for (int i = 0; i < M; i++) {
				ans = fmax(ans, fabs(L(i, i)));
			}
			return ans;
		};
	};

	/**
	 * The PDE is like
	 * d\vec{u}/dt + \mat{Ax} \partial{\vec{u}}/\partial{x}
	 *             + \mat{Ay} \partial{\vec{u}}/\partial{y}
	 *             + \mat{Az} \partial{\vec{u}}/\partial{z} = \vec{f}.
	 * The class represents \mat{Ax}, \mat{Ay} \mat{Az} with their eigensystems for a concrete model and material.
	 * @tparam TM size of vec{u}
	 * @tparam Dimensionality number of dimensions, usually from 1 to 3
	 */
	template<int TM, int Dimensionality>
	class GcmMatrices {
	public:
		static const int M = TM; // size of corresponding vector
		static const int DIMENSIONALITY = Dimensionality; // number of dimensions (stages)

		/** @return GcmMatrix along specified direction */
		const GcmMatrix<TM>& A(const int direction) const {
			return m[direction];
		};

		/** @return maximal in modulus eigenvalue of matrices from all directions */
		real getMaximalEigenvalue() const {
			real ans = 0.0;
			for (int i = 0; i < Dimensionality; i++) {
				ans = fmax(ans, A(i).getMaximalEigenvalue());
			}
			return ans;
		};

	protected:
		// TODO - logic around random basis here and/or in Grid?
		GcmMatrix<TM> m[Dimensionality];
	};
}

#endif //LIBGCM_MATRIX_HPP
