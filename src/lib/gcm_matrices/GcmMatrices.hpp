#ifndef LIBGCM_MATRIX_HPP
#define LIBGCM_MATRIX_HPP

#include "lib/linal/Vector.hpp"
#include "lib/util/Types.hpp"
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
	};

	/**
	 * The PDE is like
	 * d\vec{u}/dt + \mat{Ax} \partial{\vec{u}}/\partial{x}
	 *             + \mat{Ay} \partial{\vec{u}}/\partial{y}
	 *             + \mat{Az} \partial{\vec{u}}/\partial{z} = \vec{f}.
	 * The class represents \mat{Ax}, \mat{Ay} \mat{Az} with their eigensystems for a concrete model and material.
	 * @tparam M size of vec{u}
	 * @tparam Dimensionality number of dimensions, usually from 1 to 3
	 */
	template<int M, int Dimensionality>
	class GcmMatrices {
	public:
		static const int DIMENSIONALITY = Dimensionality; /// number of dimensions (stages)
		/** @return GcmMatrix along specified direction */
		const GcmMatrix<M>&A(const int direction) const {
			return m[direction];
		};
	protected:
		// TODO - logic around random basis here
		GcmMatrix<M> m[DIMENSIONALITY];
	};
}

#endif //LIBGCM_MATRIX_HPP
