#ifndef LIBGCM_MATRIX_HPP
#define LIBGCM_MATRIX_HPP

#include "lib/linal/Vector.hpp"
#include "lib/util/Types.hpp"
#include "lib/config.hpp"

namespace gcm {
	/**
	 * Represents a matrix in PDE along some direction with its eigensystem
	 */
	class GcmMatrix {
	public:
		linal::Matrix A; // matrix for some axis in PDE
		/* A = U1 * L * U  and so right eigenvectors are columns of the U1 */
		linal::Matrix U; // matrix of left eigenvectors (aka eigenstrings)
		linal::Matrix U1; // matrix of right eigenvectors
		linal::Matrix L; // diagonal eigenvalue matrix
	};

	/**
	 * The PDE is:
	 * d\vec{u}/dt + \mat{Ax} \partial{\vec{u}}/\partial{x} + \mat{Ay} \partial{\vec{u}}/\partial{y} = \vec{f}.
	 * The class represents \mat{Ax} and \mat{Ay} with their eigensystems for a concrete material.
	 */
	class GcmMatrices2D {
		GcmMatrix Ax;
		GcmMatrix Ay;

	public:
		real rho;
		real lambda;
		real mu;

		GcmMatrices2D(const real &rho, const real &lambda, const real &mu);

		/** @return GcmMatrix along stage direction */
		inline const GcmMatrix &A(const int stage) const;
	};
}

#endif //LIBGCM_MATRIX_HPP
