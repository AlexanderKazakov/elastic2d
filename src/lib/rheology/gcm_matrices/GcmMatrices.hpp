#ifndef LIBGCM_GCMMATRICES_HPP
#define LIBGCM_GCMMATRICES_HPP

#include <map>

#include <lib/linal/LinalRoutines.hpp>
#include <lib/util/Types.hpp>
#include <lib/util/Concepts.hpp>
#include <lib/config.hpp>

namespace gcm {
	/**
	 * Represents a matrix in PDE along some direction with its eigensystem
	 * @tparam M size of vector in PDE
	 */
	template<int M>
	class GcmMatrix {
	public:
		typedef linal::Matrix<M,M> Matrix;
		Matrix A; // matrix for some axis in PDE
		/* A = U1 * L * U  and so right eigenvectors are columns of the U1 */
		Matrix U; // matrix of left eigenvectors (aka eigenstrings)
		Matrix U1; // matrix of right eigenvectors
		linal::DiagonalMatrix<M> L; // diagonal eigenvalue matrix

		real getMaximalEigenvalue() const;
	};

	/**
	 * The PDE is like
	 * d\vec{u}/dt + \mat{Ax} \partial{\vec{u}}/\partial{x}
	 *             + \mat{Ay} \partial{\vec{u}}/\partial{y}
	 *             + \mat{Az} \partial{\vec{u}}/\partial{z} = \vec{f}.
	 * The class represents \mat{Ax}, \mat{Ay} \mat{Az} with their eigensystems for a concrete model and material.
	 * @tparam TVariables variables in equation - \vec{u}
	 * @tparam TMaterial type of material used to construct the matrices
	 */
	template<typename TVariables, class TMaterial>
	class GcmMatrices {
	public:
		typedef TMaterial Material; // type of material used to construct the matrices
		static const int M = TVariables::SIZE; // size of corresponding vector
		static const int DIMENSIONALITY = TVariables::DIMENSIONALITY; // number of dimensions (stages)
		typedef typename GcmMatrix<M>::Matrix Matrix;

		/** Map between type of wave and corresponding to that type column in matrix U1 */
		static const std::map<Waves::T, int/* index of column in U1 */> WAVE_COLUMNS;

		GcmMatrices(const TMaterial& _material) : material(_material) {
			material.constructGcmMatrices(*this);
		};

		/** @return GcmMatrix along specified direction */
		const GcmMatrix<M>& A(const int direction) const {
			return m[direction];
		};

		/** @return maximal in modulus eigenvalue of matrices from all directions */
		real getMaximalEigenvalue() const;

		Material getMaterial() const { return material; };

	protected:
		Material material;
		GcmMatrix<M> m[DIMENSIONALITY];

		friend Material;
	};
}

#endif // LIBGCM_GCMMATRICES_HPP
