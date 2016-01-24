#ifndef LIBGCM_GCMMATRICES_HPP
#define LIBGCM_GCMMATRICES_HPP

#include <map>

#include <lib/rheology/materials/OrthotropicMaterial.hpp>
#include <lib/rheology/materials/IsotropicMaterial.hpp>
#include <lib/linal/Linal.hpp>
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
	 * @tparam TM size of vec{u}
	 * @tparam Dimensionality number of dimensions, usually from 1 to 3
	 * @tparam TMaterial type of material used to construct the matrices
	 */
	template<int TM, int Dimensionality, class TMaterial>
	class GcmMatrices {
	public:
		typedef typename GcmMatrix<TM>::Matrix Matrix;
		typedef TMaterial Material; // type of material used to construct the matrices
		static const int M = TM; // size of corresponding vector
		static const int DIMENSIONALITY = Dimensionality; // number of dimensions (stages)

		/** Map between type of wave and corresponding to that type column in matrix U1 */
		static const std::map<Waves::T, int/* index of column in U1 */> WAVE_COLUMNS;

		GcmMatrices(const TMaterial& material) : material(material) {
			material.constructGcmMatrices(*this);
		};

		/** @return GcmMatrix along specified direction */
		const GcmMatrix<TM>& A(const int direction) const {
			return m[direction];
		};

		/** @return maximal in modulus eigenvalue of matrices from all directions */
		real getMaximalEigenvalue() const;

		TMaterial getMaterial() const { return material; };

	protected:

		TMaterial material;

		// TODO - logic around random basis here and/or in Grid?
		GcmMatrix<TM> m[Dimensionality];

		friend TMaterial;
	};
}

#endif // LIBGCM_GCMMATRICES_HPP
