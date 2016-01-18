#ifndef LIBGCM_ORTHOTROPICELASTIC3DGCMMATRICES_HPP
#define LIBGCM_ORTHOTROPICELASTIC3DGCMMATRICES_HPP

#include "lib/gcm_matrices/GcmMatrices.hpp"

namespace gcm {
	class OrthotropicElastic3DGcmMatrices : public GcmMatrices<9, 3> {
	public:
		real rho;
		real c11;
		real c12;
		real c13;
		real c22;
		real c23;
		real c33;
		real c44;
		real c55;
		real c66;

		/**
		 * Map between type of wave and corresponding to that type column in matrix U1
		 */
		static const std::map<Waves::T, int/* number of column in U1 */> WAVE_COLUMNS;

		OrthotropicElastic3DGcmMatrices(const real &rho, const real &lambda, const real &mu);
	};
}


#endif // LIBGCM_ORTHOTROPICELASTIC3DGCMMATRICES_HPP
