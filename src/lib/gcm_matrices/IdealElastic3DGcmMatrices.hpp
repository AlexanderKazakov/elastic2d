#ifndef LIBGCM_IDEALELASTIC3DGCMMATRICES_HPP
#define LIBGCM_IDEALELASTIC3DGCMMATRICES_HPP

#include "lib/gcm_matrices/GcmMatrices.hpp"

namespace gcm {
	class IdealElastic3DGcmMatrices : public GcmMatrices<9, 3> {
	public:
		real rho;
		real lambda;
		real mu;

		/**
		 * Map between type of wave and corresponding to that type column in matrix U1
		 */
		static const std::map<Waves::WAVE, int/* number of column in U1 */> WAVE_COLUMNS;

		IdealElastic3DGcmMatrices(const real &rho, const real &lambda, const real &mu);
	};
}


#endif // LIBGCM_IDEALELASTIC3DGCMMATRICES_HPP
