#ifndef LIBGCM_IDEALELASTIC2DGCMMATRICES_HPP
#define LIBGCM_IDEALELASTIC2DGCMMATRICES_HPP

#include "lib/gcm_matrices/GcmMatrices.hpp"

namespace gcm {
	class IdealElastic2DGcmMatrices : public GcmMatrices<5, 2> {
	public:
		real rho;
		real lambda;
		real mu;

		/**
		 * Map between type of wave and corresponding to that type column in matrix U1
		 */
		static const std::map<Waves::WAVE, int /* number of column in U1 */> WAVE_COLUMNS;

		IdealElastic2DGcmMatrices(const real &rho, const real &lambda, const real &mu);
	};
}


#endif // LIBGCM_IDEALELASTIC2DGCMMATRICES_HPP
