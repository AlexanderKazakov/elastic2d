#ifndef LIBGCM_IDEALELASTIC1DGCMMATRICES_HPP
#define LIBGCM_IDEALELASTIC1DGCMMATRICES_HPP

#include "lib/gcm_matrices/GcmMatrices.hpp"

namespace gcm {
	class Elastic1DGcmMatrices : public GcmMatrices<2, 1> {
	public:
		real rho;
		real lambda;
		real mu;

		/**
		 * Map between type of wave and corresponding to that type column in matrix U1
		 */
		static const std::map<Waves::T, int/* number of column in U1 */> WAVE_COLUMNS;

		Elastic1DGcmMatrices(const real &rho, const real &lambda, const real &mu);
	};
}


#endif // LIBGCM_IDEALELASTIC1DGCMMATRICES_HPP
