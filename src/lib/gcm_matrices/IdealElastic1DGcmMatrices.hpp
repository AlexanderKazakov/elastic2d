#ifndef LIBGCM_IDEALELASTIC1DGCMMATRICES_HPP
#define LIBGCM_IDEALELASTIC1DGCMMATRICES_HPP

#include "lib/gcm_matrices/GcmMatrices.hpp"

namespace gcm {
	class IdealElastic1DGcmMatrices : public GcmMatrices<2, 1> {
	public:
		real rho;
		real lambda;
		real mu;

		IdealElastic1DGcmMatrices(const real &rho, const real &lambda, const real &mu);
	};
}


#endif // LIBGCM_IDEALELASTIC1DGCMMATRICES_HPP
