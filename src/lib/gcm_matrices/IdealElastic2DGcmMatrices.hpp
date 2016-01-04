#ifndef LIBGCM_IDEALELASTIC2DGCMMATRICES_HPP
#define LIBGCM_IDEALELASTIC2DGCMMATRICES_HPP

#include "lib/gcm_matrices/GcmMatrices.hpp"

namespace gcm {
	class IdealElastic2DGcmMatrices : public GcmMatrices<5, 2> {
	public:
		real rho;
		real lambda;
		real mu;

		IdealElastic2DGcmMatrices(const real &rho, const real &lambda, const real &mu);
	};
}


#endif // LIBGCM_IDEALELASTIC2DGCMMATRICES_HPP
