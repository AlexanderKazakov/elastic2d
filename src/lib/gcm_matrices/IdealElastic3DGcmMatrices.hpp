#ifndef LIBGCM_IDEALELASTIC3DGCMMATRICES_HPP
#define LIBGCM_IDEALELASTIC3DGCMMATRICES_HPP

#include "lib/gcm_matrices/GcmMatrices.hpp"

namespace gcm {
	class IdealElastic3DGcmMatrices : public GcmMatrices<9, 3> {
	public:
		real rho;
		real lambda;
		real mu;

		IdealElastic3DGcmMatrices(const real &rho, const real &lambda, const real &mu);
	};
}


#endif // LIBGCM_IDEALELASTIC3DGCMMATRICES_HPP
