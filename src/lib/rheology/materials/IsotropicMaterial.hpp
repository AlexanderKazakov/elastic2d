#ifndef LIBGCM_ISOTROPICMATERIAL_HPP
#define LIBGCM_ISOTROPICMATERIAL_HPP

#include <lib/rheology/gcm_matrices/GcmMatrices.hpp>
#include <lib/rheology/variables/VelocitySigmaVariables.hpp>

namespace gcm {

	/**
	 * Isotropic elastic materials
	 */
	class IsotropicMaterial {
	public:
		real rho;

		real lambda;
		real mu;

		IsotropicMaterial();
		IsotropicMaterial(const real& _rho, const real& _lambda, const real& _mu);

		/** Fill in gcm matrices */
		void constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<1>, IsotropicMaterial>& m) const;
		void constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<2>, IsotropicMaterial>& m) const;
		void constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<3>, IsotropicMaterial>& m) const;

		/** For testing purposes */
		static IsotropicMaterial generateRandomMaterial();

	};
}

#endif // LIBGCM_ISOTROPICMATERIAL_HPP
