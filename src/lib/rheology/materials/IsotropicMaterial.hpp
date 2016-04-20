#ifndef LIBGCM_ISOTROPICMATERIAL_HPP
#define LIBGCM_ISOTROPICMATERIAL_HPP

#include <lib/rheology/gcm_matrices/GcmMatrices.hpp>
#include <lib/rheology/variables/variables.hpp>

namespace gcm {
	struct Statement;
	class OrthotropicMaterial;
	
	class IsotropicMaterial {
	public:
		real rho = 0; // density

		real lambda = 0; // Lame parameters
		real mu = 0;

		real yieldStrength = 0; // plasticity parameter

		real continualDamageParameter = 0; // parameter in continual damage equation

		IsotropicMaterial(const IsotropicMaterial& other) = default;
		IsotropicMaterial(const OrthotropicMaterial&) {
			THROW_UNSUPPORTED("This is duct tape");
		}

		IsotropicMaterial(const real rho_ = 0, const real lambda_ = 0, const real mu_ = 0,
		                  const real yieldStrength_ = 0, const real continualDamageParameter_ = 0);
		
		void initialize(const Statement& statement);

		/** Fill in gcm matrices */
		void constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<1>, IsotropicMaterial>& m) const;
		void constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<2>, IsotropicMaterial>& m) const;
		void constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<3>, IsotropicMaterial>& m) const;

		/** For testing purposes */
		static IsotropicMaterial generateRandomMaterial();

	};
}

#endif // LIBGCM_ISOTROPICMATERIAL_HPP
