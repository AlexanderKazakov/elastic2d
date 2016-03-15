#ifndef LIBGCM_ISOTROPICMATERIAL_HPP
#define LIBGCM_ISOTROPICMATERIAL_HPP

#include <lib/rheology/materials/AbstractMaterial.hpp>

namespace gcm {
	struct IsotropicMaterial : public AbstractMaterial {
		static const Materials::T ID;

		real rho = 0; // density

		real lambda = 0; // Lame parameters
		real mu = 0;

		real yieldStrength = 0; // plasticity parameter
		real continualDamageParameter = 0; // parameter in continual damage equation

		IsotropicMaterial(const IsotropicMaterial& other) = default;
		IsotropicMaterial(const real rho_ = 0, const real lambda_ = 0, const real mu_ = 0,
		                  const real yieldStrength_ = 0, const real continualDamageParameter_ = 0);
		
		/** For testing purposes */
		static IsotropicMaterial generateRandomMaterial();

	};
}

#endif // LIBGCM_ISOTROPICMATERIAL_HPP
