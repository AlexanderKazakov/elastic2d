#ifndef LIBGCM_ISOTROPICMATERIAL_HPP
#define LIBGCM_ISOTROPICMATERIAL_HPP

#include <libgcm/rheology/materials/AbstractMaterial.hpp>

namespace gcm {
struct IsotropicMaterial : public AbstractMaterial {
	static const Materials::T Type = Materials::T::ISOTROPIC;

	real rho = 0;                      ///< density

	real lambda = 0;                   ///< Lame parameters
	real mu = 0;

	real yieldStrength = 0;            ///< plasticity parameter
	real continualDamageParameter = 0; ///< parameter in continual damage equation
	real tau0 = 0;                     ///< viscosity parameter (decay time)

	IsotropicMaterial(const IsotropicMaterial& other) = default;
	IsotropicMaterial(const real rho_ = 0, const real lambda_ = 0, const real mu_ = 0,
			const real yieldStrength_ = 0, const real continualDamageParameter_ = 0,
			const int materialNumber_ = 0, const real tau0_ = 0);

	/** For testing purposes */
	static IsotropicMaterial generateRandomMaterial();

};
}

#endif // LIBGCM_ISOTROPICMATERIAL_HPP
