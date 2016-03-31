#include <lib/rheology/materials/IsotropicMaterial.hpp>

using namespace gcm;

const Materials::T IsotropicMaterial::ID = Materials::T::ISOTROPIC;

IsotropicMaterial::
IsotropicMaterial(const real rho_, const real lambda_, const real mu_,
                  const real yieldStrength_, const real continualDamageParameter_) :
	rho(rho_), lambda(lambda_), mu(mu_), yieldStrength(yieldStrength_),
	continualDamageParameter(continualDamageParameter_) { }

IsotropicMaterial IsotropicMaterial::
generateRandomMaterial() {
	const real RHO_MAX = 100.0;
	const real RHO_MIN = 0.01;
	const real LAMBDA_MAX = 1e+6;
	const real LAMBDA_MIN = 1.0;
	const real MU_MAX = 1e+6;
	const real MU_MIN = 1.0;

	real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
	real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
	real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;

	return IsotropicMaterial(rho, lambda, mu);
}


