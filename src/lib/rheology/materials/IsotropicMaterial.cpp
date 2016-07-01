#include <lib/rheology/materials/IsotropicMaterial.hpp>

using namespace gcm;

const Materials::T IsotropicMaterial::ID = Materials::T::ISOTROPIC;

IsotropicMaterial::
IsotropicMaterial(const real rho_, const real lambda_, const real mu_,
		const real yieldStrength_, const real continualDamageParameter_,
		const int materialNumber_) :
	rho(rho_), lambda(lambda_), mu(mu_), yieldStrength(yieldStrength_),
	continualDamageParameter(continualDamageParameter_) {
	
	materialNumber = materialNumber_;
}


IsotropicMaterial IsotropicMaterial::
generateRandomMaterial() {
	const real RHO_MAX = 100.0;
	const real RHO_MIN = 0.01;
	const real LAMBDA_MAX = 1e+6;
	const real LAMBDA_MIN = 1.0;
	const real MU_MAX = 1e+6;
	const real MU_MIN = 1.0;

	real rho = Utils::randomReal(RHO_MIN, RHO_MAX);
	real lambda = Utils::randomReal(LAMBDA_MIN, LAMBDA_MAX);
	real mu = Utils::randomReal(MU_MIN, MU_MAX);

	return IsotropicMaterial(rho, lambda, mu);
}


