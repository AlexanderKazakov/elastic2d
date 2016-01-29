#include <lib/rheology/materials/IsotropicMaterial.hpp>

using namespace gcm;


IsotropicMaterial::IsotropicMaterial() { }

IsotropicMaterial::IsotropicMaterial(const real _rho, const real _lambda, const real _mu) :
		rho(_rho), lambda(_lambda), mu(_mu)  { }

IsotropicMaterial::IsotropicMaterial(const real _rho, const real _lambda, const real _mu, const real _yieldStrength) :
		rho(_rho), lambda(_lambda), mu(_mu), yieldStrength(_yieldStrength) { }

void IsotropicMaterial::constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<1>, IsotropicMaterial> &m) const {
	real E = mu * (3 * lambda + 2 * mu) / (lambda + mu); // Young's modulus

	m.m[0].A.initialize({0.0, -1.0 / rho,
	                     -E, 0.0});

	m.m[0].L.initialize({sqrt(E / rho), -sqrt(E / rho)});

	m.m[0].U.initialize({-0.5, 1.0 / (2 * sqrt(E * rho)),
	                     0.5, 1.0 / (2 * sqrt(E * rho))});

	m.m[0].U1.initialize({-1.0, 1.0,
	                      sqrt(E * rho), sqrt(E * rho)});
}

void IsotropicMaterial::constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<2>, IsotropicMaterial> &m) const {
	// TODO - actually we can use orthotropic material here
	m.m[0].A.initialize({0, 0, -1.0 / rho, 0, 0,
	                     0, 0, 0, -1.0 / rho, 0,
	                     -lambda - 2.0 * mu, 0, 0, 0, 0,
	                     0, -mu, 0, 0, 0,
	                     -lambda, 0, 0, 0, 0});

	m.m[0].L.initialize(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), sqrt(mu / rho), 0});

	m.m[0].U.initialize({1.0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                     1.0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                     0, 1.0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                     0, 1.0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                     0, 0, 1.0 / (lambda + 2 * mu), 0, -1.0 / lambda});

	m.m[0].U1.initialize({0.5, 0.5, 0, 0, 0,
	                      0, 0, 0.5, 0.5, 0,
	                      0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0,
	                      0, 0, 0.5 * sqrt(mu * rho), -0.5 * sqrt(mu * rho), 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                      -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, -lambda});


	m.m[1].A.initialize({0, 0, 0, -1.0 / rho, 0,
	                     0, 0, 0, 0, -1.0 / rho,
	                     0, -lambda, 0, 0, 0,
	                     -mu, 0, 0, 0, 0,
	                     0, -lambda - 2.0 * mu, 0, 0, 0});

	m.m[1].L.initialize(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), sqrt(mu / rho), 0});

	m.m[1].U.initialize({0, 1.0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                     0, 1.0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                     1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                     1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                     0, 0, 1.0, 0, -(1.0 * lambda) / (lambda + 2 * mu)});

	m.m[1].U1.initialize({0, 0, 0.5, 0.5, 0,
	                      0.5, 0.5, 0, 0, 0,
	                      (0.5 * lambda) / sqrt((lambda + 2 * mu) / rho),
	                      -(0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), 0, 0, 1.0,
	                      0, 0, 0.5 * sqrt(mu * rho), -0.5 * sqrt(mu * rho), 0,
	                      0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0});
}

void IsotropicMaterial::constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<3>, IsotropicMaterial> &m) const {
	// TODO - actually we can use orthotropic material here
	m.m[0].A.initialize({0, 0, 0, -1.0 / rho, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
	                     -lambda - 2 * mu, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, -mu, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -mu, 0, 0, 0, 0, 0, 0,
	                     -lambda, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, 0,
	                     -lambda, 0, 0, 0, 0, 0, 0, 0, 0});

	m.m[0].L.initialize(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), -sqrt(mu / rho),
			 sqrt(mu / rho), sqrt(mu / rho), 0, 0, 0});

	m.m[0].U.initialize({1.0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0, 0, 0, 0,
	                     1.0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0, 0, 0, 0,
	                     0, 1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
	                     0, 0, 1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0,
	                     0, 1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
	                     0, 0, 1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0,
	                     0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 0, 1.0, 0, 0.0,
	                     0, 0, 0, 0, 0, 0, 0, 1.0, 0,
	                     0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 0, 0, 0, 1.0});

	m.m[0].U1.initialize({0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0.5, 0, 0.5, 0, 0, 0, 0,
	                      0, 0, 0, 0.5, 0, 0.5, 0, 0, 0,
	                      0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0, 0, 0, 0,
	                      0,
	                      0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
	                      0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                      -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 1, 0, 0,
	                      0, 0, 0, 0, 0, 0, 0, 1, 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                      -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1});


	m.m[1].A.initialize({0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, -1.0 / rho, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
	                     0, -lambda, 0, 0, 0, 0, 0, 0, 0,
	                     -mu, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, -lambda - 2 * mu, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -mu, 0, 0, 0, 0, 0, 0,
	                     0, -lambda, 0, 0, 0, 0, 0, 0, 0});

	m.m[1].L.initialize(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), -sqrt(mu / rho),
			 sqrt(mu / rho), sqrt(mu / rho), 0, 0, 0});

	m.m[1].U.initialize({0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                     0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                     1.0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
	                     0, 0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                     1.0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
	                     0, 0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                     0, 0, 0, 1.0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 0.0,
	                     0, 0, 0, 0, 0, 1.0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 1.0});

	m.m[1].U1.initialize({0, 0, 0.5, 0, 0.5, 0, 0, 0, 0,
	                      0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0, 0.5, 0, 0.5, 0, 0, 0,
	                      (0.5 * lambda) / sqrt((lambda + 2 * mu) / rho),
	                      -(0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), 0, 0, 0, 0, 1, 0, 0,
	                      0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
	                      0, 0, 0, 0, 0, 0, 0, 1, 0,
	                      0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0,
	                      0, 0, 0, 0, 0, 0,
	                      0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                      -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1});


	m.m[2].A.initialize({0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, -1.0 / rho,
	                     0, 0, -lambda, 0, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, 0,
	                     -mu, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -lambda, 0, 0, 0, 0, 0, 0,
	                     0, -mu, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -lambda - 2 * mu, 0, 0, 0, 0, 0, 0});

	m.m[2].L.initialize(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), -sqrt(mu / rho),
			 sqrt(mu / rho), sqrt(mu / rho), 0, 0, 0});

	m.m[2].U.initialize({0, 0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                     0, 0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                     1.0, 0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0,
	                     0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                     1.0, 0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0,
	                     0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                     0, 0, 0, 1.0, 0, 0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu),
	                     0, 0, 0, 0, 1.0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 1.0, 0, -(1.0 * lambda) / (lambda + 2 * mu)});

	m.m[2].U1.initialize({0, 0, 0.5, 0, 0.5, 0, 0, 0, 0,
	                      0, 0, 0, 0.5, 0, 0.5, 0, 0, 0,
	                      0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                      -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 1, 0, 0,
	                      0, 0, 0, 0, 0, 0, 0, 1, 0,
	                      0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                      -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1,
	                      0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
	                      0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0,
	                      0, 0, 0, 0, 0, 0});
}

IsotropicMaterial IsotropicMaterial::generateRandomMaterial() {
	static const real RHO_MAX = 100.0;
	static const real RHO_MIN = 0.01;
	static const real LAMBDA_MAX = 1e+6;
	static const real LAMBDA_MIN = 1.0;
	static const real MU_MAX = 1e+6;
	static const real MU_MIN = 1.0;

	real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
	real lambda = ((LAMBDA_MAX - LAMBDA_MIN) * rand()) / RAND_MAX + LAMBDA_MIN;
	real mu = ((MU_MAX - MU_MIN) * rand()) / RAND_MAX + MU_MIN;

	return IsotropicMaterial(rho, lambda, mu);
}
