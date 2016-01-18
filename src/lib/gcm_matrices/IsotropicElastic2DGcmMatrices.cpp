#include "lib/gcm_matrices/IsotropicElastic2DGcmMatrices.hpp"

using namespace gcm;

const std::map<Waves::T, int /* number of column in U1 */> IsotropicElastic2DGcmMatrices::WAVE_COLUMNS = {
		{Waves::T::P_FORWARD,   1},
		{Waves::T::P_BACKWARD,  0},
		{Waves::T::S1_FORWARD,  3},
		{Waves::T::S1_BACKWARD, 2}
};

IsotropicElastic2DGcmMatrices::IsotropicElastic2DGcmMatrices(const real &rho, const real &lambda, const real &mu) :
		rho(rho), lambda(lambda), mu(mu) {

	m[0].A = linal::Matrix<5, 5>({0, 0, -1.0 / rho, 0, 0,
	                              0, 0, 0, -1.0 / rho, 0,
	                              -lambda - 2.0 * mu, 0, 0, 0, 0,
	                              0, -mu, 0, 0, 0,
	                              -lambda, 0, 0, 0, 0});
	
	m[0].L.createDiagonal(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), sqrt(mu / rho), 0});

	m[0].U = linal::Matrix<5, 5>({1.0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                              1.0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                              0, 1.0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                              0, 1.0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                              0, 0, 1.0 / (lambda + 2 * mu), 0, -1.0 / lambda});

	m[0].U1 = linal::Matrix<5, 5>({0.5, 0.5, 0, 0, 0,
	                               0, 0, 0.5, 0.5, 0,
	                               0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0,
	                               0, 0, 0.5 * sqrt(mu * rho), -0.5 * sqrt(mu * rho), 0,
	                               (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, -lambda});


	m[1].A = linal::Matrix<5, 5>({0, 0, 0, -1.0 / rho, 0,
	                              0, 0, 0, 0, -1.0 / rho,
	                              0, -lambda, 0, 0, 0,
	                              -mu, 0, 0, 0, 0,
	                              0, -lambda - 2.0 * mu, 0, 0, 0});

	m[1].L.createDiagonal(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), sqrt(mu / rho), 0});

	m[1].U = linal::Matrix<5, 5>({0, 1.0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                              0, 1.0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                              1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                              1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                              0, 0, 1.0, 0, -(1.0 * lambda) / (lambda + 2 * mu)});

	m[1].U1 = linal::Matrix<5, 5>({0, 0, 0.5, 0.5, 0,
	                               0.5, 0.5, 0, 0, 0,
	                               (0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), -(0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), 0, 0, 1.0,
	                               0, 0, 0.5 * sqrt(mu * rho), -0.5 * sqrt(mu * rho), 0,
	                               0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0});
}