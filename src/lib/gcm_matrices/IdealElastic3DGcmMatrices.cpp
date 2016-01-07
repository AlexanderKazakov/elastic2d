#include "IdealElastic3DGcmMatrices.hpp"

using namespace gcm;

IdealElastic3DGcmMatrices::IdealElastic3DGcmMatrices(const real& rho, const real& lambda, const real& mu) :
		rho(rho), lambda(lambda), mu(mu) {

	for (int i = 0; i < DIMENSIONALITY; i++) {
		clear(m[i].A);
		clear(m[i].L);
		clear(m[i].U);
		clear(m[i].U1);
	}

	m[0].A(0, 3) = -1.0 / rho;
	m[0].A(1, 4) = -1.0 / rho;
	m[0].A(2, 5) = -1.0 / rho;
	m[0].A(3, 0) = -lambda - 2 * mu;
	m[0].A(4, 1) = -mu;
	m[0].A(5, 2) = -mu;
	m[0].A(6, 0) = -lambda;
	m[0].A(8, 0) = -lambda;

	m[0].L(0, 0) = sqrt(mu / rho);
	m[0].L(1, 1) = sqrt(mu / rho);
	m[0].L(2, 2) = -sqrt(mu / rho);
	m[0].L(3, 3) = -sqrt(mu / rho);
	m[0].L(7, 7) = sqrt((lambda + 2 * mu) / rho);
	m[0].L(8, 8) = -sqrt((lambda + 2 * mu) / rho);

	m[0].U(0, 1) = -sqrt(rho * mu)*0.5;
	m[0].U(0, 4) = 0.5;
	m[0].U(1, 2) = -sqrt(rho * mu)*0.5;
	m[0].U(1, 5) = 0.5;
	m[0].U(2, 1) = sqrt(rho * mu)*0.5;
	m[0].U(2, 4) = 0.5;
	m[0].U(3, 2) = sqrt(rho * mu)*0.5;
	m[0].U(3, 5) = 0.5;
	m[0].U(4, 7) = 1;
	m[0].U(5, 3) = -lambda / (lambda + 2 * mu);
	m[0].U(5, 6) = 1;
	m[0].U(6, 3) = -lambda / (lambda + 2 * mu);
	m[0].U(6, 8) = 1;
	m[0].U(7, 0) = -0.5 * lambda * sqrt(rho / (lambda + 2 * mu));
	m[0].U(7, 3) = 0.5 * lambda / (lambda + 2 * mu);
	m[0].U(8, 0) = 0.5 * lambda * sqrt(rho / (lambda + 2 * mu));
	m[0].U(8, 3) = 0.5 * lambda / (lambda + 2 * mu);

	m[0].U1(0, 7) = -sqrt((lambda + 2 * mu) / rho) / lambda;
	m[0].U1(0, 8) = sqrt((lambda + 2 * mu) / rho) / lambda;
	m[0].U1(1, 0) = -1.0 / sqrt(rho * mu);
	m[0].U1(1, 2) = 1.0 / sqrt(rho * mu);
	m[0].U1(2, 1) = -1.0 / sqrt(rho * mu);
	m[0].U1(2, 3) = 1.0 / sqrt(rho * mu);
	m[0].U1(3, 7) = (lambda + 2 * mu) / lambda;
	m[0].U1(3, 8) = (lambda + 2 * mu) / lambda;
	m[0].U1(4, 0) = 1;
	m[0].U1(4, 2) = 1;
	m[0].U1(5, 1) = 1;
	m[0].U1(5, 3) = 1;
	m[0].U1(6, 5) = 1;
	m[0].U1(6, 7) = 1;
	m[0].U1(6, 8) = 1;
	m[0].U1(7, 4) = 1;
	m[0].U1(8, 6) = 1;
	m[0].U1(8, 7) = 1;
	m[0].U1(8, 8) = 1;


	m[1].A(0, 4) = -1.0 / rho;
	m[1].A(1, 6) = -1.0 / rho;
	m[1].A(2, 7) = -1.0 / rho;
	m[1].A(3, 1) = -lambda;
	m[1].A(4, 0) = -mu;
	m[1].A(6, 1) = -lambda - 2 * mu;
	m[1].A(7, 2) = -mu;
	m[1].A(8, 1) = -lambda;
	
	m[1].L(0, 0) = sqrt((lambda + 2 * mu) / rho);
	m[1].L(1, 1) = -sqrt((lambda + 2 * mu) / rho);
	m[1].L(2, 2) = sqrt(mu / rho);
	m[1].L(3, 3) = sqrt(mu / rho);
	m[1].L(4, 4) = -sqrt(mu / rho);
	m[1].L(5, 5) = -sqrt(mu / rho);

	m[1].U(0, 1) = -0.5 * lambda * sqrt(rho / (lambda + 2 * mu));
	m[1].U(0, 6) = 0.5 * lambda / (lambda + 2 * mu);
	m[1].U(1, 1) = 0.5 * lambda * sqrt(rho / (lambda + 2 * mu));
	m[1].U(1, 6) = 0.5 * lambda / (lambda + 2 * mu);
	m[1].U(2, 0) = -sqrt(rho * mu)*0.5;
	m[1].U(2, 4) = 0.5;
	m[1].U(3, 2) = -sqrt(rho * mu)*0.5;
	m[1].U(3, 7) = 0.5;
	m[1].U(4, 0) = sqrt(rho * mu)*0.5;
	m[1].U(4, 4) = 0.5;
	m[1].U(5, 2) = sqrt(rho * mu)*0.5;
	m[1].U(5, 7) = 0.5;
	m[1].U(6, 3) = 1;
	m[1].U(6, 6) = -lambda / (lambda + 2 * mu);
	m[1].U(7, 5) = 1;
	m[1].U(8, 6) = -lambda / (lambda + 2 * mu);
	m[1].U(8, 8) = 1;

	m[1].U1(0, 2) = -1.0 / sqrt(rho * mu);
	m[1].U1(0, 4) = 1.0 / sqrt(rho * mu);
	m[1].U1(1, 0) = -sqrt((lambda + 2 * mu) / rho) / lambda;
	m[1].U1(1, 1) = sqrt((lambda + 2 * mu) / rho) / lambda;
	m[1].U1(2, 3) = -1.0 / sqrt(rho * mu);
	m[1].U1(2, 5) = 1.0 / sqrt(rho * mu);
	m[1].U1(3, 0) = 1;
	m[1].U1(3, 1) = 1;
	m[1].U1(3, 6) = 1;
	m[1].U1(4, 2) = 1;
	m[1].U1(4, 4) = 1;
	m[1].U1(5, 7) = 1;
	m[1].U1(6, 0) = (lambda + 2 * mu) / lambda;
	m[1].U1(6, 1) = (lambda + 2 * mu) / lambda;
	m[1].U1(7, 3) = 1;
	m[1].U1(7, 5) = 1;
	m[1].U1(8, 0) = 1;
	m[1].U1(8, 1) = 1;
	m[1].U1(8, 8) = 1;

	
	m[2].A(0, 5) = -1.0 / rho;
	m[2].A(1, 7) = -1.0 / rho;
	m[2].A(2, 8) = -1.0 / rho;
	m[2].A(3, 2) = -lambda;
	m[2].A(5, 0) = -mu;
	m[2].A(6, 2) = -lambda;
	m[2].A(7, 1) = -mu;
	m[2].A(8, 2) = -lambda - 2 * mu;
	
	m[2].L(3, 3) = sqrt((lambda + 2 * mu) / rho);
	m[2].L(4, 4) = -sqrt((lambda + 2 * mu) / rho);
	m[2].L(5, 5) = sqrt(mu / rho);
	m[2].L(6, 6) = sqrt(mu / rho);
	m[2].L(7, 7) = -sqrt(mu / rho);
	m[2].L(8, 8) = -sqrt(mu / rho);

	m[2].U(0, 6) = 1;
	m[2].U(0, 8) = -lambda / (lambda + 2 * mu);
	m[2].U(1, 3) = 1;
	m[2].U(1, 8) = -lambda / (lambda + 2 * mu);
	m[2].U(2, 4) = 1;
	m[2].U(3, 2) = -0.5 * lambda * sqrt(rho / (lambda + 2 * mu));
	m[2].U(3, 8) = 0.5 * lambda / (lambda + 2 * mu);
	m[2].U(4, 2) = 0.5 * lambda * sqrt(rho / (lambda + 2 * mu));
	m[2].U(4, 8) = 0.5 * lambda / (lambda + 2 * mu);
	m[2].U(5, 0) = 0.5;
	m[2].U(5, 5) = -0.5 / sqrt(rho * mu);
	m[2].U(6, 1) = 0.5;
	m[2].U(6, 7) = -0.5 / sqrt(rho * mu);
	m[2].U(7, 0) = 0.5;
	m[2].U(7, 5) = 0.5 / sqrt(rho * mu);
	m[2].U(8, 1) = 0.5;
	m[2].U(8, 7) = 0.5 / sqrt(rho * mu);

	m[2].U1(0, 5) = 1;
	m[2].U1(0, 7) = 1;
	m[2].U1(1, 6) = 1;
	m[2].U1(1, 8) = 1;
	m[2].U1(2, 3) = -sqrt((lambda + 2 * mu) / rho) / lambda;
	m[2].U1(2, 4) = sqrt((lambda + 2 * mu) / rho) / lambda;
	m[2].U1(3, 1) = 1;
	m[2].U1(3, 3) = 1;
	m[2].U1(3, 4) = 1;
	m[2].U1(4, 2) = 1;
	m[2].U1(5, 5) = -sqrt(rho * mu);
	m[2].U1(5, 7) = sqrt(rho * mu);
	m[2].U1(6, 0) = 1;
	m[2].U1(6, 3) = 1;
	m[2].U1(6, 4) = 1;
	m[2].U1(7, 6) = -sqrt(rho * mu);
	m[2].U1(7, 8) = sqrt(rho * mu);
	m[2].U1(8, 3) = (lambda + 2 * mu) / lambda;
	m[2].U1(8, 4) = (lambda + 2 * mu) / lambda;

}