#include "IdealElastic1DGcmMatrices.hpp"

using namespace gcm;

IdealElastic1DGcmMatrices::IdealElastic1DGcmMatrices(const real& rho, const real& lambda, const real& mu) :
		rho(rho), lambda(lambda), mu(mu) {

	real E = mu * (3 * lambda + 2 * mu) / (lambda + mu); // Young's modulus

	m[0].A(0, 0) = 0;              m[0].A(0, 1) = - 1.0 / rho;
	m[0].A(1, 0) = - E;            m[0].A(1, 1) = 0;

	m[0].L.createDiagonal( { sqrt(E / rho), - sqrt(E / rho) } );

	m[0].U = linal::Matrix<2,2>( { - 0.5, 1.0 / (2 * sqrt(E * rho)),
	                                 0.5, 1.0 / (2 * sqrt(E * rho)) } );

	m[0].U1 = linal::Matrix<2,2>( { - 1.0,         1.0,
	                                sqrt(E * rho), sqrt(E * rho) } );
}