#include <libgcm/rheology/materials/OrthotropicMaterial.hpp>
#include <libgcm/rheology/materials/IsotropicMaterial.hpp>

using namespace gcm;


OrthotropicMaterial::
OrthotropicMaterial(const IsotropicMaterial& isotropic) {
	anglesOfRotation = Real3::Zeros();
	rho = isotropic.rho;
	yieldStrength = isotropic.yieldStrength;
	continualDamageParameter = isotropic.continualDamageParameter;
	tau0 = isotropic.tau0;
	c11 = c22 = c33 = isotropic.lambda + 2 * isotropic.mu;
	c44 = c55 = c66 = isotropic.mu;
	c12 = c13 = c23 = isotropic.lambda;
}


OrthotropicMaterial::
OrthotropicMaterial(const real rho_, std::initializer_list<real> list,
                    const real yieldStrength_, 
                    const real continualDamageParameter_, const Real3 phi,
                    const real tau0_) :
		rho(rho_), yieldStrength(yieldStrength_),
		continualDamageParameter(continualDamageParameter_),
		tau0(tau0_),
		anglesOfRotation(phi) {
	assert_eq(list.size(), 9);
	std::copy(list.begin(), list.end(), c);
}


OrthotropicMaterial OrthotropicMaterial::
generateRandomMaterial(const bool rotate) {
	const real RHO_MAX = 100.0;
	const real RHO_MIN = 0.01;
	const real C_MAX = 1e+3;
	const real C_MIN = 1.0;

	real rho = Utils::randomReal(RHO_MIN, RHO_MAX);
	Real3 phi = rotate ? linal::random<Real3>(-2*M_PI, 2*M_PI) : Real3::Zeros();
	
	auto help = linal::random<linal::Matrix33>(C_MIN, C_MAX);
	// symm - random symmetric positive-definite matrix
	auto symm = linal::transposeMultiply(help, help);

	real c11 = symm(0, 0); real c12 = symm(0, 1); real c13 = symm(0, 2);
	                       real c22 = symm(1, 1); real c23 = symm(1, 2);
	                                              real c33 = symm(2, 2);

	real c44 = Utils::randomReal(C_MIN*C_MIN, C_MAX*C_MAX);
	real c55 = Utils::randomReal(C_MIN*C_MIN, C_MAX*C_MAX);
	real c66 = Utils::randomReal(C_MIN*C_MIN, C_MAX*C_MAX);

	return OrthotropicMaterial(rho, 
			{c11, c12, c13, c22, c23, c33, c44, c55, c66}, 0, 0, phi);
}


