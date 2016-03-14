#include <lib/rheology/materials/OrthotropicMaterial.hpp>
#include <lib/rheology/materials/IsotropicMaterial.hpp>
#include <lib/util/task/Task.hpp>

using namespace gcm;

const Materials::T OrthotropicMaterial::ID = Materials::T::ORTHOTROPIC;

OrthotropicMaterial::OrthotropicMaterial(const IsotropicMaterial &isotropic) {
	rho = isotropic.rho;
	yieldStrength = isotropic.yieldStrength;
	continualDamageParameter = isotropic.continualDamageParameter;
	c11 = c22 = c33 = isotropic.lambda + 2 * isotropic.mu;
	c44 = c55 = c66 = isotropic.mu;
	c12 = c13 = c23 = isotropic.lambda;
}

OrthotropicMaterial::OrthotropicMaterial(const real rho_, std::initializer_list<real> list,
                                         const real yieldStrength_, const real continualDamageParameter_) :
		rho(rho_), yieldStrength(yieldStrength_), continualDamageParameter(continualDamageParameter_) {
	int i = 0;
	for(auto& r : list) {
		c[i++] = r;
	}
}

void OrthotropicMaterial::initialize(const Statement& statement) {
	*this = statement.orthotropicMaterial;
}

OrthotropicMaterial OrthotropicMaterial::generateRandomMaterial() {
	const real RHO_MAX = 100.0;
	const real RHO_MIN = 0.01;
	const real C_MAX = 1e+3;
	const real C_MIN = 1.0;

	real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
	// Generate symmetric positive-definite matrix TODO - to linal
	linal::Matrix<3,3> help;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			help(i, j) = ((C_MAX - C_MIN) * rand()) / RAND_MAX + C_MIN;
		}
	}
	linal::Matrix<3,3> symm = help * help.transpose();

	real c11 = symm(0,0); real c12 = symm(0,1); real c13 = symm(0,2);
	                      real c22 = symm(1,1); real c23 = symm(1,2);
	                                            real c33 = symm(2,2);

	real c44 = ((C_MAX*C_MAX - C_MIN*C_MIN) * rand()) / RAND_MAX + C_MIN*C_MIN;
	real c55 = ((C_MAX*C_MAX - C_MIN*C_MIN) * rand()) / RAND_MAX + C_MIN*C_MIN;
	real c66 = ((C_MAX*C_MAX - C_MIN*C_MIN) * rand()) / RAND_MAX + C_MIN*C_MIN;

	return OrthotropicMaterial(rho, {c11, c12, c13, c22, c23, c33, c44, c55, c66});
}
