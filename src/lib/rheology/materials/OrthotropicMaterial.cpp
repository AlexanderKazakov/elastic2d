#include <lib/rheology/materials/OrthotropicMaterial.hpp>
#include <lib/rheology/materials/IsotropicMaterial.hpp>
#include <lib/util/task/Task.hpp>

using namespace gcm;


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

void OrthotropicMaterial::constructGcmMatrices(GcmMatrices<VelocitySigmaVariables<3>, OrthotropicMaterial> &m) const {
	m.m[0].A.initialize({0, 0, 0, -1.0 / rho, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
	                     -c11, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, -c66, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -c55, 0, 0, 0, 0, 0, 0,
	                     -c12, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, 0,
	                     -c13, 0, 0, 0, 0, 0, 0, 0, 0});

	m.m[0].L.initialize(
			{-sqrt(c66 / rho), sqrt(c66 / rho), -sqrt(c55 / rho), sqrt(c55 / rho), -sqrt(c11 / rho), sqrt(c11 / rho), 0,
			 0, 0});

	m.m[0].U.initialize({0, 1.0, 0, 0, 1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
	                     0, 1.0, 0, 0, -1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
	                     0, 0, 1.0, 0, 0, 1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
	                     0, 0, 1.0, 0, 0, -1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
	                     1.0, 0, 0, 1.0 / (sqrt(c11) * sqrt(rho)), 0, 0, 0, 0, 0,
	                     1.0, 0, 0, -1.0 / (sqrt(c11) * sqrt(rho)), 0, 0, 0, 0, 0,
	                     0, 0, 0, -c12 / c11, 0, 0, 1.0, 0, 0.0,
	                     0, 0, 0, 0, 0, 0, 0, 1.0, 0,
	                     0, 0, 0, -(1.0 * c13) / c11, 0, 0, 0, 0, 1.0});

	m.m[0].U1.initialize({0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
	                      0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
	                      0, 0, 0, 0, 0.5 * sqrt(c11) * sqrt(rho), -0.5 * sqrt(c11) * sqrt(rho), 0, 0, 0,
	                      0.5 * sqrt(c66) * sqrt(rho), -0.5 * sqrt(c66) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0.5 * sqrt(c55) * sqrt(rho), -0.5 * sqrt(c55) * sqrt(rho), 0, 0, 0, 0, 0,
	                      0, 0, 0, 0, (0.5 * c12 * sqrt(rho)) / sqrt(c11), -(0.5 * c12 * sqrt(rho)) / sqrt(c11), 1, 0, 0,
	                      0, 0, 0, 0, 0, 0, 0, 1, 0,
	                      0, 0, 0, 0, (0.5 * c13 * sqrt(rho)) / sqrt(c11), -(0.5 * c13 * sqrt(rho)) / sqrt(c11), 0, 0, 1});


	m.m[1].A.initialize({0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, -1.0 / rho, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
	                     0, -c12, 0, 0, 0, 0, 0, 0, 0,
	                     -c66, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, -c22, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -c44, 0, 0, 0, 0, 0, 0,
	                     0, -c23, 0, 0, 0, 0, 0, 0, 0});

	m.m[1].L.initialize(
			{-sqrt(c66 / rho), sqrt(c66 / rho), -sqrt(c44 / rho), sqrt(c44 / rho), -sqrt(c22 / rho), sqrt(c22 / rho), 0,
			 0, 0});

	m.m[1].U.initialize({1.0, 0, 0, 0, 1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
	                     1.0, 0, 0, 0, -1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
	                     0, 0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(c44) * sqrt(rho)), 0,
	                     0, 0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(c44) * sqrt(rho)), 0,
	                     0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(c22) * sqrt(rho)), 0, 0,
	                     0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(c22) * sqrt(rho)), 0, 0,
	                     0, 0, 0, 1.0, 0, 0, -c12 / c22, 0, 0.0,
	                     0, 0, 0, 0, 0, 1.0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, -(1.0 * c23) / c22, 0, 1.0});

	m.m[1].U1.initialize({0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
	                      0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
	                      0, 0, 0, 0, (0.5 * c12) / sqrt(c22 / rho), -(0.5 * c12) / sqrt(c22 / rho), 1, 0, 0,
	                      0.5 * sqrt(c66) * sqrt(rho), -0.5 * sqrt(c66) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0, 0, 0, 0, 0, 1, 0,
	                      0, 0, 0, 0, 0.5 * sqrt(c22) * sqrt(rho), -0.5 * sqrt(c22) * sqrt(rho), 0, 0, 0,
	                      0, 0, 0.5 * sqrt(c44) * sqrt(rho), -0.5 * sqrt(c44) * sqrt(rho), 0, 0, 0, 0, 0,
	                      0, 0, 0, 0, (0.5 * c23 * sqrt(rho)) / sqrt(c22), -(0.5 * c23 * sqrt(rho)) / sqrt(c22), 0, 0, 1});


	m.m[2].A.initialize({0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, -1.0 / rho,
	                     0, 0, -c13, 0, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, 0,
	                     -c55, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -c23, 0, 0, 0, 0, 0, 0,
	                     0, -c44, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -c33, 0, 0, 0, 0, 0, 0});

	m.m[2].L.initialize(
			{-sqrt(c55 / rho), sqrt(c55 / rho), -sqrt(c44 / rho), sqrt(c44 / rho), -sqrt(c33 / rho), sqrt(c33 / rho), 0,
			 0, 0});

	m.m[2].U.initialize({1.0, 0, 0, 0, 0, 1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
	                     1.0, 0, 0, 0, 0, -1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
	                     0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(c44) * sqrt(rho)), 0,
	                     0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(c44) * sqrt(rho)), 0,
	                     0, 0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(c33) * sqrt(rho)),
	                     0, 0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(c33) * sqrt(rho)),
	                     0, 0, 0, 1.0, 0, 0, 0, 0, -(1.0 * c13) / c33,
	                     0, 0, 0, 0, 1.0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 1.0, 0, -(1.0 * c23) / c33});

	m.m[2].U1.initialize({0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
	                      0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
	                      0, 0, 0, 0, (0.5 * c13 * sqrt(rho)) / sqrt(c33), -(0.5 * c13 * sqrt(rho)) / sqrt(c33), 1, 0, 0,
	                      0, 0, 0, 0, 0, 0, 0, 1, 0,
	                      0.5 * sqrt(c55) * sqrt(rho), -0.5 * sqrt(c55) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0, 0, (0.5 * c23 * sqrt(rho)) / sqrt(c33), -(0.5 * c23 * sqrt(rho)) / sqrt(c33), 0, 0, 1,
	                      0, 0, 0.5 * sqrt(c44) * sqrt(rho), -0.5 * sqrt(c44) * sqrt(rho), 0, 0, 0, 0, 0,
	                      0, 0, 0, 0, 0.5 * sqrt(c33) * sqrt(rho), -0.5 * sqrt(c33) * sqrt(rho), 0, 0, 0});
}


OrthotropicMaterial OrthotropicMaterial::generateRandomMaterial() {
	static const real RHO_MAX = 100.0;
	static const real RHO_MIN = 0.01;
	static const real C_MAX = 1e+3;
	static const real C_MIN = 1.0;

	real rho = ((RHO_MAX - RHO_MIN) * rand()) / RAND_MAX + RHO_MIN;
	// Generate symmetric positive-definite matrix
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
