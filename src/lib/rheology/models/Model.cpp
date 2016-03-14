#include <lib/rheology/models/Model.hpp>

using namespace gcm;

const MaterialsWavesMap Elastic1DModel::MATERIALS_WAVES_MAP = {
		{IsotropicMaterial::ID, {
				{Waves::T::P_FORWARD,  0},
				{Waves::T::P_BACKWARD, 1}
		 }}
};

const MaterialsWavesMap Elastic2DModel::MATERIALS_WAVES_MAP = {
		{IsotropicMaterial::ID,	{
				{Waves::T::P_FORWARD,   1},
				{Waves::T::P_BACKWARD,  0},
				{Waves::T::S1_FORWARD,  3},
				{Waves::T::S1_BACKWARD, 2}
		 }}
};

const MaterialsWavesMap Elastic3DModel::MATERIALS_WAVES_MAP = {
		{IsotropicMaterial::ID, {
				{Waves::T::P_FORWARD,   1},
				{Waves::T::P_BACKWARD,  0},
				{Waves::T::S1_FORWARD,  4},
				{Waves::T::S1_BACKWARD, 2},
				{Waves::T::S2_FORWARD,  5},
				{Waves::T::S2_BACKWARD, 3}
		 }},
		{OrthotropicMaterial::ID, {
				 {Waves::T::P_FORWARD,   5},
				 {Waves::T::P_BACKWARD,  4},
				 {Waves::T::S1_FORWARD,  1},
				 {Waves::T::S1_BACKWARD, 0},
				 {Waves::T::S2_FORWARD,  3},
				 {Waves::T::S2_BACKWARD, 2}
		 }}
};

const MaterialsWavesMap SuperDuperModel::MATERIALS_WAVES_MAP =
		Elastic3DModel::MATERIALS_WAVES_MAP;


void Elastic1DModel::constructGcmMatrices
(GcmMatricesPtr m, const PdeVector&, const IsotropicMaterial& material) {
	const real rho = material.rho;
	const real lambda = material.lambda;
	const real mu = material.mu;
	const real E = mu * (3 * lambda + 2 * mu) / (lambda + mu); // Young's modulus

	m->m[0].A.initialize({0.0, -1.0 / rho,
	                     -E, 0.0});

	m->m[0].L.initialize({sqrt(E / rho), -sqrt(E / rho)});

	m->m[0].U.initialize({-0.5, 1.0 / (2 * sqrt(E * rho)),
	                     0.5, 1.0 / (2 * sqrt(E * rho))});

	m->m[0].U1.initialize({-1.0, 1.0,
	                      sqrt(E * rho), sqrt(E * rho)});

	m->checkDecomposition();
}

void Elastic2DModel::constructGcmMatrices
(GcmMatricesPtr m, const PdeVector&, const IsotropicMaterial& material) {
	const real rho = material.rho;
	const real lambda = material.lambda;
	const real mu = material.mu;

	// TODO - actually we can use orthotropic material here
	m->m[0].A.initialize({0, 0, -1.0 / rho, 0, 0,
	                     0, 0, 0, -1.0 / rho, 0,
	                     -lambda - 2.0 * mu, 0, 0, 0, 0,
	                     0, -mu, 0, 0, 0,
	                     -lambda, 0, 0, 0, 0});

	m->m[0].L.initialize(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), sqrt(mu / rho), 0});

	m->m[0].U.initialize({1.0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                     1.0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                     0, 1.0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                     0, 1.0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                     0, 0, 1.0 / (lambda + 2 * mu), 0, -1.0 / lambda});

	m->m[0].U1.initialize({0.5, 0.5, 0, 0, 0,
	                      0, 0, 0.5, 0.5, 0,
	                      0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0,
	                      0, 0, 0.5 * sqrt(mu * rho), -0.5 * sqrt(mu * rho), 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, -lambda});


	m->m[1].A.initialize({0, 0, 0, -1.0 / rho, 0,
	                     0, 0, 0, 0, -1.0 / rho,
	                     0, -lambda, 0, 0, 0,
	                     -mu, 0, 0, 0, 0,
	                     0, -lambda - 2.0 * mu, 0, 0, 0});

	m->m[1].L.initialize(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), sqrt(mu / rho), 0});

	m->m[1].U.initialize({0, 1.0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                     0, 1.0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                     1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                     1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                     0, 0, 1.0, 0, -(1.0 * lambda) / (lambda + 2 * mu)});

	m->m[1].U1.initialize({0, 0, 0.5, 0.5, 0,
	                      0.5, 0.5, 0, 0, 0,
	                      (0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), -(0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), 0, 0, 1.0,
	                      0, 0, 0.5 * sqrt(mu * rho), -0.5 * sqrt(mu * rho), 0,
	                      0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0});

	m->checkDecomposition();
}

void Elastic3DModel::constructGcmMatrices
(GcmMatricesPtr m, const PdeVector&, const IsotropicMaterial& material) {
	const real rho = material.rho;
	const real lambda = material.lambda;
	const real mu = material.mu;

	// TODO - actually we can use orthotropic material here
	m->m[0].A.initialize({0, 0, 0, -1.0 / rho, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
	                     -lambda - 2 * mu, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, -mu, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -mu, 0, 0, 0, 0, 0, 0,
	                     -lambda, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, 0,
	                     -lambda, 0, 0, 0, 0, 0, 0, 0, 0});

	m->m[0].L.initialize(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), -sqrt(mu / rho),
			 sqrt(mu / rho), sqrt(mu / rho), 0, 0, 0});

	m->m[0].U.initialize({1.0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0, 0, 0, 0,
	                     1.0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0, 0, 0, 0,
	                     0, 1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
	                     0, 0, 1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0,
	                     0, 1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
	                     0, 0, 1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0,
	                     0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 0, 1.0, 0, 0.0,
	                     0, 0, 0, 0, 0, 0, 0, 1.0, 0,
	                     0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 0, 0, 0, 1.0});

	m->m[0].U1.initialize({0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0.5, 0, 0.5, 0, 0, 0, 0,
	                      0, 0, 0, 0.5, 0, 0.5, 0, 0, 0,
	                      0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
	                      0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 1, 0, 0,
	                      0, 0, 0, 0, 0, 0, 0, 1, 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1});


	m->m[1].A.initialize({0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, -1.0 / rho, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
	                     0, -lambda, 0, 0, 0, 0, 0, 0, 0,
	                     -mu, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, -lambda - 2 * mu, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -mu, 0, 0, 0, 0, 0, 0,
	                     0, -lambda, 0, 0, 0, 0, 0, 0, 0});

	m->m[1].L.initialize(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), -sqrt(mu / rho),
			 sqrt(mu / rho), sqrt(mu / rho), 0, 0, 0});

	m->m[1].U.initialize({0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                     0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                     1.0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
	                     0, 0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                     1.0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
	                     0, 0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                     0, 0, 0, 1.0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 0.0,
	                     0, 0, 0, 0, 0, 1.0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 1.0});

	m->m[1].U1.initialize({0, 0, 0.5, 0, 0.5, 0, 0, 0, 0,
	                      0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0, 0.5, 0, 0.5, 0, 0, 0,
	                      (0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), -(0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), 0, 0, 0, 0, 1, 0, 0,
	                      0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
	                      0, 0, 0, 0, 0, 0, 0, 1, 0,
	                      0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0, 0, 0, 0, 0,
	                      0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1});


	m->m[2].A.initialize({0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, -1.0 / rho,
	                     0, 0, -lambda, 0, 0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 0, 0, 0,
	                     -mu, 0, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -lambda, 0, 0, 0, 0, 0, 0,
	                     0, -mu, 0, 0, 0, 0, 0, 0, 0,
	                     0, 0, -lambda - 2 * mu, 0, 0, 0, 0, 0, 0});

	m->m[2].L.initialize(
			{-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(mu / rho), -sqrt(mu / rho),
			 sqrt(mu / rho), sqrt(mu / rho), 0, 0, 0});

	m->m[2].U.initialize({0, 0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                     0, 0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                     1.0, 0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0,
	                     0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                     1.0, 0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0,
	                     0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                     0, 0, 0, 1.0, 0, 0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu),
	                     0, 0, 0, 0, 1.0, 0, 0, 0, 0,
	                     0, 0, 0, 0, 0, 0, 1.0, 0, -(1.0 * lambda) / (lambda + 2 * mu)});

	m->m[2].U1.initialize({0, 0, 0.5, 0, 0.5, 0, 0, 0, 0,
	                      0, 0, 0, 0.5, 0, 0.5, 0, 0, 0,
	                      0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 1, 0, 0,
	                      0, 0, 0, 0, 0, 0, 0, 1, 0,
	                      0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
	                      (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1,
	                      0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
	                      0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0, 0, 0, 0, 0});

	m->checkDecomposition();
}


void Elastic3DModel::constructGcmMatrices
(GcmMatricesPtr m, const PdeVector&, const OrthotropicMaterial& material) {
#define COPY_FROM_MATERIAL(var) const real var = material.var
	COPY_FROM_MATERIAL(rho);
	COPY_FROM_MATERIAL(c11);
	COPY_FROM_MATERIAL(c12);
	COPY_FROM_MATERIAL(c13);
	COPY_FROM_MATERIAL(c22);
	COPY_FROM_MATERIAL(c23);
	COPY_FROM_MATERIAL(c33);
	COPY_FROM_MATERIAL(c44);
	COPY_FROM_MATERIAL(c55);
	COPY_FROM_MATERIAL(c66);
#undef COPY_FROM_MATERIAL

	m->m[0].A.initialize({0, 0, 0, -1.0 / rho, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
						 -c11, 0, 0, 0, 0, 0, 0, 0, 0,
						 0, -c66, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, -c55, 0, 0, 0, 0, 0, 0,
						 -c12, 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0, 0,
						 -c13, 0, 0, 0, 0, 0, 0, 0, 0});

	m->m[0].L.initialize(
			{-sqrt(c66 / rho), sqrt(c66 / rho), -sqrt(c55 / rho), sqrt(c55 / rho), -sqrt(c11 / rho), sqrt(c11 / rho), 0,
			 0, 0});

	m->m[0].U.initialize({0, 1.0, 0, 0, 1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
						 0, 1.0, 0, 0, -1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
						 0, 0, 1.0, 0, 0, 1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
						 0, 0, 1.0, 0, 0, -1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
						 1.0, 0, 0, 1.0 / (sqrt(c11) * sqrt(rho)), 0, 0, 0, 0, 0,
						 1.0, 0, 0, -1.0 / (sqrt(c11) * sqrt(rho)), 0, 0, 0, 0, 0,
						 0, 0, 0, -c12 / c11, 0, 0, 1.0, 0, 0.0,
						 0, 0, 0, 0, 0, 0, 0, 1.0, 0,
						 0, 0, 0, -(1.0 * c13) / c11, 0, 0, 0, 0, 1.0});

	m->m[0].U1.initialize({0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
						  0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0.5 * sqrt(c11) * sqrt(rho), -0.5 * sqrt(c11) * sqrt(rho), 0, 0, 0,
						  0.5 * sqrt(c66) * sqrt(rho), -0.5 * sqrt(c66) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0.5 * sqrt(c55) * sqrt(rho), -0.5 * sqrt(c55) * sqrt(rho), 0, 0, 0, 0, 0,
						  0, 0, 0, 0, (0.5 * c12 * sqrt(rho)) / sqrt(c11), -(0.5 * c12 * sqrt(rho)) / sqrt(c11), 1, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 1, 0,
						  0, 0, 0, 0, (0.5 * c13 * sqrt(rho)) / sqrt(c11), -(0.5 * c13 * sqrt(rho)) / sqrt(c11), 0, 0, 1});


	m->m[1].A.initialize({0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, -1.0 / rho, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
						 0, -c12, 0, 0, 0, 0, 0, 0, 0,
						 -c66, 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0, 0,
						 0, -c22, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, -c44, 0, 0, 0, 0, 0, 0,
						 0, -c23, 0, 0, 0, 0, 0, 0, 0});

	m->m[1].L.initialize(
			{-sqrt(c66 / rho), sqrt(c66 / rho), -sqrt(c44 / rho), sqrt(c44 / rho), -sqrt(c22 / rho), sqrt(c22 / rho), 0,
			 0, 0});

	m->m[1].U.initialize({1.0, 0, 0, 0, 1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
						 1.0, 0, 0, 0, -1.0 / (sqrt(c66) * sqrt(rho)), 0, 0, 0, 0,
						 0, 0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(c44) * sqrt(rho)), 0,
						 0, 0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(c44) * sqrt(rho)), 0,
						 0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(c22) * sqrt(rho)), 0, 0,
						 0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(c22) * sqrt(rho)), 0, 0,
						 0, 0, 0, 1.0, 0, 0, -c12 / c22, 0, 0.0,
						 0, 0, 0, 0, 0, 1.0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, -(1.0 * c23) / c22, 0, 1.0});

	m->m[1].U1.initialize({0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
						  0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, (0.5 * c12) / sqrt(c22 / rho), -(0.5 * c12) / sqrt(c22 / rho), 1, 0, 0,
						  0.5 * sqrt(c66) * sqrt(rho), -0.5 * sqrt(c66) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 1, 0,
						  0, 0, 0, 0, 0.5 * sqrt(c22) * sqrt(rho), -0.5 * sqrt(c22) * sqrt(rho), 0, 0, 0,
						  0, 0, 0.5 * sqrt(c44) * sqrt(rho), -0.5 * sqrt(c44) * sqrt(rho), 0, 0, 0, 0, 0,
						  0, 0, 0, 0, (0.5 * c23 * sqrt(rho)) / sqrt(c22), -(0.5 * c23 * sqrt(rho)) / sqrt(c22), 0, 0, 1});


	m->m[2].A.initialize({0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
						 0, 0, 0, 0, 0, 0, 0, 0, -1.0 / rho,
						 0, 0, -c13, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0, 0,
						 -c55, 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, -c23, 0, 0, 0, 0, 0, 0,
						 0, -c44, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, -c33, 0, 0, 0, 0, 0, 0});

	m->m[2].L.initialize(
			{-sqrt(c55 / rho), sqrt(c55 / rho), -sqrt(c44 / rho), sqrt(c44 / rho), -sqrt(c33 / rho), sqrt(c33 / rho), 0,
			 0, 0});

	m->m[2].U.initialize({1.0, 0, 0, 0, 0, 1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
						 1.0, 0, 0, 0, 0, -1.0 / (sqrt(c55) * sqrt(rho)), 0, 0, 0,
						 0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(c44) * sqrt(rho)), 0,
						 0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(c44) * sqrt(rho)), 0,
						 0, 0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(c33) * sqrt(rho)),
						 0, 0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(c33) * sqrt(rho)),
						 0, 0, 0, 1.0, 0, 0, 0, 0, -(1.0 * c13) / c33,
						 0, 0, 0, 0, 1.0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 1.0, 0, -(1.0 * c23) / c33});

	m->m[2].U1.initialize({0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
						  0, 0, 0, 0, (0.5 * c13 * sqrt(rho)) / sqrt(c33), -(0.5 * c13 * sqrt(rho)) / sqrt(c33), 1, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 1, 0,
						  0.5 * sqrt(c55) * sqrt(rho), -0.5 * sqrt(c55) * sqrt(rho), 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, (0.5 * c23 * sqrt(rho)) / sqrt(c33), -(0.5 * c23 * sqrt(rho)) / sqrt(c33), 0, 0, 1,
						  0, 0, 0.5 * sqrt(c44) * sqrt(rho), -0.5 * sqrt(c44) * sqrt(rho), 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0.5 * sqrt(c33) * sqrt(rho), -0.5 * sqrt(c33) * sqrt(rho), 0, 0, 0});

	m->checkDecomposition();
}

