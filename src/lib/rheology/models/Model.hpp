#ifndef LIBGCM_MODEL_HPP
#define LIBGCM_MODEL_HPP

#include <lib/rheology/materials/materials.hpp>
#include <lib/rheology/variables/variables.hpp>
#include <lib/rheology/ode/ode.hpp>
#include <lib/rheology/correctors/correctors.hpp>
#include <lib/numeric/gcm/GcmMatrices.hpp>


namespace gcm {

/**
 * Map between wave types and corresponding columns in matrix U1 in GcmMatrix.
 */
typedef std::map<Waves::T, int>                      WavesEigenvectorsMap;
/**
 * Map between material types and corresponding WavesEigenvectorsMap for concrete Model.
 * Note that in concrete Model for concrete Material the order of eigenvalues
 * in GcmMatrices has to be the same for all spatial directions.
 */
typedef std::map<Materials::T, WavesEigenvectorsMap> MaterialsWavesMap;

/**
 * The list of rheology models
 */

class Elastic1DModel {
public:
	static const int DIMENSIONALITY = 1;

	typedef VelocitySigmaVariables<DIMENSIONALITY> PdeVariables;
	typedef DummyOde                               InternalOde;
	typedef DummyCorrector                         Corrector;
	typedef typename PdeVariables::PdeVector       PdeVector;

	static const int PDE_SIZE = PdeVector::SIZE;

	typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY>  GCM_MATRICES;
	typedef typename InternalOde::Variables        OdeVariables;
	typedef std::shared_ptr<GCM_MATRICES>          GcmMatricesPtr;
	typedef std::shared_ptr<const GCM_MATRICES>    ConstGcmMatricesPtr;

	static const MaterialsWavesMap MATERIALS_WAVES_MAP;

	template<typename ... Args>
	static void constructGcmMatrices(GcmMatricesPtr m,
	                                 std::shared_ptr<const IsotropicMaterial> material,
	                                 const Args& ...);

};


class Elastic2DModel {
public:
	static const int DIMENSIONALITY = 2;

	typedef VelocitySigmaVariables<DIMENSIONALITY> PdeVariables;
	typedef DummyOde                               InternalOde;
	typedef DummyCorrector                         Corrector;
	typedef typename PdeVariables::PdeVector       PdeVector;

	static const int PDE_SIZE = PdeVector::SIZE;
	
	typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY>  GCM_MATRICES;
	typedef typename InternalOde::Variables        OdeVariables;
	typedef std::shared_ptr<GCM_MATRICES>          GcmMatricesPtr;
	typedef std::shared_ptr<const GCM_MATRICES>    ConstGcmMatricesPtr;

	static const MaterialsWavesMap MATERIALS_WAVES_MAP;

	template<typename ... Args>
	static void constructGcmMatrices(GcmMatricesPtr m,
	                                 std::shared_ptr<const IsotropicMaterial> material,
	                                 const Args& ...);

	template<typename ... Args>
	static void constructGcmMatrices(GcmMatricesPtr m,
	                                 std::shared_ptr<const OrthotropicMaterial> material,
	                                 const Args& ...);

};


class Elastic3DModel {
public:
	static const int DIMENSIONALITY = 3;

	typedef VelocitySigmaVariables<DIMENSIONALITY> PdeVariables;
	typedef DummyOde                               InternalOde;
	typedef DummyCorrector                         Corrector;
	typedef typename PdeVariables::PdeVector       PdeVector;

	static const int PDE_SIZE = PdeVector::SIZE;
	
	typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY>  GCM_MATRICES;
	typedef typename InternalOde::Variables        OdeVariables;
	typedef std::shared_ptr<GCM_MATRICES>          GcmMatricesPtr;
	typedef std::shared_ptr<const GCM_MATRICES>    ConstGcmMatricesPtr;

	static const MaterialsWavesMap MATERIALS_WAVES_MAP;

	template<typename ... Args>
	static void constructGcmMatrices(GcmMatricesPtr m,
	                                 std::shared_ptr<const IsotropicMaterial> material,
	                                 const Args& ...);

	template<typename ... Args>
	static void constructGcmMatrices(GcmMatricesPtr m,
	                                 std::shared_ptr<const OrthotropicMaterial> material,
	                                 const Args& ...);

};


class SuperDuperModel {
public:
	static const int DIMENSIONALITY = 3;

	typedef VelocitySigmaVariables<DIMENSIONALITY> PdeVariables;
	typedef ContinualDamageOde                     InternalOde;
	typedef IdealPlasticFlowCorrector              Corrector;
	typedef typename PdeVariables::PdeVector       PdeVector;

	static const int PDE_SIZE = PdeVector::SIZE;
	 
	typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY>  GCM_MATRICES;
	typedef typename InternalOde::Variables        OdeVariables;
	typedef std::shared_ptr<GCM_MATRICES>          GcmMatricesPtr;
	typedef std::shared_ptr<const GCM_MATRICES>    ConstGcmMatricesPtr;

	static const MaterialsWavesMap MATERIALS_WAVES_MAP;

	template<typename TMaterialPtr, typename ... Args>
	static void constructGcmMatrices(GcmMatricesPtr m, const TMaterialPtr material,
	                                 const Args& ... args) {
		Elastic3DModel::constructGcmMatrices(m, material, args ...);
	}

};


template<typename ... Args>
void Elastic1DModel::
constructGcmMatrices(GcmMatricesPtr m, std::shared_ptr<const IsotropicMaterial> material,
                     const Args& ...) {
	const real rho = material->rho;
	const real lambda = material->lambda;
	const real mu = material->mu;
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


template<typename ... Args>
void Elastic2DModel::
constructGcmMatrices(GcmMatricesPtr m, std::shared_ptr<const IsotropicMaterial> material,
                     const Args& ...) {
	const real rho = material->rho;
	const real lambda = material->lambda;
	const real mu = material->mu;

	// TODO - actually we can use orthotropic material here
	m->m[0].A.initialize({0, 0, -1.0 / rho, 0, 0,
	                      0, 0, 0, -1.0 / rho, 0,
	                      -lambda - 2.0 * mu, 0, 0, 0, 0,
	                      0, -mu, 0, 0, 0,
	                      -lambda, 0, 0, 0, 0});

	m->m[0].L.initialize(
	        {-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(
	                 mu / rho), sqrt(mu / rho), 0});

	m->m[0].U.initialize({1.0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                      1.0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
	                      0, 1.0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                      0, 1.0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                      0, 0, 1.0 / (lambda + 2 * mu), 0, -1.0 / lambda});

	m->m[0].U1.initialize({0.5, 0.5, 0, 0, 0,
	                       0, 0, 0.5, 0.5, 0,
	                       0.5 *
	                       sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(
	                               rho * (lambda + 2 * mu)), 0, 0, 0,
	                       0, 0, 0.5 * sqrt(mu * rho), -0.5 * sqrt(mu * rho), 0,
	                       (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                       -(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, -lambda});


	m->m[1].A.initialize({0, 0, 0, -1.0 / rho, 0,
	                      0, 0, 0, 0, -1.0 / rho,
	                      0, -lambda, 0, 0, 0,
	                      -mu, 0, 0, 0, 0,
	                      0, -lambda - 2.0 * mu, 0, 0, 0});

	m->m[1].L.initialize(
	        {-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(
	                 mu / rho), sqrt(mu / rho), 0});

	m->m[1].U.initialize({0, 1.0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                      0, 1.0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))),
	                      1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
	                      1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
	                      0, 0, 1.0, 0, -(1.0 * lambda) / (lambda + 2 * mu)});

	m->m[1].U1.initialize({0, 0, 0.5, 0.5, 0,
	                       0.5, 0.5, 0, 0, 0,
	                       (0.5 * lambda) / sqrt(
	                               (lambda + 2 * mu) / rho), -(0.5 * lambda) / sqrt(
	                               (lambda + 2 * mu) / rho), 0, 0, 1.0,
	                       0, 0, 0.5 * sqrt(mu * rho), -0.5 * sqrt(mu * rho), 0,
	                       0.5 *
	                       sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(
	                               rho * (lambda + 2 * mu)), 0, 0, 0});

	m->checkDecomposition();
}


template<typename ... Args>
void Elastic2DModel::
constructGcmMatrices(GcmMatricesPtr /*m*/, std::shared_ptr<const OrthotropicMaterial> /*material*/,
                     const Args& ...) {
	THROW_UNSUPPORTED("TODO");
}


template<typename ... Args>
void Elastic3DModel::
constructGcmMatrices(GcmMatricesPtr m, std::shared_ptr<const IsotropicMaterial> material,
                     const Args& ...) {
	const real rho = material->rho;
	const real lambda = material->lambda;
	const real mu = material->mu;

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
	        {-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(
	                 mu / rho), -sqrt(mu / rho),
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
	                       0.5 *
	                       sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(
	                               rho * (lambda + 2 * mu)), 0, 0, 0, 0, 0, 0, 0,
	                       0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
	                       0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
	                       (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                       -(0.5 * sqrt(rho) * lambda) / sqrt(
	                               lambda + 2 * mu), 0, 0, 0, 0, 1, 0, 0,
	                       0, 0, 0, 0, 0, 0, 0, 1, 0,
	                       (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                       -(0.5 * sqrt(rho) * lambda) / sqrt(
	                               lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1});


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
	        {-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(
	                 mu / rho), -sqrt(mu / rho),
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
	                       (0.5 * lambda) / sqrt(
	                               (lambda + 2 * mu) / rho), -(0.5 * lambda) / sqrt(
	                               (lambda + 2 * mu) / rho), 0, 0, 0, 0, 1, 0, 0,
	                       0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
	                       0, 0, 0, 0, 0, 0, 0, 1, 0,
	                       0.5 *
	                       sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(
	                               rho * (lambda + 2 * mu)), 0, 0, 0, 0, 0, 0, 0,
	                       0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
	                       (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                       -(0.5 * sqrt(rho) * lambda) / sqrt(
	                               lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1});


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
	        {-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(
	                 mu / rho), -sqrt(mu / rho),
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
	                       (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                       -(0.5 * sqrt(rho) * lambda) / sqrt(
	                               lambda + 2 * mu), 0, 0, 0, 0, 1, 0, 0,
	                       0, 0, 0, 0, 0, 0, 0, 1, 0,
	                       0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
	                       (0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
	                       -(0.5 * sqrt(rho) * lambda) / sqrt(
	                               lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1,
	                       0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
	                       0.5 *
	                       sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(
	                               rho * (lambda + 2 * mu)), 0, 0, 0, 0, 0, 0, 0});

	m->checkDecomposition();
}


template<typename ... Args>
void Elastic3DModel::
constructGcmMatrices(GcmMatricesPtr m, std::shared_ptr<const OrthotropicMaterial> material,
                     const Args& ...) {
#define COPY_FROM_MATERIAL(var) const real var = material->var
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
	        {-sqrt(c66 / rho), sqrt(c66 / rho), -sqrt(c55 / rho), sqrt(c55 / rho),
	         -sqrt(c11 / rho), sqrt(c11 / rho), 0,
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
	                       0, 0, 0, 0, 0.5 * sqrt(c11) * sqrt(rho), -0.5 * sqrt(c11) * sqrt(
	                               rho), 0, 0, 0,
	                       0.5 * sqrt(c66) * sqrt(rho), -0.5 * sqrt(c66) * sqrt(
	                               rho), 0, 0, 0, 0, 0, 0, 0,
	                       0, 0, 0.5 * sqrt(c55) * sqrt(rho), -0.5 * sqrt(c55) * sqrt(
	                               rho), 0, 0, 0, 0, 0,
	                       0, 0, 0, 0, (0.5 * c12 * sqrt(rho)) / sqrt(c11),
	                       -(0.5 * c12 * sqrt(rho)) / sqrt(c11), 1, 0, 0,
	                       0, 0, 0, 0, 0, 0, 0, 1, 0,
	                       0, 0, 0, 0, (0.5 * c13 * sqrt(rho)) / sqrt(c11),
	                       -(0.5 * c13 * sqrt(rho)) / sqrt(c11), 0, 0, 1});


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
	        {-sqrt(c66 / rho), sqrt(c66 / rho), -sqrt(c44 / rho), sqrt(c44 / rho),
	         -sqrt(c22 / rho), sqrt(c22 / rho), 0,
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
	                       0, 0, 0, 0, (0.5 * c12) / sqrt(c22 / rho),
	                       -(0.5 * c12) / sqrt(c22 / rho), 1, 0, 0,
	                       0.5 * sqrt(c66) * sqrt(rho), -0.5 * sqrt(c66) * sqrt(
	                               rho), 0, 0, 0, 0, 0, 0, 0,
	                       0, 0, 0, 0, 0, 0, 0, 1, 0,
	                       0, 0, 0, 0, 0.5 * sqrt(c22) * sqrt(rho), -0.5 * sqrt(c22) * sqrt(
	                               rho), 0, 0, 0,
	                       0, 0, 0.5 * sqrt(c44) * sqrt(rho), -0.5 * sqrt(c44) * sqrt(
	                               rho), 0, 0, 0, 0, 0,
	                       0, 0, 0, 0, (0.5 * c23 * sqrt(rho)) / sqrt(c22),
	                       -(0.5 * c23 * sqrt(rho)) / sqrt(c22), 0, 0, 1});


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
	        {-sqrt(c55 / rho), sqrt(c55 / rho), -sqrt(c44 / rho), sqrt(c44 / rho),
	         -sqrt(c33 / rho), sqrt(c33 / rho), 0,
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
	                       0, 0, 0, 0, (0.5 * c13 * sqrt(rho)) / sqrt(c33),
	                       -(0.5 * c13 * sqrt(rho)) / sqrt(c33), 1, 0, 0,
	                       0, 0, 0, 0, 0, 0, 0, 1, 0,
	                       0.5 * sqrt(c55) * sqrt(rho), -0.5 * sqrt(c55) * sqrt(
	                               rho), 0, 0, 0, 0, 0, 0, 0,
	                       0, 0, 0, 0, (0.5 * c23 * sqrt(rho)) / sqrt(c33),
	                       -(0.5 * c23 * sqrt(rho)) / sqrt(c33), 0, 0, 1,
	                       0, 0, 0.5 * sqrt(c44) * sqrt(rho), -0.5 * sqrt(c44) * sqrt(
	                               rho), 0, 0, 0, 0, 0,
	                       0, 0, 0, 0, 0.5 * sqrt(c33) * sqrt(rho), -0.5 * sqrt(c33) * sqrt(
	                               rho), 0, 0, 0});

	m->checkDecomposition();
}


}

#endif // LIBGCM_MODEL_HPP
