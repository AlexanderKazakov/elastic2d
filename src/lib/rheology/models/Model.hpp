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

	static const int PDE_SIZE = PdeVector::M;

	typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY> GCM_MATRICES;
	typedef typename InternalOde::Variables       OdeVariables;
	typedef std::shared_ptr<GCM_MATRICES>         GcmMatricesPtr;
	typedef std::shared_ptr<const GCM_MATRICES>   ConstGcmMatricesPtr;

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

	static const int PDE_SIZE = PdeVector::M;

	typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY> GCM_MATRICES;
	typedef typename InternalOde::Variables       OdeVariables;
	typedef std::shared_ptr<GCM_MATRICES>         GcmMatricesPtr;
	typedef std::shared_ptr<const GCM_MATRICES>   ConstGcmMatricesPtr;

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

	static const int PDE_SIZE = PdeVector::M;

	typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY> GCM_MATRICES;
	typedef typename InternalOde::Variables       OdeVariables;
	typedef std::shared_ptr<GCM_MATRICES>         GcmMatricesPtr;
	typedef std::shared_ptr<const GCM_MATRICES>   ConstGcmMatricesPtr;

	static const MaterialsWavesMap MATERIALS_WAVES_MAP;

	template<typename ... Args>
	static void constructGcmMatrices(GcmMatricesPtr m,
	                                 std::shared_ptr<const IsotropicMaterial> material,
	                                 const Args& ...);

	template<typename ... Args>
	static void constructGcmMatrices(GcmMatricesPtr m,
	                                 std::shared_ptr<const OrthotropicMaterial> material,
	                                 const Args& ...);

private:
	template<typename ... Args>
	static void constructRotated(GcmMatricesPtr m,
	                             std::shared_ptr<const OrthotropicMaterial> material,
	                             const Args& ...);
	
	template<typename ... Args>
	static void constructNotRotated(GcmMatricesPtr m,
	                                std::shared_ptr<const OrthotropicMaterial> material,
	                                const Args& ...);

	static linal::VECTOR<3, long double> constructEigenvaluesPolynomial(
			ConstGcmMatricesPtr m, const int s);
	
	static std::vector<linal::VECTOR<9, long double>> findEigenvectors(
		const long double l, const linal::Matrix<9, 9>& A,
		const int stage, const int numberOfVectorsToSearch);
	
	static std::vector<linal::VECTOR<9, long double>> findEigenstrings(
		const long double l, const linal::Matrix<9, 9>& A,
		const int stage, const int numberOfStringsToSearch);
	
	static void getColumnsWithRho(const int stage, int& i, int& j, int& k) {
	/// indices of columns which contains (- 1 / rho) in increasing order
		switch (stage) {
		case 0 :
			i = 3; j = 4; k = 5;
			break;
		case 1 :
			i = 4; j = 6; k = 7;
			break;
		case 2 :
			i = 5; j = 7; k = 8;
			break;
		default:
			THROW_INVALID_ARG("Invalid stage number");
		}
	}
	
	static void getZeroColumns(const int stage, int& i, int& j, int& k) {
	/// indices of columns with all zeros in increasing order
		switch (stage) {
		case 0 :
			i = 6; j = 7; k = 8;
			break;
		case 1 :
			i = 3; j = 5; k = 8;
			break;
		case 2 :
			i = 3; j = 4; k = 6;
			break;
		default:
			THROW_INVALID_ARG("Invalid stage number");
		}
	}
};


class SuperDuperModel {
public:
	static const int DIMENSIONALITY = 3;

	typedef VelocitySigmaVariables<DIMENSIONALITY> PdeVariables;
	typedef ContinualDamageOde                     InternalOde;
	typedef IdealPlasticFlowCorrector              Corrector;
	typedef typename PdeVariables::PdeVector       PdeVector;

	static const int PDE_SIZE = PdeVector::M;

	typedef GcmMatrices<PDE_SIZE, DIMENSIONALITY> GCM_MATRICES;
	typedef typename InternalOde::Variables       OdeVariables;
	typedef std::shared_ptr<GCM_MATRICES>         GcmMatricesPtr;
	typedef std::shared_ptr<const GCM_MATRICES>   ConstGcmMatricesPtr;

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

	m->m[0].A = {
		0.0, -1.0 / rho,
		-E, 0.0
	};

	m->m[0].L = {sqrt(E / rho), -sqrt(E / rho)};

	m->m[0].U = {
		-0.5, 1.0 / (2 * sqrt(E * rho)),
		0.5, 1.0 / (2 * sqrt(E * rho))
	};

	m->m[0].U1 = {
		-1.0, 1.0,
		sqrt(E * rho), sqrt(E * rho)
	};

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
	m->m[0].A = {
		0, 0, -1.0 / rho, 0, 0,
		0, 0, 0, -1.0 / rho, 0,
		-lambda - 2.0 * mu, 0, 0, 0, 0,
		0, -mu, 0, 0, 0,
		-lambda, 0, 0, 0, 0
	};

	m->m[0].L = {
		-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho),
		-sqrt(mu / rho), sqrt(mu / rho), 0
	};

	m->m[0].U = {
		1.0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
		1.0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
		0, 1.0, 0, 1.0 / (sqrt(mu * rho)), 0,
		0, 1.0, 0, -1.0 / (sqrt(mu * rho)), 0,
		0, 0, 1.0 / (lambda + 2 * mu), 0, -1.0 / lambda
	};

	m->m[0].U1 = {
		0.5, 0.5, 0, 0, 0,
		0, 0, 0.5, 0.5, 0,
		0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0,
		0, 0, 0.5 * sqrt(mu * rho), -0.5 * sqrt(mu * rho), 0,
		(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
		-(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, -lambda
	};


	m->m[1].A = {
		0, 0, 0, -1.0 / rho, 0,
		0, 0, 0, 0, -1.0 / rho,
		0, -lambda, 0, 0, 0,
		-mu, 0, 0, 0, 0,
		0, -lambda - 2.0 * mu, 0, 0, 0
	};

	m->m[1].L =
	{
		-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(
		        mu / rho), sqrt(mu / rho), 0
	};

	m->m[1].U = {
		0, 1.0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))),
		0, 1.0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))),
		1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
		1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
		0, 0, 1.0, 0, -(1.0 * lambda) / (lambda + 2 * mu)
	};

	m->m[1].U1 = {
		0, 0, 0.5, 0.5, 0,
		0.5, 0.5, 0, 0, 0,
		(0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), -(0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), 0, 0, 1.0,
		0, 0, 0.5 * sqrt(mu * rho), -0.5 * sqrt(mu * rho), 0,
		0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0
	};

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
	m->m[0].A = {
		0, 0, 0, -1.0 / rho, 0, 0, 0, 0, 0,
		0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
		0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
		-lambda - 2 * mu, 0, 0, 0, 0, 0, 0, 0, 0,
		0, -mu, 0, 0, 0, 0, 0, 0, 0,
		0, 0, -mu, 0, 0, 0, 0, 0, 0,
		-lambda, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0,
		-lambda, 0, 0, 0, 0, 0, 0, 0, 0
	};

	m->m[0].L =
	{
		-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(
		        mu / rho), -sqrt(mu / rho),
		sqrt(mu / rho), sqrt(mu / rho), 0, 0, 0
	};

	m->m[0].U = {
		1.0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0, 0, 0, 0,
		1.0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0, 0, 0, 0,
		0, 1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
		0, 0, 1.0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0,
		0, 1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
		0, 0, 1.0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0,
		0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 0, 1.0, 0, 0.0,
		0, 0, 0, 0, 0, 0, 0, 1.0, 0,
		0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 0, 0, 0, 1.0
	};

	m->m[0].U1 = {
		0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0.5, 0, 0.5, 0, 0, 0, 0,
		0, 0, 0, 0.5, 0, 0.5, 0, 0, 0,
		0.5 * sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
		0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
		(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
		-(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0,
		(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
		-(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1
	};


	m->m[1].A = {
		0, 0, 0, 0, -1.0 / rho, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, -1.0 / rho, 0, 0,
		0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
		0, -lambda, 0, 0, 0, 0, 0, 0, 0,
		-mu, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, -lambda - 2 * mu, 0, 0, 0, 0, 0, 0, 0,
		0, 0, -mu, 0, 0, 0, 0, 0, 0,
		0, -lambda, 0, 0, 0, 0, 0, 0, 0
	};

	m->m[1].L =
	{
		-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), -sqrt(
		        mu / rho), -sqrt(mu / rho),
		sqrt(mu / rho), sqrt(mu / rho), 0, 0, 0
	};

	m->m[1].U = {
		0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
		0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))), 0, 0,
		1.0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
		0, 0, 1.0, 0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
		1.0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0, 0,
		0, 0, 1.0, 0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
		0, 0, 0, 1.0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 0.0,
		0, 0, 0, 0, 0, 1.0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu), 0, 1.0
	};

	m->m[1].U1 = {
		0, 0, 0.5, 0, 0.5, 0, 0, 0, 0,
		0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0.5, 0, 0.5, 0, 0, 0,
		(0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), -(0.5 * lambda) / sqrt((lambda + 2 * mu) / rho), 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0,
		0.5 *sqrt(rho * (lambda + 2 * mu)), -0.5 * sqrt(rho * (lambda + 2 * mu)), 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0.5 * sqrt(mu * rho), 0, -0.5 * sqrt(mu * rho), 0, 0, 0,
		(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu),
		-(0.5 * sqrt(rho) * lambda) / sqrt(lambda + 2 * mu), 0, 0, 0, 0, 0, 0, 1
	};


	m->m[2].A = {
		0, 0, 0, 0, 0, -1.0 / rho, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, -1.0 / rho, 0,
		0, 0, 0, 0, 0, 0, 0, 0, -1.0 / rho,
		0, 0, -lambda, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0,
		-mu, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, -lambda, 0, 0, 0, 0, 0, 0,
		0, -mu, 0, 0, 0, 0, 0, 0, 0,
		0, 0, -lambda - 2 * mu, 0, 0, 0, 0, 0, 0
	};

	m->m[2].L =
	{
		-sqrt((lambda + 2 * mu) / rho), sqrt((lambda + 2 * mu) / rho), 
		-sqrt(mu / rho), -sqrt(mu / rho),
		sqrt(mu / rho), sqrt(mu / rho), 0, 0, 0
	};

	m->m[2].U = {
		0, 0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(rho * (lambda + 2 * mu))),
		0, 0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(rho * (lambda + 2 * mu))),
		1.0, 0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0, 0, 0,
		0, 1.0, 0, 0, 0, 0, 0, 1.0 / (sqrt(mu * rho)), 0,
		1.0, 0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0, 0, 0,
		0, 1.0, 0, 0, 0, 0, 0, -1.0 / (sqrt(mu * rho)), 0,
		0, 0, 0, 1.0, 0, 0, 0, 0, -(1.0 * lambda) / (lambda + 2 * mu),
		0, 0, 0, 0, 1.0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 1.0, 0, -(1.0 * lambda) / (lambda + 2 * mu)
	};

	m->m[2].U1 = {
		0, 0, 0.5, 0, 0.5, 0, 0, 0, 0,
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
		        rho * (lambda + 2 * mu)), 0, 0, 0, 0, 0, 0, 0
	};

	m->checkDecomposition();
}


template<typename ... Args>
void Elastic3DModel::
constructGcmMatrices(GcmMatricesPtr m, std::shared_ptr<const OrthotropicMaterial> material,
                     const Args& ... args) {

	if (material->anglesOfRotation == Real3::Zeros()) {
		constructNotRotated(m, material, args...);
		m->checkDecomposition();
	} else {
		constructRotated(m, material, args...);
	}
}


}

#endif // LIBGCM_MODEL_HPP
