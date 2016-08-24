#ifndef LIBGCM_ODE_HPP
#define LIBGCM_ODE_HPP

#include <cmath>

#include <lib/rheology/variables/GetSetter.hpp>
#include <lib/util/Types.hpp>
#include <lib/util/task/Task.hpp>

namespace gcm {


/**
 * Class to perform correctors, ordinary evolution equations and other
 * actions not connected with waves propagation
 */
template<typename TMesh>
class AbstractOde {
public:
	virtual void apply(TMesh& mesh, const real timeStep) = 0;
};



/**
 * The simplest Maxwell viscosity model
 */
template<typename TMesh>
class MaxwellViscosityOde : public AbstractOde<TMesh> {
public:
	virtual void apply(TMesh& mesh, const real timeStep) override {
		for (auto& iter : mesh) {
			mesh._pdeVars(iter).setSigma(mesh.pdeVars(iter).getSigma() *
					exp(-timeStep / mesh.material(iter)->tau0));
		}
	}
};


/**
 * The simplest continual damage model
 */
template<typename TMesh>
class ContinualDamageOde : public AbstractOde<TMesh> {
public:
	virtual void apply(TMesh& mesh, const real timeStep) override {
		for (auto& iter : mesh) {
			const real dHiDt = fabs(mesh.pdeVars(iter).getPressure());
			mesh._ode(iter) += timeStep * dHiDt;
		}
	}
};


/**
 * The simplest plasticity flow model corrector
 */
template<typename TMesh>
class IdealPlasticFlowCorrector : public AbstractOde<TMesh> {
	static const int D = TMesh::DIMENSIONALITY;
	typedef linal::DiagonalMatrix<D> DiagD;
public:
	virtual void apply(TMesh& mesh, const real) override {
		for (auto& iter : mesh) {
			real pressure = mesh.pdeVars(iter).getPressure();
			real J2 = mesh.pdeVars(iter).getJ2();
			
			/// Correction parameter
			real x = J2 / mesh.material(iter)->yieldStrength;
			
			if (x > 1) {
				DiagD spherical = DiagD::Identity() * pressure;
				mesh._pdeVars.setSigma(
						mesh.pdeVars.getSigma() + (1.0 / x - 1) * spherical);
			}
		}
	}
};



template<typename TMesh>
class OdeFactory {
public:
	
	static std::shared_ptr<AbstractOde<TMesh>> create(const Odes::T type) {
		switch (type) {
			case Odes::T::MAXWELL_VISCOSITY:
				return std::make_shared<MaxwellViscosityOde<TMesh>>();
//			case Odes::T::CONTINUAL_DAMAGE:
//				return std::make_shared<ContinualDamageOde<TMesh>>();
//			case Odes::T::IDEAL_PLASTIC_FLOW:
//				return std::make_shared<IdealPlasticFlowCorrector<TMesh>>();
			default:
				THROW_UNSUPPORTED("Unknown or unsupported ODE type");
		}
	}
	
};


}

#endif // LIBGCM_ODE_HPP
