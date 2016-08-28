#ifndef LIBGCM_ODE_HPP
#define LIBGCM_ODE_HPP

#include <cmath>

#include <libgcm/grid/AbstractGrid.hpp>
#include <libgcm/util/task/Task.hpp>


namespace gcm {


/**
 * Class to perform correctors, ordinary evolution equations and other
 * actions not connected with waves propagation
 */
class AbstractOde {
public:
	virtual void apply(AbstractGrid& mesh_, const real timeStep) = 0;
};



/**
 * The simplest Maxwell viscosity model
 */
template<typename TMesh>
class MaxwellViscosityOde : public AbstractOde {
public:
	virtual void apply(AbstractGrid& mesh_, const real timeStep) override {
		TMesh& mesh = dynamic_cast<TMesh&>(mesh_);
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
class ContinualDamageOde : public AbstractOde {
public:
	virtual void apply(AbstractGrid& mesh_, const real timeStep) override {
		TMesh& mesh = dynamic_cast<TMesh&>(mesh_);
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
class IdealPlasticFlowCorrector : public AbstractOde {
	static const int D = TMesh::DIMENSIONALITY;
	typedef linal::DiagonalMatrix<D> DiagD;
public:
	virtual void apply(AbstractGrid& mesh_, const real) override {
		TMesh& mesh = dynamic_cast<TMesh&>(mesh_);
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


}

#endif // LIBGCM_ODE_HPP
