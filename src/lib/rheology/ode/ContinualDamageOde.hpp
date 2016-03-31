#ifndef LIBGCM_CONTINUALDAMAGEODE_HPP
#define LIBGCM_CONTINUALDAMAGEODE_HPP

#include <cmath>

#include <lib/rheology/variables/GetSetter.hpp>
#include <lib/util/Types.hpp>
#include <lib/util/task/Task.hpp>

namespace gcm {
/**
 * The simplest continual damage model
 */
class ContinualDamageOde {
public:
	static const bool NonTrivial = true;

	struct Variables {
		real hi = 0;
	};

	static const std::map<PhysicalQuantities::T, GetSetter<Variables> > QUANTITIES;
	static real GetHi(const Variables& variablesToGetFrom) { return variablesToGetFrom.hi; }
	static void SetHi(const real& value, Variables& variablesToSetTo) { variablesToSetTo.hi =
		                                                                    value; }

	ContinualDamageOde(const Statement&) { }

	template<typename TNodePtr>
	void nextStep(TNodePtr node, const real timeStep) {
		node->_ode().hi += timeStep * dHiDt(node, timeStep);
	}

	/**
	 * This is stupid model just for demonstration
	 * @return time-derivative of damage measure
	 */
	template<typename TNodePtr>
	real dHiDt(const TNodePtr node, const real) {
		return fabs(node->pde().getPressure()) * node->material()->continualDamageParameter;
	}

};


}

#endif // LIBGCM_CONTINUALDAMAGEODE_HPP
