#ifndef LIBGCM_CONTINUALDAMAGEODE_HPP
#define LIBGCM_CONTINUALDAMAGEODE_HPP

#include <cmath>

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

		static const std::map<PhysicalQuantities::T, GetSetter<Variables>> QUANTITIES;
		static real GetHi(const Variables& variablesToGetFrom) { return variablesToGetFrom.hi; };
		static void SetHi(const real& value, Variables& variablesToSetTo) { variablesToSetTo.hi = value; };

		void initialize(const Task &task) {
			parameter = task.continualDamageParameter;
			assert_gt(parameter, 0.0);
		};

		/**
		 * Calculate new ODE values and place it to current.
		 * Given ODE values from current will be moved to previous.
		 * Given ODE values from previous will be discarded.
		 */
		template<typename PDEVariables>
		void nextStep(Variables& odeVariables, PDEVariables& pdeVariables, const real timeStep) {
			odeVariables.hi += timeStep * dHiDt(odeVariables, pdeVariables, timeStep);
		}

		/**
		 * This is stupid model just for demonstration
		 * @return time-derivative of damage measure
		 */
		template<typename PDEVariables>
		real dHiDt(Variables& odeVariables, PDEVariables& pdeVariables, const real timeStep) {
			SUPPRESS_WUNUSED(odeVariables);
			SUPPRESS_WUNUSED(timeStep);
			return fabs(pdeVariables.getPressure()) * parameter;
		}

	protected:
		real parameter = 0; // parameter in evolution equation

	};
};

#endif // LIBGCM_CONTINUALDAMAGEODE_HPP
