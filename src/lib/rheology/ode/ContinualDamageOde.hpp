#ifndef LIBGCM_CONTINUALDAMAGEODE_HPP
#define LIBGCM_CONTINUALDAMAGEODE_HPP

#include <type_traits>
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
			parameter = task.material.continualDamageParameter;
		};

		/**
		 * Calculate new ODE values and place it to current.
		 * Given ODE values from current will be moved to previous.
		 * Given ODE values from previous will be discarded.
		 */
		template<typename TNode>
		void nextStep(TNode& current, TNode& previous, const real timeStep) {
			static_assert(std::is_same<typename TNode::Model::OdeVariables, Variables>::value,
			              "ODE called for another's variables");

			real currentHi = current.hi;
			current.hi += timeStep * dHiDt(current, previous, timeStep);
			previous.hi = currentHi;
		}

		/**
		 * This is stupid model just for demonstration
		 * @return time-derivative of damage measure
		 */
		template<typename TNode>
		real dHiDt(TNode& current, TNode& previous, const real timeStep) {
			SUPPRESS_WUNUSED(timeStep);
			SUPPRESS_WUNUSED(previous);
			return fabs(current.u.getPressure()) * parameter;
		}

	protected:
		real parameter = 0; // parameter in evolution equation

	};
};

#endif // LIBGCM_CONTINUALDAMAGEODE_HPP
