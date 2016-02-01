#ifndef LIBGCM_DUMMYODE_HPP
#define LIBGCM_DUMMYODE_HPP

#include <cmath>

#include <lib/rheology/variables/GetSetter.hpp>
#include <lib/util/Types.hpp>
#include <lib/util/task/Task.hpp>

namespace gcm {
	/**
	 * Dummy (stub) for cases when internal ODE is not present
	 */
	class DummyOde {
	public:
		static const bool NonTrivial = false;
		struct Variables { };

		static const std::map<PhysicalQuantities::T, GetSetter<Variables>> QUANTITIES;

		void initialize(const Task &task) {
			SUPPRESS_WUNUSED(task);
		};

		template<typename TNode>
		void nextStep(TNode& current, TNode& previous, const real timeStep) {
			SUPPRESS_WUNUSED(current);
			SUPPRESS_WUNUSED(previous);
			SUPPRESS_WUNUSED(timeStep);
		};
	};
};

#endif // LIBGCM_DUMMYODE_HPP
