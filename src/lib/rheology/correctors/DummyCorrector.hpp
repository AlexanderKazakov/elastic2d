#ifndef LIBGCM_DUMMYCORRECTOR_HPP
#define LIBGCM_DUMMYCORRECTOR_HPP

#include <lib/util/task/Task.hpp>

namespace gcm {
	/**
	 * Dummy (stub) corrector for the case when there aren't correctors
	 */
	class DummyCorrector {
	public:
		static const bool NonTrivial = false;

		void initialize(const Task&) { };

		template<typename TNode>
		void apply(TNode&) const { };
	};
};

#endif // LIBGCM_DUMMYCORRECTOR_HPP
